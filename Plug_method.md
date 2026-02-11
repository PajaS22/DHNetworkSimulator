
# Plug method for thermal dynamics

This project solves heat transport in a district heating network using a **plug-flow (parcel) method**.
Each pipe contains a queue of discrete **plugs of water**, where every plug has:

- temperature `T` [°C]
- mass `m` [kg]

Plugs advect through pipes according to the current mass flows, exchange heat with the environment via a simple heat-loss model, and (optionally) lose heat at loads according to the load power demand.

The implementation lives primarily in:

- `src/types.jl` (type `Plug` and pipe storage `plugs_f` / `plugs_b`)
- `src/simulation.jl` (`time_step_thermal_dynamics_forward!`, `time_step_thermal_dynamics_backward!`)

## Why plugs?

The approach is based on a **quasi-dynamic** assumption:

- **Hydraulics** (mass flow distribution) is assumed to reach steady state “instantaneously” compared to
- **Thermal dynamics**, which are dominated by advection (transport of hot water) and slow heat losses.

As a result, each simulation time step does:

1. compute steady-state mass flows (`steady_state_hydronynamics!`)
2. transport heat by moving plugs
   1. forward (supply) direction: from producer to loads
   2. heat extraction at loads: using a power-demand model (often based on outdoor-temperature compensation, sometimes called *equithermal regulation*)
   3. backward (return) direction: from loads back to the producer

## State representation

Each pipe edge represents two pipes and therefore stores two independent plug queues:

- `plugs_f`: plugs moving in the **forward/supply** direction
- `plugs_b`: plugs moving in the **backward/return** direction

Within a pipe, plugs are treated as **non-mixing** parcels (no axial mixing). Mixing is handled explicitly only where the network topology merges/splits flow.

## One simulation time step (high level)

For each time step of length $\Delta t$:

1. **Insert new hot plugs at the producer** into each outgoing supply pipe.
2. **Forward pass (supply)**: move plugs from producer to loads using current mass flows.
3. **Loads**: compute required power $P(T_a)$ from ambient temperature and cool the arriving plug accordingly.
4. **Backward pass (return)**: push cooled plugs back to the producer, mixing at junctions as needed.
5. **Heat losses to ambient**: cool all plugs remaining in pipes (supply and return) using an exponential heat-loss model.

The overall simulation loop is orchestrated by `run_simulation`.

---

## Forward pass (supply): advection and splitting

### 1) Plug injection at the source

For each outgoing edge from the producer, a new plug is created with:

$$ m_{in} = \dot m\,\Delta t $$

and temperature equal to the producer outlet temperature for that step. The plug is appended to that pipe’s forward queue.

### 2) Plug advection through a pipe

For a pipe with mass flow $\dot m$, the algorithm computes the mass that must exit the pipe in this step:

$$ m_{out} = \dot m\,\Delta t $$

It then pops plugs from the **front** of the queue until the exiting mass is reached. If a plug is larger than the remaining mass to exit, it is **split** into an exiting part and a remaining part.

This is implemented by `collect_exiting_water_plugs!`.

### 3) Junction splitting (tree)

When a node has multiple children, the exiting plug mass is split among outgoing edges proportional to edge mass flows:

$$ m_{child} = m_{plug}\,\frac{\dot m_{child}}{\sum_k \dot m_k} $$

Each child receives a new plug with the same temperature and the computed mass. These plugs are appended to the child edges’ forward queues.

### 4) Leaf handling (load inlet plug)

At a leaf node, all plugs that arrive during the step are combined into a single representative plug using a **mass-weighted average**:

$$ T_{avg} = \frac{\sum_i T_i m_i}{\sum_i m_i} $$

This is implemented by `combine_plugs`.

The resulting plug is interpreted as the **supply plug entering the load** for this time step.

---

## Load model: consuming heat

Each load node computes power demand based on ambient temperature $T_a$ (we use polynomial of 2nd order approximation)

$$ P = P(T_a) = p_1 + p_2 T_a + p_3 T_a^2 $$

The entering plug is cooled by energy extraction over the step:

$$ \Delta T = \frac{P\,\Delta t}{m\,c_p} $$

so the return-side plug temperature becomes $T - \Delta T$.

To avoid unphysical results (like cooling the plug to lower temperature than is inside the building), the implementation clamps return temperature to a configured minimal value (`MINIMAL_RETURN_TEMPERATURE = 25.0`).

---

## Backward pass (return): advection and merging

The backward pass pushes the cooled plugs from loads back to the producer using the return queues `plugs_b`.

### Leaf injection

At each load (leaf), the cooled plug is pushed into the parent edge’s return queue.

### Junction merging

At an internal node, return plugs are collected from each child return edge (again using `collect_exiting_water_plugs!`).
These multiple plug sequences must be merged into a single sequence going into the parent.

The code uses `merge_water_plug_vectors!`, which mixes the plugs according to the individual flows in merging pipes.

At the producer root, the merged plug sequence is combined to a single plug representing the **return temperature entering the producer** for that step.

---

## Heat loss to ambient (dissipation to atmosphere)

After advection, both in forward and backward pass, the remaining plugs in pipes are cooled to account for heat loss to the environment.

For each plug, temperature is updated using an exponential model:

$$
T_{next} = T_a + (T - T_a)\exp\left(-\frac{\Delta t}{\rho c_p A R}\right)
$$

where:

- $\rho$ is water density (`WATER_DENSITY`)
- $c_p$ is specific heat (`WATER_SPECIFIC_HEAT`)
- $A$ is cross-sectional area of the pipe ($\pi(d/2)^2$)
- $R$ is pipe thermal resistance ... this is given by pipe insulation and usually differs for forward and backward pass


---

## Practical notes and limitations

- The current network traversal assumes a **directed tree** (each non-root node has one parent).
- There is **no axial mixing** inside a pipe: plugs only merge when explicitly combined (e.g., at reporting points or junction return merging).
- The plug representation is simplified over time by merging consecutive plugs with nearly identical temperature (`merge_same_temperature_plugs!`).
- Stability and realism depend on choosing a reasonable $\Delta t$ relative to flows and pipe volumes.

