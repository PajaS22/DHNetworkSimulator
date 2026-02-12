# Network nodes and edges

This page documents the core data types used to represent a district heating network.

The package models a network as a directed graph where:

- vertices are *nodes* (`ProducerNode`, `JunctionNode`, `LoadNode`), and
- edges are *pipes* (`InsulatedPipe`) carrying water plugs.

## Conventions

- Temperatures are in °C.
- Mass flow is in kg/s.
- Pipe lengths are in m, inner diameters in m.
- Power demand curves in `LoadNode` use coefficients in kW.

!!! note "`missing` means not initialized"
	Many fields use `missing` to represent “not set” or “not computed yet”. For example, `mass_flow` is typically filled during the
	steady-state hydraulics step.

## Nodes

```@docs
NodeType
NodeCommon
JunctionNode
LoadNode
ProducerNode
EmptyNode
```

## Edges
```@docs
EdgeType
EmptyEdge
Plug
InsulatedPipe
PipeParams
```