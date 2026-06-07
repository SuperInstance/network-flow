# network-flow

**Maximum flow (Edmonds-Karp), minimum s-t cut, residual graphs, and multi-commodity flow for network optimization.**

[![Rust](https://img.shields.io/badge/rust-1.70%2B-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

## Overview

Network flow problems are among the most fundamental problems in combinatorial optimization and operations research. Given a directed graph with edge capacities, the **maximum flow problem** asks: what is the maximum amount of flow that can be sent from a source to a sink while respecting capacity constraints?

The **max-flow min-cut theorem** (a cornerstone of graph theory) establishes that the maximum flow equals the minimum cut capacity — the smallest total capacity of edges whose removal disconnects source from sink.

Applications include:

- **Transportation networks** — maximum throughput of goods through a supply chain
- **Bipartite matching** — maximum matching via reduction to max flow
- **Image segmentation** — min cut defines optimal foreground/background partition
- **Project selection** — optimal subset of projects under constraints
- **Network reliability** — connectivity analysis under edge failures

## Features

- **`FlowNetwork`** — Directed graph with edge capacities, supporting parallel edges
- **`FordFulkerson`** — Max flow via BFS augmenting paths (Edmonds-Karp: O(VE²))
- **`MinCut`** — Minimum s-t cut from max flow residual graph
- **`ResidualGraph`** — Compute and validate residual capacities
- **`MultiCommodityFlow`** — Multiple source-sink commodity pairs with demands

## Installation

```toml
[dependencies]
network-flow = "0.1.0"
```

## Quick Start

```rust
use network_flow::*;

// Build a simple network
let mut net = FlowNetwork::new(4);
net.add_edge(0, 1, 5.0); // source -> A, capacity 5
net.add_edge(0, 2, 5.0); // source -> B, capacity 5
net.add_edge(1, 3, 5.0); // A -> sink, capacity 5
net.add_edge(2, 3, 5.0); // B -> sink, capacity 5

// Compute max flow
let (max_flow, flow_matrix) = FordFulkerson::max_flow(&net, 0, 3);
assert_eq!(max_flow, 10.0);

// Compute min cut
let (cut_value, reachable, unreachable, cut_edges) = MinCut::compute(&net, 0, 3);
assert_eq!(cut_value, 10.0);
```

## Flow Network Construction

```rust
use network_flow::*;

let mut net = FlowNetwork::new(5);

// Add directed edges
net.add_edge(0, 1, 10.0);
net.add_edge(0, 2, 8.0);
net.add_edge(1, 2, 5.0);
net.add_edge(1, 3, 7.0);
net.add_edge(2, 4, 10.0);
net.add_edge(3, 4, 10.0);

// Query capacity
println!("Capacity 0→1: {}", net.capacity(0, 1));

// Get neighbors
println!("Neighbors of 0: {:?}", net.neighbors(0));

// Parallel edges accumulate capacity
net.add_edge(0, 1, 3.0); // now capacity 0→1 = 13.0
```

## Min Cut

The minimum cut partitions nodes into two sets: those reachable from source in the residual graph, and those not. The cut edges go from the reachable set to the unreachable set.

```rust
use network_flow::*;

let mut net = FlowNetwork::new(6);
net.add_edge(0, 1, 16.0);
net.add_edge(0, 2, 13.0);
net.add_edge(1, 2, 10.0);
net.add_edge(1, 3, 12.0);
net.add_edge(2, 1, 4.0);
net.add_edge(2, 4, 14.0);
net.add_edge(3, 2, 9.0);
net.add_edge(3, 5, 20.0);
net.add_edge(4, 3, 7.0);
net.add_edge(4, 5, 4.0);

let (cut_val, reachable, unreachable, cut_edges) = MinCut::compute(&net, 0, 5);
println!("Min cut value: {}", cut_val);
println!("Reachable: {:?}", reachable);
println!("Cut edges: {:?}", cut_edges);
```

## Residual Graph

```rust
use network_flow::*;

let mut net = FlowNetwork::new(2);
net.add_edge(0, 1, 10.0);

let (_, flow) = FordFulkerson::max_flow(&net, 0, 1);
let residual = ResidualGraph::compute(&net, &flow);

// Check flow validity
assert!(ResidualGraph::is_valid_flow(&net, &flow, 0, 1));
```

## Multi-Commodity Flow

Route multiple commodities through a shared network with limited capacity:

```rust
use network_flow::*;

let mut net = FlowNetwork::new(4);
net.add_edge(0, 2, 5.0);
net.add_edge(1, 2, 5.0);
net.add_edge(2, 3, 6.0); // bottleneck

let mcf = MultiCommodityFlow::new(net);
let commodities = vec![
    (0, 3, 5.0), // commodity 1: 0→3, demand 5
    (1, 3, 5.0), // commodity 2: 1→3, demand 5
];

let (results, total_flow) = mcf.route(&commodities);
println!("Total flow routed: {:.1} / {:.1}", total_flow, 10.0);
// First commodity gets 5, second gets 1 (bottleneck: 6)
```

## Algorithm Details

### Edmonds-Karp (BFS Ford-Fulkerson)

The algorithm repeatedly finds shortest augmenting paths (by BFS) in the residual graph and pushes maximum possible flow through them.

1. Initialize flow to 0 on all edges
2. While an augmenting path exists (BFS in residual graph):
   - Find the bottleneck capacity along the path
   - Augment: increase forward flow, decrease backward flow
   - Update residual capacities
3. Return total flow

**Time complexity**: O(V · E²) — guaranteed by BFS finding shortest paths.

### Min Cut

After max flow is computed:
1. Build the residual graph
2. BFS from source in the residual graph
3. Reachable nodes form the S-side of the cut
4. Unreachable nodes form the T-side
5. Cut edges: edges from S-side to T-side with positive original capacity

### Multi-Commodity Flow

Uses a greedy sequential approach:
1. Route commodities one at a time
2. After routing each commodity, reduce remaining capacity
3. Later commodities may get less flow due to earlier routing

This is a heuristic; optimal multi-commodity flow is NP-hard in general.

## API Reference

| Type | Key Methods | Description |
|------|-------------|-------------|
| `FlowNetwork` | `new`, `add_edge`, `capacity`, `neighbors` | Network construction |
| `FordFulkerson` | `max_flow` | Max flow + flow matrix |
| `MinCut` | `compute` | Min cut value, partition, cut edges |
| `ResidualGraph` | `compute`, `is_valid_flow` | Residual capacity analysis |
| `MultiCommodityFlow` | `new`, `route` | Multi-commodity routing |

## Theoretical Background

### Max-Flow Min-Cut Theorem
For any flow network with source s and sink t:
```
max flow = min cut capacity
```

This is one of the most elegant results in combinatorial optimization. The proof shows that:
1. max flow ≤ min cut (every flow must cross the cut)
2. When no augmenting path exists, the cut found equals the flow

### Integrality
If all capacities are integers, the maximum flow is also integer-valued. This enables applications to combinatorial problems (matching, disjoint paths).

## Performance

| Operation | Complexity |
|-----------|-----------|
| Max flow (Edmonds-Karp) | O(V · E²) |
| Min cut | O(V · E²) (dominated by max flow) |
| Residual graph | O(V²) |
| Multi-commodity | O(k · V · E²) where k = commodities |

## License

MIT License. See [LICENSE](LICENSE) for details.
