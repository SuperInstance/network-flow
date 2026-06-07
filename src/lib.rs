//! # network-flow
//!
//! Maximum flow (Ford-Fulkerson/Edmonds-Karp), minimum s-t cut, residual graphs,
//! and multi-commodity flow for network optimization problems.
//!
//! ## Overview
//!
//! Network flow problems are fundamental in combinatorial optimization. Given a directed
//! graph with edge capacities, the **maximum flow problem** asks: what is the maximum
//! amount of "flow" that can be sent from a source to a sink without exceeding edge capacities?
//!
//! The **min-cut max-flow theorem** states that the maximum flow equals the minimum cut capacity.
//!
//! ## Core Types
//!
//! - [`FlowNetwork`] — directed graph with edge capacities
//! - [`FordFulkerson`] — max flow via augmenting paths (BFS-based = Edmonds-Karp)
//! - [`MinCut`] — minimum s-t cut from max flow
//! - [`ResidualGraph`] — compute residual capacities after flow assignment
//! - [`MultiCommodityFlow`] — multiple source-sink pairs with demands

/// A flow network: directed graph with edge capacities.
#[derive(Clone, Debug)]
pub struct FlowNetwork {
    /// Number of nodes
    n: usize,
    /// Adjacency list: adjacency[u] = list of (v, capacity, flow)
    /// Stored as capacity matrix for simplicity
    capacity: Vec<Vec<f64>>,
}

impl FlowNetwork {
    /// Create a new flow network with `n` nodes (no edges).
    pub fn new(n: usize) -> Self {
        Self {
            n,
            capacity: vec![vec![0.0; n]; n],
        }
    }

    /// Add a directed edge from `u` to `v` with the given capacity.
    pub fn add_edge(&mut self, u: usize, v: usize, capacity: f64) {
        assert!(u < self.n && v < self.n, "Node index out of bounds");
        assert!(capacity >= 0.0, "Capacity must be non-negative");
        self.capacity[u][v] += capacity;
    }

    /// Add an undirected edge (both directions).
    pub fn add_undirected_edge(&mut self, u: usize, v: usize, capacity: f64) {
        self.add_edge(u, v, capacity);
        self.add_edge(v, u, capacity);
    }

    /// Number of nodes.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Returns true if the network has no nodes.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Get capacity of edge (u, v).
    pub fn capacity(&self, u: usize, v: usize) -> f64 {
        self.capacity[u][v]
    }

    /// Get the capacity matrix.
    pub fn capacity_matrix(&self) -> &[Vec<f64>] {
        &self.capacity
    }

    /// Get neighbors of node `u` (nodes with positive capacity from u).
    pub fn neighbors(&self, u: usize) -> Vec<usize> {
        (0..self.n)
            .filter(|&v| self.capacity[u][v] > 0.0)
            .collect()
    }
}

/// Max flow computation using the Edmonds-Karp algorithm (BFS-based Ford-Fulkerson).
///
/// Time complexity: O(V * E²) where V = nodes, E = edges.
pub struct FordFulkerson;

impl FordFulkerson {
    /// Compute maximum flow from `source` to `sink` in the given network.
    ///
    /// Returns the max flow value and the flow matrix.
    pub fn max_flow(network: &FlowNetwork, source: usize, sink: usize) -> (f64, Vec<Vec<f64>>) {
        let n = network.len();
        let mut residual = network.capacity_matrix().to_vec();
        let mut flow = vec![vec![0.0; n]; n];
        let mut total_flow = 0.0;

        loop {
            // BFS to find augmenting path
            let path = Self::bfs_path(&residual, source, sink, n);
            match path {
                Some(path) => {
                    // Find minimum residual capacity along path
                    let mut min_cap = f64::INFINITY;
                    for i in 0..path.len() - 1 {
                        let u = path[i];
                        let v = path[i + 1];
                        min_cap = min_cap.min(residual[u][v]);
                    }

                    // Augment flow
                    for i in 0..path.len() - 1 {
                        let u = path[i];
                        let v = path[i + 1];
                        flow[u][v] += min_cap;
                        flow[v][u] -= min_cap;
                        residual[u][v] -= min_cap;
                        residual[v][u] += min_cap;
                    }

                    total_flow += min_cap;
                }
                None => break,
            }
        }

        (total_flow, flow)
    }

    /// BFS to find a shortest augmenting path in the residual graph.
    fn bfs_path(residual: &[Vec<f64>], source: usize, sink: usize, n: usize) -> Option<Vec<usize>> {
        let mut parent = vec![None; n];
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();

        visited[source] = true;
        queue.push_back(source);

        while let Some(u) = queue.pop_front() {
            if u == sink {
                // Reconstruct path
                let mut path = Vec::new();
                let mut current = sink;
                path.push(current);
                while let Some(p) = parent[current] {
                    path.push(p);
                    current = p;
                }
                path.reverse();
                return Some(path);
            }

            for v in 0..n {
                if !visited[v] && residual[u][v] > 1e-10 {
                    visited[v] = true;
                    parent[v] = Some(u);
                    queue.push_back(v);
                }
            }
        }

        None
    }
}

/// Minimum s-t cut computation.
///
/// After computing max flow, the min cut is found by identifying all nodes
/// reachable from the source in the residual graph.
pub struct MinCut;

impl MinCut {
    /// Compute the minimum s-t cut.
    ///
    /// Returns (cut_value, reachable_nodes, unreachable_nodes, cut_edges).
    pub fn compute(network: &FlowNetwork, source: usize, sink: usize) -> (f64, Vec<usize>, Vec<usize>, Vec<(usize, usize)>) {
        let (max_flow_val, _) = FordFulkerson::max_flow(network, source, sink);

        // Build residual graph from flow
        let n = network.len();

        // Recompute residual from max_flow
        let (_, flow) = FordFulkerson::max_flow(network, source, sink);
        let mut residual = vec![vec![0.0; n]; n];
        for u in 0..n {
            for v in 0..n {
                residual[u][v] = (network.capacity(u, v) - flow[u][v]).max(0.0);
            }
        }

        // BFS from source in residual graph
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        visited[source] = true;
        queue.push_back(source);

        while let Some(u) = queue.pop_front() {
            for v in 0..n {
                if !visited[v] && residual[u][v] > 1e-10 {
                    visited[v] = true;
                    queue.push_back(v);
                }
            }
        }

        let reachable: Vec<usize> = (0..n).filter(|&i| visited[i]).collect();
        let unreachable: Vec<usize> = (0..n).filter(|&i| !visited[i]).collect();

        // Find cut edges (from reachable to unreachable)
        let mut cut_edges = Vec::new();
        for &u in &reachable {
            for &v in &unreachable {
                if network.capacity(u, v) > 0.0 {
                    cut_edges.push((u, v));
                }
            }
        }

        (max_flow_val, reachable, unreachable, cut_edges)
    }
}

/// Residual graph computation.
pub struct ResidualGraph;

impl ResidualGraph {
    /// Compute the residual graph given a flow network and a flow assignment.
    ///
    /// Residual capacity: r(u,v) = c(u,v) - f(u,v)
    pub fn compute(network: &FlowNetwork, flow: &[Vec<f64>]) -> Vec<Vec<f64>> {
        let n = network.len();
        let mut residual = vec![vec![0.0; n]; n];
        for u in 0..n {
            for v in 0..n {
                let cap = network.capacity(u, v);
                let f = flow[u][v];
                // Forward edge residual
                residual[u][v] = (cap - f).max(0.0);
                // Backward edge (already handled by flow matrix having negative entries)
            }
        }
        residual
    }

    /// Check if the given flow is valid: capacity constraints and flow conservation.
    pub fn is_valid_flow(network: &FlowNetwork, flow: &[Vec<f64>], source: usize, sink: usize) -> bool {
        let n = network.len();

        // Check capacity constraints (forward flow only)
        for u in 0..n {
            for v in 0..n {
                let f = flow[u][v];
                // Net flow from u to v should not exceed capacity
                if f > network.capacity(u, v) + 1e-10 {
                    return false;
                }
            }
        }

        // Check flow conservation (except source and sink)
        for v in 0..n {
            if v == source || v == sink {
                continue;
            }
            let inflow: f64 = (0..n).map(|u| flow[u][v].max(0.0)).sum();
            let outflow: f64 = (0..n).map(|u| flow[v][u].max(0.0)).sum();
            if (inflow - outflow).abs() > 1e-6 {
                return false;
            }
        }

        true
    }
}

/// Multi-commodity flow: multiple source-sink pairs each with a demand.
///
/// Uses a greedy sequential approach: route each commodity one at a time
/// through the remaining capacity.
pub struct MultiCommodityFlow {
    network: FlowNetwork,
}

impl MultiCommodityFlow {
    /// Create a new multi-commodity flow solver for the given network.
    pub fn new(network: FlowNetwork) -> Self {
        Self { network }
    }

    /// Attempt to route all commodities through the network.
    ///
    /// Returns a vector of (flow_value, flow_matrix) for each commodity,
    /// and the total flow routed.
    pub fn route(&self, commodities: &[(usize, usize, f64)]) -> (Vec<(f64, Vec<Vec<f64>>)>, f64) {
        let n = self.network.len();
        let mut remaining_capacity = self.network.capacity_matrix().to_vec();
        let mut results = Vec::new();
        let mut total = 0.0;

        for &(source, sink, demand) in commodities {
            // Create a temporary network with remaining capacity
            let mut temp = FlowNetwork::new(n);
            for u in 0..n {
                for v in 0..n {
                    if remaining_capacity[u][v] > 0.0 {
                        temp.add_edge(u, v, remaining_capacity[u][v]);
                    }
                }
            }

            // Compute max flow for this commodity
            let (flow_val, flow) = FordFulkerson::max_flow(&temp, source, sink);

            // Route at most the demand
            let routed = flow_val.min(demand);

            // Scale flow to match routed amount
            let scale = if flow_val > 1e-10 { routed / flow_val } else { 0.0 };
            let mut scaled_flow = vec![vec![0.0; n]; n];
            for u in 0..n {
                for v in 0..n {
                    scaled_flow[u][v] = flow[u][v] * scale;
                    remaining_capacity[u][v] -= scaled_flow[u][v];
                }
            }

            results.push((routed, scaled_flow));
            total += routed;
        }

        (results, total)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_max_flow() {
        // 0 --10--> 1 --10--> 2
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        net.add_edge(1, 2, 10.0);
        let (flow_val, _) = FordFulkerson::max_flow(&net, 0, 2);
        assert!((flow_val - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_bottleneck_flow() {
        // 0 --10--> 1 --5--> 2
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        net.add_edge(1, 2, 5.0);
        let (flow_val, _) = FordFulkerson::max_flow(&net, 0, 2);
        assert!((flow_val - 5.0).abs() < 1e-6);
    }

    #[test]
    fn test_parallel_paths() {
        // 0 --5--> 1 --5--> 3
        // 0 --5--> 2 --5--> 3
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 5.0);
        net.add_edge(0, 2, 5.0);
        net.add_edge(1, 3, 5.0);
        net.add_edge(2, 3, 5.0);
        let (flow_val, _) = FordFulkerson::max_flow(&net, 0, 3);
        assert!((flow_val - 10.0).abs() < 1e-6);
    }

    #[test]
    fn test_no_path() {
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        // No edge 1->2
        let (flow_val, _) = FordFulkerson::max_flow(&net, 0, 2);
        assert!((flow_val - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_min_cut() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 5.0);
        net.add_edge(0, 2, 5.0);
        net.add_edge(1, 3, 5.0);
        net.add_edge(2, 3, 5.0);
        let (cut_val, reachable, _, cut_edges) = MinCut::compute(&net, 0, 3);
        assert!((cut_val - 10.0).abs() < 1e-6);
        assert!(reachable.contains(&0));
    }

    #[test]
    fn test_residual_graph() {
        let mut net = FlowNetwork::new(2);
        net.add_edge(0, 1, 10.0);
        let flow = vec![vec![7.0, 7.0], vec![-7.0, -7.0]];
        let residual = ResidualGraph::compute(&net, &flow);
        assert!((residual[0][1] - 3.0).abs() < 1e-6);
    }

    #[test]
    fn test_valid_flow() {
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 10.0);
        net.add_edge(1, 2, 10.0);
        let (_, flow) = FordFulkerson::max_flow(&net, 0, 2);
        assert!(ResidualGraph::is_valid_flow(&net, &flow, 0, 2));
    }

    #[test]
    fn test_network_neighbors() {
        let mut net = FlowNetwork::new(3);
        net.add_edge(0, 1, 5.0);
        net.add_edge(0, 2, 3.0);
        let neighbors = net.neighbors(0);
        assert_eq!(neighbors.len(), 2);
        assert!(neighbors.contains(&1));
        assert!(neighbors.contains(&2));
    }

    #[test]
    fn test_network_capacity() {
        let mut net = FlowNetwork::new(2);
        net.add_edge(0, 1, 7.5);
        assert!((net.capacity(0, 1) - 7.5).abs() < 1e-10);
    }

    #[test]
    fn test_multi_commodity() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 1, 5.0);
        net.add_edge(0, 2, 5.0);
        net.add_edge(1, 3, 5.0);
        net.add_edge(2, 3, 5.0);
        let mcf = MultiCommodityFlow::new(net);
        let commodities = vec![(0, 3, 10.0)];
        let (results, total) = mcf.route(&commodities);
        assert!((total - 10.0).abs() < 1e-6);
        assert_eq!(results.len(), 1);
    }

    #[test]
    fn test_multi_commodity_competing() {
        let mut net = FlowNetwork::new(4);
        net.add_edge(0, 2, 5.0);
        net.add_edge(1, 2, 5.0);
        net.add_edge(2, 3, 6.0); // bottleneck
        let mcf = MultiCommodityFlow::new(net);
        let commodities = vec![(0, 3, 5.0), (1, 3, 5.0)];
        let (_, total) = mcf.route(&commodities);
        // First gets 5, second gets 1 (bottleneck 6 - 5 = 1)
        assert!((total - 6.0).abs() < 1e-6);
    }

    #[test]
    fn test_empty_network() {
        let net = FlowNetwork::new(0);
        assert!(net.is_empty());
    }

    #[test]
    fn test_add_edge_accumulates() {
        let mut net = FlowNetwork::new(2);
        net.add_edge(0, 1, 3.0);
        net.add_edge(0, 1, 4.0);
        assert!((net.capacity(0, 1) - 7.0).abs() < 1e-10);
    }
}
