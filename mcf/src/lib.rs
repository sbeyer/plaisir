#[cfg(test)]
mod tests {
    #[test]
    fn petgraph_works() {
        let mut graph = petgraph::Graph::<i32, i32, petgraph::Directed>::new();
        let s = graph.add_node(2);
        let v = graph.add_node(0);
        let t = graph.add_node(-2);
        graph.extend_with_edges(&[(s, v), (v, t)]);
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 2);
    }
}
