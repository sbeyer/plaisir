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

    #[test]
    fn minilp_works() {
        let mut problem = minilp::Problem::new(minilp::OptimizationDirection::Maximize);
        let x = problem.add_var(1.0, (0.0, f64::INFINITY));
        let y = problem.add_var(2.0, (0.0, 3.0));

        problem.add_constraint(&[(x, 1.0), (y, 1.0)], minilp::ComparisonOp::Le, 4.0);
        problem.add_constraint(&[(x, 2.0), (y, 1.0)], minilp::ComparisonOp::Ge, 2.0);

        let solution = problem.solve().unwrap();
        assert_eq!(solution.objective(), 7.0);
        assert_eq!(solution[x], 1.0);
        assert_eq!(solution[y], 3.0);
    }
}
