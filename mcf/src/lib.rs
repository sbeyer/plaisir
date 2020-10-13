use thiserror::Error;

#[derive(Error, Debug, PartialEq)]
pub enum MCFError {
    #[error("problem is not feasible")]
    Infeasible,
    #[error("internal error")]
    Internal(#[source] minilp::Error),
}

#[derive(Clone, Debug)]
pub struct EdgeValues {
    lower_bound: f64,
    upper_bound: f64,
    cost: f64,
}

impl EdgeValues {
    fn new(lb: f64, ub: f64, cost: f64) -> Self {
        Self {
            lower_bound: lb,
            upper_bound: ub,
            cost: cost,
        }
    }
}

impl Default for EdgeValues {
    fn default() -> Self {
        // another default is probably better, or none at all, but this is a start
        Self::new(0.0, f64::INFINITY, 1.0)
    }
}

type Instance = petgraph::Graph<f64, EdgeValues, petgraph::Directed>;

#[derive(Debug)]
pub struct Solution {
    flow: Vec<f64>,
    cost: f64,
}

pub fn run(instance: &Instance) -> Result<Solution, MCFError> {
    let mut lp = minilp::Problem::new(minilp::OptimizationDirection::Minimize);
    let mut vars = Vec::with_capacity(instance.edge_count());

    for edge in instance.edge_indices() {
        let val = &instance[edge];
        vars.push(lp.add_var(val.cost, (val.lower_bound, val.upper_bound)));
        // XXX: assumes indices of edges are 0..m-1
    }

    for node in instance.node_indices() {
        use petgraph::visit::EdgeRef;

        let mut lhs = minilp::LinearExpr::empty();

        for edge in instance.edges_directed(node, petgraph::EdgeDirection::Incoming) {
            let edge_index = edge.id().index();
            lhs.add(vars[edge_index], 1.0);
        }

        for edge in instance.edges_directed(node, petgraph::EdgeDirection::Outgoing) {
            let edge_index = edge.id().index();
            lhs.add(vars[edge_index], -1.0);
        }

        lp.add_constraint(lhs, minilp::ComparisonOp::Eq, -instance[node])
    }

    let result = lp.solve();

    match result {
        Ok(solution) => {
            let mut flow = Vec::<f64>::with_capacity(instance.edge_count());
            for var in vars {
                flow.push(solution[var]);
            }

            let cost = solution.objective();

            Ok(Solution {
                flow: flow,
                cost: cost,
            })
        }
        Err(minilp::Error::Infeasible) => Err(MCFError::Infeasible),
        Err(error) => Err(error).map_err(|source| MCFError::Internal(source)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_mcf_on_single_node() {
        let mut graph = Instance::new();
        graph.add_node(0.0);

        let solution = run(&graph).unwrap();

        assert_eq!(solution.cost, 0.0);
    }

    #[test]
    fn computes_mcf_on_easy_instance() {
        let mut graph = Instance::new();
        let s = graph.add_node(2.0);
        let u = graph.add_node(0.0);
        let v = graph.add_node(0.0);
        let t = graph.add_node(-2.0);
        let su = graph.add_edge(s, u, EdgeValues::new(0.0, 2.0, 11.0));
        let ut = graph.add_edge(u, t, EdgeValues::new(0.0, 2.0, 11.0));
        let st = graph.add_edge(s, t, EdgeValues::new(0.0, 2.0, 12.0));
        let sv = graph.add_edge(s, v, EdgeValues::new(0.0, 5.0, 2.0));
        let vt = graph.add_edge(v, t, EdgeValues::new(0.0, 1.0, 6.0));

        let solution = run(&graph).unwrap();

        assert_eq!(solution.cost, 20.0);
        assert_eq!(solution.flow[su.index()], 0.0);
        assert_eq!(solution.flow[ut.index()], 0.0);
        assert_eq!(solution.flow[st.index()], 1.0);
        assert_eq!(solution.flow[sv.index()], 1.0);
        assert_eq!(solution.flow[vt.index()], 1.0);
    }

    #[test]
    fn errors_on_infeasible_instance_due_to_disconnected_nodes() {
        let mut graph = Instance::new();
        graph.add_node(1.0);
        graph.add_node(-1.0);

        let result = run(&graph);

        assert!(result.is_err());
        assert_eq!(result.err(), Some(MCFError::Infeasible));
    }

    #[test]
    fn errors_on_infeasible_instance_due_to_imbalanced_supply_and_demand() {
        let mut graph = Instance::new();
        let s = graph.add_node(2.0);
        let t = graph.add_node(-1.0);
        graph.add_edge(s, t, EdgeValues::new(0.0, 2.0, 1.0));

        let result = run(&graph);

        assert!(result.is_err());
        assert_eq!(result.err(), Some(MCFError::Infeasible));
    }

    #[test]
    fn errors_on_infeasible_instance_due_to_capacity_constraints() {
        let mut graph = Instance::new();
        let s = graph.add_node(2.0);
        let t = graph.add_node(-2.0);
        graph.add_edge(s, t, EdgeValues::new(0.0, 1.0, 1.0));

        let result = run(&graph);

        assert!(result.is_err());
        assert_eq!(result.err(), Some(MCFError::Infeasible));
    }
}
