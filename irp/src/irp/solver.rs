use super::*;

type NodesVec = Vec<petgraph::graph::NodeIndex>;

struct BranchAndBound {
    // The underlying graph for the mcf algorithm
    graph: mcf::Instance,

    // Location nodes indexed by location and day
    v: Vec<NodesVec>,
}

impl BranchAndBound {
    fn solve(problem: Problem) {
        let mut instance = Self::new();
        instance.load(problem);
    }

    fn new() -> Self {
        let mut graph = mcf::Instance::new();
        let mut v = Vec::<NodesVec>::new();

        Self {
            graph: graph,
            v: v,
        }
    }

    fn load(&mut self, problem: Problem) {
        let mut flow_imbalance = 0.0;

        self.init_location_nodes_vector(problem.depot.start_level, problem.num_days);
        flow_imbalance += problem.depot.start_level as f64;

        for i in 0..(problem.num_nodes - 1) {
            assert_eq!(problem.customers[i].id, i + 1);
            self.init_location_nodes_vector(problem.customers[i].start_level, problem.num_days);
            flow_imbalance += problem.customers[i].start_level as f64;
        }

        // TODO
    }

    fn init_location_nodes_vector(&mut self, start_level: i32, num_days: usize) {
        let mut v = NodesVec::new();
        v.push(self.graph.add_node(start_level as f64));

        for _ in 0..num_days {
            v.push(self.graph.add_node(0.0));
        }

        self.v.push(v);
    }
}

pub fn solve(problem: Problem) {
    let solution = BranchAndBound::solve(problem);
}
