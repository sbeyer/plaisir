use super::*;

type Node = petgraph::graph::NodeIndex;

struct AuxiliaryMCFInstance {
    // The underlying graph for the mcf algorithm
    graph: mcf::Instance,

    // Location nodes indexed by location and day
    v: Vec<Vec<Node>>,

    // Nodes for daily production and consumption
    daily: Vec<Vec<Node>>,
}

impl AuxiliaryMCFInstance {
    fn new() -> Self {
        let graph = mcf::Instance::new();
        let v = Vec::<Vec<Node>>::new();
        let daily = Vec::<Vec<Node>>::new();

        Self {
            graph: graph,
            v: v,
            daily: daily,
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

        // Add free edges between depot and customers... temporary solution to compute lower bound
        for customer in 1..problem.num_nodes {
            for _vehicle in 0..problem.num_vehicles {
                for day in 0..problem.num_days {
                    self.graph.add_edge(
                        self.depot(day),
                        self.customer(day, customer),
                        mcf::FlowValues::new(0.0, problem.capacity.into(), 0.0),
                    );
                }
                for day in 0..(problem.num_days - 1) {
                    self.graph.add_edge(
                        self.customer(day, customer),
                        self.depot(day + 1),
                        mcf::FlowValues::new(0.0, problem.capacity.into(), 0.0),
                    );
                }
            }
        }

        // Add daily production
        let mut daily_production = Vec::<Node>::new();
        let daily_production_value = problem.depot.daily_production as f64;

        daily_production.push(self.graph.add_node(0.0)); // dummy node, only necessary for indexing

        for day in 1..problem.num_days {
            flow_imbalance += daily_production_value;
            let node = self.graph.add_node(daily_production_value);
            daily_production.push(node);
            self.graph.add_edge(
                node,
                self.depot(day),
                mcf::FlowValues::new(daily_production_value, daily_production_value, 0.0),
            );
        }
        self.daily.push(daily_production);

        // TODO: daily consumption
        // TODO: super-sink with flow_imbalance
        // TODO
    }

    pub fn depot(&self, day: usize) -> Node {
        self.v[0][day]
    }

    pub fn customer(&self, day: usize, idx: usize) -> Node {
        assert!(idx > 0);
        self.v[idx][day]
    }

    fn init_location_nodes_vector(&mut self, start_level: i32, num_days: usize) {
        let mut v = Vec::<Node>::new();
        v.push(self.graph.add_node(start_level as f64));

        for _ in 0..num_days {
            v.push(self.graph.add_node(0.0));
        }

        self.v.push(v);
    }
}

struct BranchAndBound {
    //instance: AuxiliaryMCFInstance,
}

impl BranchAndBound {
    fn solve(problem: Problem) {
        let mut instance = AuxiliaryMCFInstance::new();
        instance.load(problem);
    }
}

pub fn solve(problem: Problem) {
    let solution = BranchAndBound::solve(problem);
}
