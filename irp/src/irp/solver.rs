use super::*;

type Node = petgraph::graph::NodeIndex;

struct AuxiliaryMCFInstance {
    // The underlying graph for the mcf algorithm
    graph: mcf::Instance,

    // Location nodes indexed by location and day
    customers: Vec<Vec<Node>>,

    // Depot source nodes indexed by day
    depot_source: Vec<Node>,

    // Depot target nodes indexed by day
    depot_target: Vec<Node>,
}

impl AuxiliaryMCFInstance {
    fn empty() -> Self {
        Self {
            graph: mcf::Instance::new(),
            customers: Vec::<Vec<Node>>::new(),
            depot_source: Vec::<Node>::new(),
            depot_target: Vec::<Node>::new(),
        }
    }

    fn new(problem: Problem) -> Self {
        let builder = AuxiliaryMCFInstanceBuilder::new(problem);
        builder.build()
    }
}

struct AuxiliaryMCFInstanceBuilder {
    problem: Problem,
    instance: AuxiliaryMCFInstance,
    flow_imbalance: f64,
}

impl AuxiliaryMCFInstanceBuilder {
    fn new(problem: Problem) -> Self {
        AuxiliaryMCFInstanceBuilder {
            problem: problem,
            instance: AuxiliaryMCFInstance::empty(),
            flow_imbalance: 0.0,
        }
    }

    fn build(mut self) -> AuxiliaryMCFInstance {
        self.add_location_nodes();
        self.add_free_edges_between_depot_and_customers();
        self.add_daily_production();
        self.add_daily_consumption();
        self.add_overnight_edges();
        self.add_super_sink();

        self.instance
    }

    fn new_node(&mut self, value: i32) -> Node {
        let value = value as f64;
        self.flow_imbalance += value;
        self.instance.graph.add_node(value)
    }

    fn new_location_nodes_vector(&mut self, start_level: i32) -> Vec<Node> {
        let mut v = Vec::<Node>::new();
        v.push(self.new_node(start_level));

        for _ in 0..self.problem.num_days {
            v.push(self.new_node(0));
        }

        debug_assert!(v.len() == self.problem.num_days + 1);

        v
    }

    fn add_location_nodes(&mut self) {
        self.instance.depot_source = self.new_location_nodes_vector(self.problem.depot.start_level);
        self.instance.depot_target = self.new_location_nodes_vector(0);
        // we don't need a depot target node for the last day, but keep it for simplicity

        for i in 0..(self.problem.num_nodes - 1) {
            let v = self.new_location_nodes_vector(self.problem.customers[i].start_level);
            self.instance.customers.push(v);
        }
    }

    // temporary solution to compute lower bound... TODO
    fn add_free_edges_between_depot_and_customers(&mut self) {
        for customer in 0..(self.problem.num_nodes - 1) {
            for _vehicle in 0..self.problem.num_vehicles {
                for day in 0..self.problem.num_days {
                    self.instance.graph.add_edge(
                        self.instance.depot_source[day],
                        self.instance.customers[customer][day],
                        mcf::FlowValues::new(0.0, self.problem.capacity.into(), 0.0),
                    );
                    self.instance.graph.add_edge(
                        self.instance.customers[customer][day],
                        self.instance.depot_target[day],
                        mcf::FlowValues::new(0.0, self.problem.capacity.into(), 0.0),
                    );
                }
            }
        }
    }

    fn add_daily_production(&mut self) {
        let daily_production_value = self.problem.depot.daily_production;

        for day in 1..self.problem.num_days {
            let node = self.new_node(daily_production_value);
            self.instance.graph.add_edge(
                node,
                self.instance.depot_source[day],
                mcf::FlowValues::new(
                    daily_production_value.into(),
                    daily_production_value.into(),
                    0.0,
                ),
            );
        }
    }

    fn add_daily_consumption(&mut self) {
        for customer in 0..(self.problem.num_nodes - 1) {
            let daily_consumption_value = self.problem.customers[customer].daily_consumption;

            for day in 0..self.problem.num_days {
                let node = self.new_node(-daily_consumption_value);
                self.instance.graph.add_edge(
                    self.instance.customers[customer][day],
                    node,
                    mcf::FlowValues::new(
                        daily_consumption_value.into(),
                        daily_consumption_value.into(),
                        0.0,
                    ),
                );
            }
        }
    }

    fn add_overnight_edges(&mut self) {
        // TODO
    }

    fn add_super_sink(&mut self) {
        // TODO
    }
}

struct BranchAndBound {
    //instance: AuxiliaryMCFInstance,
}

impl BranchAndBound {
    fn solve(problem: Problem) {
        let mut _instance = AuxiliaryMCFInstance::new(problem);

        // TODO
    }
}

pub fn solve(problem: Problem) {
    BranchAndBound::solve(problem)
}
