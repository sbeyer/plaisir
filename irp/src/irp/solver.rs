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
    fn empty() -> Self {
        let graph = mcf::Instance::new();
        let v = Vec::<Vec<Node>>::new();
        let daily = Vec::<Vec<Node>>::new();

        Self {
            graph: graph,
            v: v,
            daily: daily,
        }
    }

    fn new(problem: Problem) -> Self {
        let builder = AuxiliaryMCFInstanceBuilder::new(problem);
        builder.build()
    }

    pub fn depot(&self, day: usize) -> Node {
        self.v[0][day]
    }

    pub fn customer(&self, day: usize, idx: usize) -> Node {
        debug_assert!(idx > 0, "the depot (index 0) is not a customer");
        self.v[idx][day]
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

        // TODO: super-sink with flow_imbalance
        // TODO

        self.instance
    }

    fn new_node(&mut self, value: i32) -> Node {
        let value = value as f64;
        self.flow_imbalance += value;
        self.instance.graph.add_node(value)
    }

    fn init_location_nodes_vector(&mut self, start_level: i32) {
        let mut v = Vec::<Node>::new();
        v.push(self.new_node(start_level));

        for _ in 0..self.problem.num_days {
            v.push(self.new_node(0));
        }

        self.instance.v.push(v);
    }

    fn add_location_nodes(&mut self) {
        self.init_location_nodes_vector(self.problem.depot.start_level);

        for i in 0..(self.problem.num_nodes - 1) {
            debug_assert_eq!(self.problem.customers[i].id, i + 1);
            self.init_location_nodes_vector(self.problem.customers[i].start_level);
        }
    }

    // temporary solution to compute lower bound... TODO
    fn add_free_edges_between_depot_and_customers(&mut self) {
        for customer in 1..self.problem.num_nodes {
            for _vehicle in 0..self.problem.num_vehicles {
                for day in 0..self.problem.num_days {
                    self.instance.graph.add_edge(
                        self.instance.depot(day),
                        self.instance.customer(day, customer),
                        mcf::FlowValues::new(0.0, self.problem.capacity.into(), 0.0),
                    );
                }
                for day in 0..(self.problem.num_days - 1) {
                    self.instance.graph.add_edge(
                        self.instance.customer(day, customer),
                        self.instance.depot(day + 1),
                        mcf::FlowValues::new(0.0, self.problem.capacity.into(), 0.0),
                    );
                }
            }
        }
    }

    fn add_daily_production(&mut self) {
        let mut daily_production = Vec::<Node>::new();
        let daily_production_value = self.problem.depot.daily_production;

        daily_production.push(self.new_node(0)); // dummy node, only necessary for indexing

        for day in 1..self.problem.num_days {
            let node = self.new_node(daily_production_value);
            daily_production.push(node);
            self.instance.graph.add_edge(
                node,
                self.instance.depot(day),
                mcf::FlowValues::new(
                    daily_production_value.into(),
                    daily_production_value.into(),
                    0.0,
                ),
            );
        }
        self.instance.daily.push(daily_production);
    }

    fn add_daily_consumption(&mut self) {
        for customer in 1..self.problem.num_nodes {
            let mut daily_consumption = Vec::<Node>::new();
            let daily_consumption_value = self.problem.customers[customer - 1].daily_consumption;

            for day in 0..self.problem.num_days {
                let node = self.new_node(-daily_consumption_value);
                daily_consumption.push(node);
                self.instance.graph.add_edge(
                    self.instance.customer(day, customer),
                    node,
                    mcf::FlowValues::new(
                        daily_consumption_value.into(),
                        daily_consumption_value.into(),
                        0.0,
                    ),
                );
            }
            self.instance.daily.push(daily_consumption);
        }
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
