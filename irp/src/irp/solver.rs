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
        self.add_edges_between_depot_and_customers();
        self.add_routes_between_customers();
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

        for i in 0..self.problem.num_customers {
            let v = self.new_location_nodes_vector(self.problem.customers[i].start_level);
            self.instance.customers.push(v);
        }
    }

    fn add_edges_between_depot_and_customers(&mut self) {
        let capacity = self.problem.capacity as f64;
        for customer in 0..self.problem.num_customers {
            for _vehicle in 0..self.problem.num_vehicles {
                for day in 0..self.problem.num_days {
                    let distance = self
                        .problem
                        .depot
                        .position
                        .distance(&self.problem.customers[customer].position)
                        as f64;
                    let initial_cost = distance / capacity;
                    self.instance.graph.add_edge(
                        self.instance.depot_source[day],
                        self.instance.customers[customer][day],
                        mcf::FlowValues::new(0.0, capacity, initial_cost),
                    );
                    self.instance.graph.add_edge(
                        self.instance.customers[customer][day],
                        self.instance.depot_target[day],
                        mcf::FlowValues::new(0.0, capacity, initial_cost),
                    );
                }
            }
        }
    }

    fn add_routes_between_customers(&mut self) {
        let capacity = self.problem.capacity as f64;
        for source in 0..self.problem.num_customers {
            for target in 0..self.problem.num_customers {
                if source != target {
                    let distance = self.problem.customers[source]
                        .position
                        .distance(&self.problem.customers[target].position)
                        as f64;
                    let initial_cost = distance / capacity;
                    for day in 0..self.problem.num_days {
                        for _vehicle in 0..self.problem.num_vehicles {
                            self.instance.graph.add_edge(
                                self.instance.customers[source][day],
                                self.instance.customers[target][day],
                                mcf::FlowValues::new(0.0, capacity, initial_cost),
                            );
                        }
                    }
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
        for customer in 0..self.problem.num_customers {
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
        for day in 0..self.problem.num_days {
            self.instance.graph.add_edge(
                self.instance.depot_target[day],
                self.instance.depot_source[day + 1],
                mcf::FlowValues::new_unconstrained(self.problem.depot.daily_cost),
            );

            for customer in 0..self.problem.num_customers {
                self.instance.graph.add_edge(
                    self.instance.customers[customer][day],
                    self.instance.customers[customer][day + 1],
                    mcf::FlowValues::new_unconstrained(self.problem.customers[customer].daily_cost),
                );
            }
        }
    }

    fn add_super_sink(&mut self) {
        let day = self.problem.num_days;
        let node = self.instance.graph.add_node(-self.flow_imbalance);

        self.instance.graph.add_edge(
            self.instance.depot_source[day],
            node,
            mcf::FlowValues::new_unconstrained(0.0),
        );

        for customer in 0..self.problem.num_customers {
            self.instance.graph.add_edge(
                self.instance.customers[customer][day],
                node,
                mcf::FlowValues::new_unconstrained(0.0),
            );
        }
    }
}

struct BranchAndBound {
    //instance: AuxiliaryMCFInstance,
}

impl BranchAndBound {
    fn solve(problem: Problem) {
        let instance = AuxiliaryMCFInstance::new(problem);

        let result = mcf::run(&instance.graph);

        match result {
            Ok(solution) => {
                println!("MCF solution is: {}", solution.cost);
            }
            Err(mcf::Error::Infeasible) => {
                println!("Instance is infeasible!");
            }
            Err(error) => {
                panic!(error);
            }
        }

        // TODO
    }
}

pub fn solve(problem: Problem) {
    BranchAndBound::solve(problem)
}
