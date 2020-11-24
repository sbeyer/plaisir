use super::*;

struct Variables<'a> {
    problem: &'a Problem,
    route: Vec<minilp::Variable>,
    carry: Vec<Vec<Vec<minilp::Variable>>>,
    deliver: Vec<Vec<minilp::Variable>>,
    inventory: Vec<Vec<minilp::Variable>>,
}

impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut minilp::Problem) -> Self {
        let mut vars = Variables {
            problem: &problem,
            route: Vec::new(),
            carry: Vec::with_capacity(problem.num_days),
            deliver: Vec::with_capacity(problem.num_days),
            inventory: Vec::with_capacity(problem.num_days),
        };

        // route variables
        for _ in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                for j in 0..=problem.num_customers {
                    vars.route
                        .push(lp.add_var(problem.distance(i, j).into(), (0.0, 1.0)));
                }
            }
        }

        // carry variables
        for t in 0..problem.num_days {
            vars.carry
                .push(Vec::with_capacity(problem.num_customers + 1));
            for i in 0..=problem.num_customers {
                vars.carry[t].push(Vec::with_capacity(problem.num_customers));
                for j in 0..=problem.num_customers {
                    if i != j {
                        vars.carry[t][i].push(lp.add_var(0.0, (0.0, problem.capacity.into())));
                    }
                }
            }
        }

        // inventory variables
        for t in 0..problem.num_days {
            vars.inventory
                .push(Vec::with_capacity(problem.num_customers + 1));
            for i in 0..=problem.num_customers {
                vars.inventory[t].push(lp.add_var(problem.daily_cost(i), problem.level_bounds(i)));
            }
        }

        // deliver variables
        for t in 0..problem.num_days {
            vars.deliver.push(Vec::with_capacity(problem.num_customers));
            for _ in 1..=problem.num_customers {
                vars.deliver[t].push(lp.add_var(0.0, (0.0, problem.capacity.into())));
            }
        }

        vars
    }

    fn route(&self, t: usize, i: usize, j: usize) -> minilp::Variable {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);

        let targetsize = self.problem.num_customers + 1;
        debug_assert!(j < targetsize);
        let sourcesize = targetsize * targetsize;
        debug_assert!(i * targetsize + j < sourcesize);

        self.route[t * sourcesize + i * targetsize + j]
    }

    fn carry(&self, t: usize, i: usize, j: usize) -> minilp::Variable {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);
        debug_assert!(i != j);

        self.carry[t][i][if j > i { j - 1 } else { j }]
    }

    fn inventory(&self, t: usize, i: usize) -> minilp::Variable {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);

        self.inventory[t][i]
    }

    fn deliver(&self, t: usize, i: usize) -> minilp::Variable {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i >= 1);
        debug_assert!(i <= self.problem.num_customers);

        self.deliver[t][i - 1]
    }
}

struct BranchAndBound {}

impl BranchAndBound {
    fn solve(problem: Problem) {
        let mut lp = minilp::Problem::new(minilp::OptimizationDirection::Minimize);
        let vars = Variables::new(&problem, &mut lp);

        // at most m vehicles for the routing
        for t in 0..problem.num_days {
            let mut lhs = minilp::LinearExpr::empty();
            for j in 1..=problem.num_customers {
                lhs.add(vars.route(t, 0, j), 1.0);
            }

            lp.add_constraint(lhs, minilp::ComparisonOp::Le, problem.num_vehicles as f64);
            //lp.add_constraint(lhs, minilp::ComparisonOp::Eq, -instance[node]) // geq M
        }

        // route flow at arriving zone of each customer
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = minilp::LinearExpr::empty();
                for j in 0..=problem.num_customers {
                    let coeff = if i != j { 1.0 } else { -1.0 };
                    lhs.add(vars.route(t, j, i), coeff);
                }

                lp.add_constraint(lhs, minilp::ComparisonOp::Eq, 0.0);
            }
        }

        // route flow at leaving zone of each customer
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = minilp::LinearExpr::empty();
                for j in 0..=problem.num_customers {
                    let coeff = if i != j { -1.0 } else { 1.0 };
                    lhs.add(vars.route(t, i, j), coeff);
                }

                lp.add_constraint(lhs, minilp::ComparisonOp::Eq, 0.0);
            }
        }

        // glue: disable carry flow if we have no route flow
        for t in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                for j in 0..=problem.num_customers {
                    if j != i {
                        let mut lhs = minilp::LinearExpr::empty();
                        lhs.add(vars.route(t, i, j), problem.capacity as f64);
                        lhs.add(vars.carry(t, i, j), -1.0);

                        lp.add_constraint(lhs, minilp::ComparisonOp::Ge, 0.0);
                    }
                }
            }
        }

        // carry and deliver flow
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = minilp::LinearExpr::empty();
                for j in 0..=problem.num_customers {
                    if j != i {
                        lhs.add(vars.carry(t, j, i), 1.0); // incoming carry
                    }
                }
                for j in 1..=problem.num_customers {
                    if j != i {
                        lhs.add(vars.carry(t, i, j), -1.0); // outgoing carry
                    }
                }
                lhs.add(vars.deliver(t, i), -1.0); // deliver to customer

                lp.add_constraint(lhs, minilp::ComparisonOp::Eq, 0.0);
            }
        }

        // inventory flow for depot
        for t in 0..problem.num_days {
            let mut lhs = minilp::LinearExpr::empty();
            for j in 1..=problem.num_customers {
                lhs.add(vars.carry(t, 0, j), -1.0); // outgoing carry
            }
            lhs.add(vars.inventory(t, 0), -1.0); // outgoing inventory
            let mut value = -problem.daily_level_change(0);

            if t == 0 {
                value -= problem.start_level(0);
            } else {
                lhs.add(vars.inventory(t - 1, 0), 1.0); // incoming inventory
            }

            lp.add_constraint(lhs, minilp::ComparisonOp::Eq, value);
        }

        // inventory flow for customers
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = minilp::LinearExpr::empty();
                lhs.add(vars.deliver(t, i), 1.0); // delivered
                lhs.add(vars.inventory(t, i), -1.0); // outgoing inventory
                let mut value = -problem.daily_level_change(i);

                if t == 0 {
                    value -= problem.start_level(i);
                } else {
                    lhs.add(vars.inventory(t - 1, i), 1.0); // incoming inventory
                }

                lp.add_constraint(lhs, minilp::ComparisonOp::Eq, value);
            }
        }

        // finalize flow; not necessary
        /*{
            let mut lhs = minilp::LinearExpr::empty();
            let t = problem.num_days - 1;
            for i in 0..=problem.num_customers {
                lhs.add(vars.inventory(t, i), 1.0); // incoming inventory
            }
            let mut value = 0 as f64;
            for i in 0..=problem.num_customers {
                value += problem.daily_level_change(i);
            }
            value *= t as f64;
            for i in 0..=problem.num_customers {
                value += problem.start_level(i);
            }

            lp.add_constraint(lhs, minilp::ComparisonOp::Eq, value);
        }*/

        let result = lp.solve();

        match result {
            Ok(solution) => {
                println!("MIP solution is: {}", solution.objective());
            }
            Err(minilp::Error::Infeasible) => {
                println!("Instance is infeasible!");
            }
            Err(error) => {
                panic!(error);
            }
        }
    }
}

pub fn solve(problem: Problem) {
    BranchAndBound::solve(problem)
}
