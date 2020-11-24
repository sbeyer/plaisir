use super::*;

struct Variables<'a> {
    problem: &'a Problem,
    route: Vec<minilp::Variable>,
    carry: Vec<minilp::Variable>,
    deliver: Vec<minilp::Variable>,
    inventory: Vec<minilp::Variable>,
}

impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut minilp::Problem) -> Self {
        let mut vars = Variables {
            problem: &problem,
            route: Vec::new(),
            carry: Vec::new(),
            deliver: Vec::new(),
            inventory: Vec::new(),
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

        // TODO: more variables, more constraints

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
