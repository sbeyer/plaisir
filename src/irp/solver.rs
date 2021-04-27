use super::*;
use grb::prelude as gurobi;

struct Variables<'a> {
    problem: &'a Problem,
    route: Vec<Vec<Vec<gurobi::Var>>>,
    carry: Vec<Vec<Vec<gurobi::Var>>>,
    deliver: Vec<Vec<gurobi::Var>>,
    inventory: Vec<Vec<gurobi::Var>>,
}

impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut gurobi::Model) -> Self {
        let mut vars = Variables {
            problem: &problem,
            route: Vec::with_capacity(problem.num_days),
            carry: Vec::with_capacity(problem.num_days),
            deliver: Vec::with_capacity(problem.num_days),
            inventory: Vec::with_capacity(problem.num_days),
        };

        // route variables
        for t in 0..problem.num_days {
            vars.route
                .push(Vec::with_capacity(problem.num_customers + 1));
            for i in 0..=problem.num_customers {
                vars.route[t].push(Vec::with_capacity(problem.num_customers));
                for j in 0..=problem.num_customers {
                    if i != j {
                        let name = String::from(format!("r_{}_{}_{}", t, i, j));
                        let coeff = problem.distance(i, j).into();
                        let bounds = (0.0, 1.0);
                        //let var = grb::add_binvar!(lp, name: &name, obj: coeff).unwrap();
                        let var = lp
                            .add_var(
                                &name,
                                gurobi::VarType::Binary,
                                coeff,
                                bounds.0,
                                bounds.1,
                                std::iter::empty(),
                            )
                            .unwrap();
                        vars.route[t][i].push(var);
                    }
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
                        let name = String::from(format!("c_{}_{}_{}", t, i, j));
                        let coeff = 0.0;
                        let bounds = (0.0, problem.capacity.into());
                        let var = lp
                            .add_var(
                                &name,
                                gurobi::VarType::Continuous,
                                coeff,
                                bounds.0,
                                bounds.1,
                                std::iter::empty(),
                            )
                            .unwrap();
                        vars.carry[t][i].push(var);
                    }
                }
            }
        }

        // inventory variables
        for t in 0..problem.num_days {
            vars.inventory
                .push(Vec::with_capacity(problem.num_customers + 1));
            for i in 0..=problem.num_customers {
                let name = String::from(format!("i_{}_{}", t, i));
                let coeff = problem.daily_cost(i);
                let bounds = problem.level_bounds(i);
                let var = lp
                    .add_var(
                        &name,
                        gurobi::VarType::Continuous,
                        coeff,
                        bounds.0,
                        bounds.1,
                        std::iter::empty(),
                    )
                    .unwrap();
                vars.inventory[t].push(var);
            }
        }

        // deliver variables
        for t in 0..problem.num_days {
            vars.deliver.push(Vec::with_capacity(problem.num_customers));
            for i in 1..=problem.num_customers {
                let name = String::from(format!("d_{}_{}", t, i));
                let coeff = 0.0;
                let bounds = (0.0, problem.capacity.into());
                let var = lp
                    .add_var(
                        &name,
                        gurobi::VarType::Continuous,
                        coeff,
                        bounds.0,
                        bounds.1,
                        std::iter::empty(),
                    )
                    .unwrap();
                vars.deliver[t].push(var);
            }
        }

        vars
    }

    fn route(&self, t: usize, i: usize, j: usize) -> gurobi::Var {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);
        debug_assert!(i != j);

        self.route[t][i][if j > i { j - 1 } else { j }]
    }

    fn carry(&self, t: usize, i: usize, j: usize) -> gurobi::Var {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);
        debug_assert!(i != j);

        self.carry[t][i][if j > i { j - 1 } else { j }]
    }

    fn inventory(&self, t: usize, i: usize) -> gurobi::Var {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);

        self.inventory[t][i]
    }

    fn deliver(&self, t: usize, i: usize) -> gurobi::Var {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i >= 1);
        debug_assert!(i <= self.problem.num_customers);

        self.deliver[t][i - 1]
    }
}

struct Delivery {
    quantity: usize,
    customer: usize,
}

struct Solution {
    routes: Vec<Vec<Vec<Delivery>>>,
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (t, routes) in self.routes.iter().enumerate() {
            writeln!(f, "Day {}", t + 1)?;
            for (route_idx, route) in routes.iter().enumerate() {
                write!(f, "Route {}: ", route_idx + 1)?;
                for route_stop in route.iter() {
                    write!(f, "{} ", route_stop.customer)?;
                    if route_stop.quantity != 0 {
                        write!(f, "( {} ) ", route_stop.quantity)?;
                    }
                    write!(f, "- ")?;
                }
                writeln!(f, "{}", route[0].customer)?;
            }
        }
        Ok(())
    }
}

struct Solver<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    lp: gurobi::Model,
}

impl<'a> Solver<'a> {
    fn solve(problem: &'a Problem) -> grb::Result<Self> {
        let mut lp = gurobi::Model::new("irp")?;
        lp.set_objective(0, gurobi::ModelSense::Minimize)?;
        let vars = Variables::new(&problem, &mut lp);

        // at most m vehicles for the routing
        for t in 0..problem.num_days {
            let mut lhs = grb::expr::LinExpr::new();
            for j in 1..=problem.num_customers {
                lhs.add_term(1.0, vars.route(t, 0, j));
            }

            lp.add_constr("Rmv", grb::c!(lhs <= problem.num_vehicles))?;
        }

        // route flow node-disjointness (at most one visit)
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                for j in 0..=problem.num_customers {
                    if i != j {
                        lhs.add_term(1.0, vars.route(t, j, i));
                    }
                }

                lp.add_constr("Rnd", grb::c!(lhs <= 1.0))?;
            }
        }

        // route flow conservation
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                for j in 0..=problem.num_customers {
                    if i != j {
                        lhs.add_term(1.0, vars.route(t, j, i));
                        lhs.add_term(-1.0, vars.route(t, i, j));
                    }
                }

                lp.add_constr("Rfc", grb::c!(lhs == 0.0))?;
            }
        }

        // glue: disable carry flow if we have no route flow
        for t in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                for j in 0..=problem.num_customers {
                    if j != i {
                        let mut lhs = grb::expr::LinExpr::new();
                        lhs.add_term(problem.capacity as f64, vars.route(t, i, j));
                        lhs.add_term(-1.0, vars.carry(t, i, j));

                        lp.add_constr("Gcr", grb::c!(lhs >= 0.0))?;
                    }
                }
            }
        }

        // don't carry too much
        let max_amount = problem.capacity as f64 * problem.num_vehicles as f64;
        for t in 0..problem.num_days {
            let mut lhs = grb::expr::LinExpr::new();
            for j in 1..=problem.num_customers {
                lhs.add_term(1.0, vars.carry(t, 0, j));
            }
            lp.add_constr("Clim", grb::c!(lhs <= max_amount))?;
        }

        // carry and deliver flow
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                for j in 0..=problem.num_customers {
                    if j != i {
                        lhs.add_term(1.0, vars.carry(t, j, i)); // incoming carry
                    }
                }
                for j in 1..=problem.num_customers {
                    if j != i {
                        lhs.add_term(-1.0, vars.carry(t, i, j)); // outgoing carry
                    }
                }
                lhs.add_term(-1.0, vars.deliver(t, i)); // deliver to customer

                lp.add_constr("CDf", grb::c!(lhs == 0.0))?;
            }
        }

        // inventory flow for depot
        for t in 0..problem.num_days {
            let mut lhs = grb::expr::LinExpr::new();
            for j in 1..=problem.num_customers {
                lhs.add_term(-1.0, vars.carry(t, 0, j)); // outgoing carry
            }
            lhs.add_term(-1.0, vars.inventory(t, 0)); // outgoing inventory
            let mut value = -problem.daily_level_change(0);

            if t == 0 {
                value -= problem.start_level(0);
            } else {
                lhs.add_term(1.0, vars.inventory(t - 1, 0)); // incoming inventory
            }

            lp.add_constr("Ifd", grb::c!(lhs == value))?;
        }

        // inventory flow for customers
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                lhs.add_term(1.0, vars.deliver(t, i)); // delivered
                lhs.add_term(-1.0, vars.inventory(t, i)); // outgoing inventory
                let mut value = -problem.daily_level_change(i);

                if t == 0 {
                    value -= problem.start_level(i);
                } else {
                    lhs.add_term(1.0, vars.inventory(t - 1, i)); // incoming inventory
                }

                lp.add_constr("Ifc", grb::c!(lhs == value))?;
            }
        }

        lp.optimize()?;

        Ok(Self {
            problem: &problem,
            vars,
            lp,
        })
    }

    fn print_raw_variable(&self, var: &gurobi::Var) -> grb::Result<()> {
        let name = self.lp.get_obj_attr(grb::attr::VarName, &var)?;
        let value = self.lp.get_obj_attr(grb::attr::X, &var)?;
        if value > 0. {
            println!("  - {}: {}", name, value);
        }
        Ok(())
    }

    fn print_raw_solution(&self) -> grb::Result<()> {
        let status = self.lp.status()?;
        println!("MIP solution status: {:?}", status);

        let objective = self.lp.get_attr(gurobi::attr::ObjVal)?;
        println!("MIP solution value: {}", objective);

        // print raw solution:
        for t in 0..self.problem.num_days {
            for i in 0..=self.problem.num_customers {
                for j in 0..=self.problem.num_customers {
                    if i != j {
                        self.print_raw_variable(&self.vars.route(t, i, j))?;
                        self.print_raw_variable(&self.vars.carry(t, i, j))?;
                    }
                }
            }

            for i in 0..=self.problem.num_customers {
                self.print_raw_variable(&self.vars.inventory(t, i))?;
            }

            for i in 1..=self.problem.num_customers {
                self.print_raw_variable(&self.vars.deliver(t, i))?;
            }
        }

        Ok(())
    }

    fn get_delivery_amount(&self, t: usize, source: usize, target: usize) -> grb::Result<usize> {
        let var_route = self.vars.route(t, source, target);
        let delivered = self.lp.get_obj_attr(grb::attr::X, &var_route)?;
        if delivered > 0. {
            let var_deliver = self.vars.deliver(t, target);
            let quantity = self.lp.get_obj_attr(grb::attr::X, &var_deliver)?;
            Ok(quantity.round() as usize)
        } else {
            Ok(0)
        }
    }

    fn get_solution(&self) -> grb::Result<Solution> {
        let mut sol = Solution {
            routes: Vec::with_capacity(self.problem.num_days),
        };
        for t in 0..self.problem.num_days {
            sol.routes.push(Vec::new());
            for _ in 0..self.problem.num_vehicles {
                sol.routes[t].push(vec![Delivery {
                    quantity: 0,
                    customer: 0,
                }]);
            }

            let mut visited = vec![false; self.problem.num_customers + 1];
            visited[0] = true;
            for route in 0..self.problem.num_vehicles {
                let mut found_something = false;
                loop {
                    let i = sol.routes[t][route].last().unwrap().customer;
                    let mut found = false;
                    for j in 1..=self.problem.num_customers {
                        if i != j && !visited[j] {
                            let delivered = self.get_delivery_amount(t, i, j)?;
                            if delivered > 0 {
                                visited[j] = true;
                                sol.routes[t][route].push(Delivery {
                                    quantity: delivered,
                                    customer: j.into(),
                                });
                                found = true;
                                found_something = true;
                                break;
                            }
                        }
                    }

                    if !found {
                        break;
                    }
                }

                if !found_something {
                    break;
                }
            }
        }

        Ok(sol)
    }
}

pub fn solve(problem: Problem) {
    let solver = Solver::solve(&problem).unwrap();
    solver.print_raw_solution().unwrap();
    let solution = solver.get_solution().unwrap();
    println!("{}", solution);
}
