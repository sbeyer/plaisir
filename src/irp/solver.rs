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

struct BranchAndBound {}

impl BranchAndBound {
    fn solve(problem: Problem) -> grb::Result<()> {
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

        let status = lp.status()?;
        println!("MIP solution status: {:?}", status);

        let objective = lp.get_attr(gurobi::attr::ObjVal)?;
        println!("MIP solution value: {}", objective);

        Ok(())
    }
}

pub fn solve(problem: Problem) {
    BranchAndBound::solve(problem).unwrap();
}
