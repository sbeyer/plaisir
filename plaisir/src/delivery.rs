use crate::problem::*;
use std::cmp::Ordering;

const PRINT_VARIABLE_VALUES: bool = false;
const MIP_EPSILON: f64 = 1e-7;

#[derive(Debug, Eq)]
pub struct Delivery {
    pub quantity: usize,
    pub customer: usize,
}

impl PartialEq for Delivery {
    fn eq(&self, other: &Self) -> bool {
        self.quantity == other.quantity
    }
}

impl PartialOrd for Delivery {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Delivery {
    fn cmp(&self, other: &Self) -> Ordering {
        self.quantity.cmp(&other.quantity)
    }
}

/// Vector indexed by day, vehicle and customer containing the number of delivered items
#[derive(Debug)]
pub struct Deliveries(Vec<Vec<Vec<usize>>>);

impl Deliveries {
    pub fn new(problem: &Problem) -> Self {
        Deliveries(vec![
            vec![
                vec![0; problem.num_customers];
                problem.num_vehicles
            ];
            problem.num_days
        ])
    }

    pub fn set(&mut self, t: usize, v: usize, i: usize, quantity: usize) {
        self.0[t][v][i - 1] = quantity;
    }

    pub fn get(&self, t: usize, v: usize, i: usize) -> usize {
        self.0[t][v][i - 1]
    }

    pub fn change_vehicle(&mut self, t: usize, from_v: usize, to_v: usize, i: usize) {
        let quantity = self.get(t, from_v, i);
        self.set(t, from_v, i, 0);
        self.set(t, to_v, i, quantity);
    }

    pub fn get_all_delivered_customers(&self, t: usize, v: usize) -> Vec<usize> {
        self.0[t][v]
            .iter()
            .enumerate()
            .filter(|(_, quantity)| *quantity > &0)
            .map(|(i, _)| i + 1)
            .collect()
    }
}

pub struct Solver<'a> {
    problem: &'a Problem,
    model: grb::Model,
    vars: Vec<Vec<Vec<grb::Var>>>,
}

impl<'a> Solver<'a> {
    pub fn new(env: &grb::Env, problem: &'a Problem) -> grb::Result<Self> {
        let mut model = grb::Model::with_env("deliveries", env)?;

        model.set_objective(0, grb::ModelSense::Minimize)?;

        let inventory_vars = problem
            .all_days()
            .map(|t| {
                problem
                    .all_sites()
                    .map(|i| {
                        let site = problem.site(i);
                        let name = format!("i_{t}_{i}");
                        let coeff = site.cost();
                        let bounds = site.level_bounds();
                        model
                            .add_var(
                                &name,
                                grb::VarType::Continuous,
                                coeff,
                                bounds.0,
                                bounds.1 + site.level_change(),
                                std::iter::empty(),
                            )
                            .unwrap()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let vars = problem
            .all_days()
            .map(|t| {
                problem
                    .all_vehicles()
                    .map(|v| {
                        problem
                            .all_customers()
                            .map(|i| {
                                let name = format!("d_{t}_{v}_{i}");
                                let coeff = 0.0;
                                let bounds = (0.0, 0.0);
                                model
                                    .add_var(
                                        &name,
                                        grb::VarType::Continuous,
                                        coeff,
                                        bounds.0,
                                        bounds.1,
                                        std::iter::empty(),
                                    )
                                    .unwrap()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        let mut solver = Solver {
            problem,
            model,
            vars,
        };
        // Use the original constraints from Solver here, too:

        // ensure vehicle capacity is not exceeded
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                let mut lhs = grb::expr::LinExpr::new();
                for i in problem.all_customers() {
                    lhs.add_term(1.0, solver.var(t, v, i));
                }

                solver.model.add_constr(
                    &format!("VC_{t}_{v}"),
                    grb::c!(lhs <= problem.capacity as f64),
                )?;
            }
        }

        // inventory flow for depot
        for t in problem.all_days() {
            let mut lhs = grb::expr::LinExpr::new();
            for v in problem.all_vehicles() {
                for j in problem.all_customers() {
                    lhs.add_term(-1.0, solver.var(t, v, j));
                }
            }
            lhs.add_term(-1.0, inventory_vars[t][0]); // outgoing inventory
            let depot = problem.site(0);
            let mut value = -depot.level_change();

            if t == 0 {
                value -= depot.level_start();
            } else {
                lhs.add_term(1.0, inventory_vars[t - 1][0]); // incoming inventory
            }

            solver
                .model
                .add_constr(&format!("Ifd_{t}"), grb::c!(lhs == value))?;
        }

        // inventory flow for customers
        for t in problem.all_days() {
            for i in problem.all_customers() {
                let mut lhs = grb::expr::LinExpr::new();
                for v in problem.all_vehicles() {
                    lhs.add_term(1.0, solver.var(t, v, i)); // delivered
                }
                lhs.add_term(-1.0, inventory_vars[t][i]); // outgoing inventory
                let customer = problem.site(i);
                let mut value = -customer.level_change();

                if t == 0 {
                    value -= customer.level_start();
                } else {
                    lhs.add_term(1.0, inventory_vars[t - 1][i]); // incoming inventory
                }

                solver
                    .model
                    .add_constr(&format!("Ifc_{t}_{i}"), grb::c!(lhs == value))?;
            }
        }

        Ok(solver)
    }

    pub fn var(&self, t: usize, v: usize, i: usize) -> grb::Var {
        self.vars[t][v][i - 1]
    }

    /// Sets whether the delivery at day `t` with vehicle `v` to customer `i` is active or not
    pub fn set_status(&mut self, t: usize, v: usize, i: usize, active: bool) -> grb::Result<()> {
        self.model.set_obj_attr(
            grb::attr::UB,
            &self.var(t, v, i),
            if active {
                self.problem.capacity as f64
            } else {
                0.0
            },
        )
    }

    /// Sets whether all deliveries are active or not using a function
    pub fn set_all_statuses<F>(&mut self, is_active: F) -> grb::Result<()>
    where
        F: Fn(usize, usize, usize) -> bool,
    {
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                for i in self.problem.all_customers() {
                    self.set_status(t, v, i, is_active(t, v, i))?;
                }
            }
        }
        Ok(())
    }

    /// Solve Minimum-Cost Flow LP to improve deliveries
    pub fn solve(&mut self) -> grb::Result<Option<Deliveries>> {
        let mut deliveries = Deliveries::new(self.problem);

        self.model.optimize()?;

        let status = self.model.status()?;
        if status == grb::Status::Optimal {
            if PRINT_VARIABLE_VALUES {
                eprintln!("#### MIP solution status: {status:?}");

                let objective = self.model.get_attr(grb::attr::ObjVal)?;
                eprintln!("#### MIP solution value: {objective}");

                for delivery_vars_per_day in self.vars.iter() {
                    for delivery_vars_per_customer in delivery_vars_per_day.iter() {
                        for var in delivery_vars_per_customer.iter() {
                            let name = self.model.get_obj_attr(grb::attr::VarName, var)?;
                            let value = self.model.get_obj_attr(grb::attr::X, var)?;
                            if value > MIP_EPSILON {
                                eprintln!("####   - {name}: {value}");
                            }
                        }
                    }
                }
            }

            if status == grb::Status::Optimal {
                for t in self.problem.all_days() {
                    for v in self.problem.all_vehicles() {
                        for i in self.problem.all_customers() {
                            let var = self.var(t, v, i);
                            let value = self.model.get_obj_attr(grb::attr::X, &var)?;

                            deliveries.set(t, v, i, value.round() as usize);
                        }
                    }
                }
            }

            Ok(Some(deliveries))
        } else {
            Ok(None)
        }
    }
}
