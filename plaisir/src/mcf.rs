use crate::problem::*;

pub struct McfSubproblem<'a> {
    problem: &'a Problem,
    pub model: grb::Model,
    pub delivery_vars: Vec<Vec<Vec<grb::Var>>>,
}

impl<'a> McfSubproblem<'a> {
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
                        let name = format!("i_{}_{}", t, i);
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

        let delivery_vars = problem
            .all_days()
            .map(|t| {
                problem
                    .all_vehicles()
                    .map(|v| {
                        problem
                            .all_customers()
                            .map(|i| {
                                let name = format!("d_{}_{}_{}", t, v, i);
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

        let mut mcf = McfSubproblem {
            problem,
            model,
            delivery_vars,
        };
        // Use the original constraints from Solver here, too:

        // ensure vehicle capacity is not exceeded
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                let mut lhs = grb::expr::LinExpr::new();
                for i in problem.all_customers() {
                    lhs.add_term(1.0, mcf.delivery_var(t, v, i));
                }

                mcf.model.add_constr(
                    &format!("VC_{}_{}", t, v),
                    grb::c!(lhs <= problem.capacity as f64),
                )?;
            }
        }

        // inventory flow for depot
        for t in problem.all_days() {
            let mut lhs = grb::expr::LinExpr::new();
            for v in problem.all_vehicles() {
                for j in problem.all_customers() {
                    lhs.add_term(-1.0, mcf.delivery_var(t, v, j));
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

            mcf.model
                .add_constr(&format!("Ifd_{}", t), grb::c!(lhs == value))?;
        }

        // inventory flow for customers
        for t in problem.all_days() {
            for i in problem.all_customers() {
                let mut lhs = grb::expr::LinExpr::new();
                for v in problem.all_vehicles() {
                    lhs.add_term(1.0, mcf.delivery_var(t, v, i)); // delivered
                }
                lhs.add_term(-1.0, inventory_vars[t][i]); // outgoing inventory
                let customer = problem.site(i);
                let mut value = -customer.level_change();

                if t == 0 {
                    value -= customer.level_start();
                } else {
                    lhs.add_term(1.0, inventory_vars[t - 1][i]); // incoming inventory
                }

                mcf.model
                    .add_constr(&format!("Ifc_{}_{}", t, i), grb::c!(lhs == value))?;
            }
        }

        Ok(mcf)
    }

    pub fn delivery_var(&self, t: usize, v: usize, i: usize) -> grb::Var {
        self.delivery_vars[t][v][i - 1]
    }

    /// Sets whether the delivery at day `t` with vehicle `v` to customer `i` is active or not
    pub fn set_delivery_status(
        &mut self,
        t: usize,
        v: usize,
        i: usize,
        active: bool,
    ) -> grb::Result<()> {
        self.model.set_obj_attr(
            grb::attr::UB,
            &self.delivery_var(t, v, i),
            if active {
                self.problem.capacity as f64
            } else {
                0.0
            },
        )
    }
}
