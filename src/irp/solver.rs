use super::*;
use grb::prelude as gurobi;
use std::time;

struct Variables<'a> {
    problem: &'a Problem,
    variables: Vec<gurobi::Var>,
    route_range: (usize, usize),
    deliver_range: (usize, usize),
    carry_range: (usize, usize),
    inventory_range: (usize, usize),
}

impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut gurobi::Model) -> Self {
        // Number of variables necessary for...
        //  # route: problem.num_days * (problem.num_customers + 1) * problem.num_customers
        //  # carry: problem.num_days * (problem.num_customers + 1) * problem.num_customers
        //  # inventory: problem.num_days * (problem.num_customers + 1)
        //  # deliver: problem.num_days * problem.num_customers
        // This is in sum:
        let num_variables =
            ((2 * problem.num_customers + 4) * problem.num_customers + 1) * problem.num_days;
        let mut vars = Variables {
            problem: &problem,
            variables: Vec::with_capacity(num_variables),
            route_range: (0, num_variables),
            deliver_range: (0, num_variables),
            carry_range: (0, num_variables),
            inventory_range: (0, num_variables),
        };

        // route variables
        for t in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                for j in 0..=problem.num_customers {
                    if i != j {
                        let name = format!("r_{}_{}_{}", t, i, j);
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
                        debug_assert!(vars.variables.len() == vars.route_index(t, i, j));
                        vars.variables.push(var);
                    }
                }
            }
        }
        vars.route_range.1 = vars.variables.len();

        // carry variables
        vars.carry_range.0 = vars.route_range.1;
        for t in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                for j in 0..=problem.num_customers {
                    if i != j {
                        let name = format!("c_{}_{}_{}", t, i, j);
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
                        debug_assert!(vars.variables.len() == vars.carry_index(t, i, j));
                        vars.variables.push(var);
                    }
                }
            }
        }
        vars.carry_range.1 = vars.variables.len();

        // inventory variables
        vars.inventory_range.0 = vars.carry_range.1;
        for t in 0..problem.num_days {
            for i in 0..=problem.num_customers {
                let name = format!("i_{}_{}", t, i);
                let coeff = problem.daily_cost(i);
                let bounds = problem.level_bounds(i);
                let var = lp
                    .add_var(
                        &name,
                        gurobi::VarType::Continuous,
                        coeff,
                        bounds.0,
                        bounds.1 + problem.daily_level_change(i),
                        std::iter::empty(),
                    )
                    .unwrap();
                debug_assert!(vars.variables.len() == vars.inventory_index(t, i));
                vars.variables.push(var);
            }
        }
        vars.inventory_range.1 = vars.variables.len();

        // deliver variables
        vars.deliver_range.0 = vars.inventory_range.1;
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let name = format!("d_{}_{}", t, i);
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
                debug_assert!(vars.variables.len() == vars.deliver_index(t, i));
                vars.variables.push(var);
            }
        }
        vars.deliver_range.1 = vars.variables.len();

        vars
    }

    fn route_index(&self, t: usize, i: usize, j: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);
        debug_assert!(i != j);

        let j_block_size = self.problem.num_customers;
        let i_j_block_size = (self.problem.num_customers + 1) * j_block_size;
        let offset = self.route_range.0 + t * i_j_block_size + i * j_block_size;
        let result = offset + (if j > i { j - 1 } else { j });
        debug_assert!(result < self.route_range.1);

        result
    }

    fn route(&self, t: usize, i: usize, j: usize) -> gurobi::Var {
        self.variables[self.route_index(t, i, j)]
    }

    fn carry_index(&self, t: usize, i: usize, j: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);
        debug_assert!(j <= self.problem.num_customers);
        debug_assert!(i != j);

        let j_block_size = self.problem.num_customers;
        let i_j_block_size = (self.problem.num_customers + 1) * j_block_size;
        let offset = self.carry_range.0 + t * i_j_block_size + i * j_block_size;
        let result = offset + (if j > i { j - 1 } else { j });
        debug_assert!(result < self.carry_range.1);

        result
    }

    fn carry(&self, t: usize, i: usize, j: usize) -> gurobi::Var {
        self.variables[self.carry_index(t, i, j)]
    }

    fn inventory_index(&self, t: usize, i: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);

        let i_block_size = self.problem.num_customers + 1;
        let offset = self.inventory_range.0 + t * i_block_size;
        let result = offset + i;
        debug_assert!(result < self.inventory_range.1);

        result
    }

    fn inventory(&self, t: usize, i: usize) -> gurobi::Var {
        self.variables[self.inventory_index(t, i)]
    }

    fn deliver_index(&self, t: usize, i: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(i >= 1);
        debug_assert!(i <= self.problem.num_customers);

        let i_block_size = self.problem.num_customers;
        let offset = self.deliver_range.0 + t * i_block_size;
        let result = offset + i - 1;
        debug_assert!(result < self.deliver_range.1);

        result
    }

    fn deliver(&self, t: usize, i: usize) -> gurobi::Var {
        self.variables[self.deliver_index(t, i)]
    }
}

struct Delivery {
    quantity: usize,
    customer: usize,
}

type Routes = Vec<Vec<Vec<Delivery>>>;

struct Solution {
    routes: Routes,
    cost_transportation: f64,
    cost_inventory_depot: f64,
    cost_inventory_customers: f64,
    cost_total: f64,
    processor: String,
    time: time::Duration,
}

impl Solution {
    fn new(problem: &Problem, routes: Routes, time: time::Duration, cpu: String) -> Self {
        let mut sol = Solution {
            routes,
            cost_transportation: 0.,
            cost_inventory_depot: 0.,
            cost_inventory_customers: 0.,
            cost_total: 0.,
            processor: cpu,
            time,
        };

        // transportation cost
        for day_routes in sol.routes.iter() {
            for route in day_routes.iter() {
                let mut tour: Vec<usize> = route.iter().map(|x| x.customer).collect();
                if tour.len() > 1 {
                    tour.push(0);
                    for path in tour.windows(2) {
                        sol.cost_transportation += problem.distance(path[0], path[1]) as f64
                    }
                }
            }
        }

        // inventory cost
        let mut inventory: Vec<f64> = (0..=problem.num_customers)
            .map(|x| problem.start_level(x))
            .collect();
        let mut cost_inventory = vec![0.; problem.num_customers + 1];

        for day_routes in sol.routes.iter() {
            // step one: deliveries
            for route in day_routes.iter() {
                for Delivery { quantity, customer } in route.iter() {
                    inventory[*customer] += *quantity as f64;
                    inventory[0] -= *quantity as f64;
                }
            }

            // step two: daily change (production at depot, consumption at customers)
            for i in 0..=problem.num_customers {
                inventory[i] += problem.daily_level_change(i);
            }

            // update inventory costs
            for i in 0..=problem.num_customers {
                cost_inventory[i] += problem.daily_cost(i) * inventory[i]
            }
        }

        sol.cost_inventory_depot = cost_inventory[0];
        sol.cost_inventory_customers = cost_inventory[1..].iter().sum();

        // total cost
        sol.cost_total =
            sol.cost_transportation + sol.cost_inventory_depot + sol.cost_inventory_customers;

        sol
    }
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // Routes
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

        // Costs
        writeln!(f, "{}", self.cost_transportation)?;
        writeln!(f, "{}", self.cost_inventory_customers)?;
        writeln!(f, "{}", self.cost_inventory_depot)?;
        writeln!(f, "{}", self.cost_total)?;

        // Meta
        writeln!(f, "{}", self.processor)?;
        writeln!(f, "{}", self.time.as_secs())?;

        Ok(())
    }
}

struct Solver<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    lp: gurobi::Model,
    start_time: time::Instant,
}

impl<'a> Solver<'a> {
    fn solve(problem: &'a Problem) -> grb::Result<Self> {
        let start_time = time::Instant::now();
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
            start_time,
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

    fn get_solution_routes(&self) -> grb::Result<Routes> {
        let mut routes = Vec::with_capacity(self.problem.num_days);

        for t in 0..self.problem.num_days {
            routes.push(Vec::new());
            for _ in 0..self.problem.num_vehicles {
                routes[t].push(vec![Delivery {
                    quantity: 0,
                    customer: 0,
                }]);
            }

            let mut visited = vec![false; self.problem.num_customers + 1];
            visited[0] = true;
            for route in 0..self.problem.num_vehicles {
                let mut found_something = false;
                loop {
                    let i = routes[t][route].last().unwrap().customer;
                    let mut found = false;
                    for j in 1..=self.problem.num_customers {
                        if i != j && !visited[j] {
                            let delivered = self.get_delivery_amount(t, i, j)?;
                            if delivered > 0 {
                                visited[j] = true;
                                routes[t][route].push(Delivery {
                                    quantity: delivered,
                                    customer: j,
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

        Ok(routes)
    }
}

pub fn solve(problem: Problem, cpu: String) {
    let solver = Solver::solve(&problem).unwrap();
    solver.print_raw_solution().unwrap();
    let routes = solver.get_solution_routes().unwrap();
    let solution = Solution::new(&problem, routes, solver.start_time.elapsed(), cpu);
    println!("{}", solution);
}
