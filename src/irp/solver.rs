use super::*;
use grb::prelude as gurobi;
use std::time;

struct Variables<'a> {
    problem: &'a Problem,
    variables: Vec<gurobi::Var>,
    route_range: (usize, usize),
    deliver_range: (usize, usize),
    visit_range: (usize, usize),
    inventory_range: (usize, usize),
}

impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut gurobi::Model) -> Self {
        let num_variables_route = problem.num_days
            * problem.num_vehicles
            * ((problem.num_customers + 1) * problem.num_customers / 2);
        let num_variables_visit =
            problem.num_days * problem.num_vehicles * (problem.num_customers + 1);
        let num_variables_inventory = problem.num_days * (problem.num_customers + 1);
        let num_variables_deliver = problem.num_days * problem.num_vehicles * problem.num_customers;
        let num_variables = num_variables_route
            + num_variables_visit
            + num_variables_inventory
            + num_variables_deliver;
        let mut vars = Variables {
            problem: &problem,
            variables: Vec::with_capacity(num_variables),
            route_range: (0, num_variables),
            deliver_range: (0, num_variables),
            visit_range: (0, num_variables),
            inventory_range: (0, num_variables),
        };

        // route variables
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    for j in i + 1..=problem.num_customers {
                        let name = format!("r_{}_{}_{}_{}", t, v, i, j);
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
                        debug_assert_eq!(vars.variables.len(), vars.route_index(t, v, i, j));
                        vars.variables.push(var);
                    }
                }
            }
        }
        vars.route_range.1 = vars.variables.len();
        debug_assert_eq!(num_variables_route, vars.route_range.1 - vars.route_range.0);

        // visit variables
        vars.visit_range.0 = vars.route_range.1;
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    let name = format!("v_{}_{}_{}", t, v, i);
                    let coeff = 0.0;
                    let bounds = (0.0, 1.0);
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
                    debug_assert_eq!(vars.variables.len(), vars.visit_index(t, v, i));
                    vars.variables.push(var);
                }
            }
        }
        vars.visit_range.1 = vars.variables.len();
        debug_assert_eq!(num_variables_visit, vars.visit_range.1 - vars.visit_range.0);

        // inventory variables
        vars.inventory_range.0 = vars.visit_range.1;
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
                debug_assert_eq!(vars.variables.len(), vars.inventory_index(t, i));
                vars.variables.push(var);
            }
        }
        vars.inventory_range.1 = vars.variables.len();
        debug_assert_eq!(
            num_variables_inventory,
            vars.inventory_range.1 - vars.inventory_range.0
        );

        // deliver variables
        vars.deliver_range.0 = vars.inventory_range.1;
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 1..=problem.num_customers {
                    let name = format!("d_{}_{}_{}", t, v, i);
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
                    debug_assert_eq!(vars.variables.len(), vars.deliver_index(t, v, i));
                    vars.variables.push(var);
                }
            }
        }
        vars.deliver_range.1 = vars.variables.len();
        debug_assert_eq!(
            num_variables_deliver,
            vars.deliver_range.1 - vars.deliver_range.0
        );

        vars
    }

    fn route_index(&self, t: usize, v: usize, i: usize, j: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(v < self.problem.num_vehicles);
        debug_assert!(i < j);
        let n = self.problem.num_customers;
        debug_assert!(j <= n);

        let i_j_block_size = n * (n + 1) / 2;
        let v_i_j_block_size = self.problem.num_vehicles * i_j_block_size;
        let offset = self.route_range.0 + t * v_i_j_block_size + v * i_j_block_size;
        let i_offset = i * (2 * n - i + 1) / 2;
        let i_j_offset = i_offset + j - i - 1;

        let result = offset + i_j_offset;
        debug_assert!(result < self.route_range.1);

        result
    }

    fn route(&self, t: usize, v: usize, i: usize, j: usize) -> gurobi::Var {
        let index = if i < j {
            self.route_index(t, v, i, j)
        } else {
            self.route_index(t, v, j, i)
        };
        self.variables[index]
    }

    fn visit_index(&self, t: usize, v: usize, i: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(v < self.problem.num_vehicles);
        debug_assert!(i <= self.problem.num_customers);

        let i_block_size = self.problem.num_customers + 1;
        let v_i_block_size = self.problem.num_vehicles * i_block_size;
        let offset = self.visit_range.0 + t * v_i_block_size + v * i_block_size;
        let result = offset + i;
        debug_assert!(result < self.visit_range.1);

        result
    }

    fn visit(&self, t: usize, v: usize, i: usize) -> gurobi::Var {
        self.variables[self.visit_index(t, v, i)]
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

    fn deliver_index(&self, t: usize, v: usize, i: usize) -> usize {
        debug_assert!(t < self.problem.num_days);
        debug_assert!(v < self.problem.num_vehicles);
        debug_assert!(i >= 1);
        debug_assert!(i <= self.problem.num_customers);

        let i_block_size = self.problem.num_customers;
        let v_i_block_size = self.problem.num_vehicles * i_block_size;
        let offset = self.deliver_range.0 + t * v_i_block_size + v * i_block_size;
        let result = offset + i - 1;
        debug_assert!(result < self.deliver_range.1);

        result
    }

    fn deliver(&self, t: usize, v: usize, i: usize) -> gurobi::Var {
        self.variables[self.deliver_index(t, v, i)]
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
            #[allow(clippy::needless_range_loop)]
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
        writeln!(f, "{:.2}", self.cost_inventory_customers)?;
        writeln!(f, "{:.2}", self.cost_inventory_depot)?;
        writeln!(f, "{:.2}", self.cost_total)?;

        // Meta
        writeln!(f, "{}", self.processor)?;
        writeln!(f, "{}", self.time.as_millis() as f64 * 1e-3)?;

        Ok(())
    }
}

struct SolverData<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    varnames: Vec<String>,
    start_time: time::Instant,
    cpu: String,
    ncalls: usize,
}

impl<'a> SolverData<'a> {
    const EPSILON: f64 = 1e-7;

    fn new(problem: &'a Problem, lp: &mut gurobi::Model, cpu: String) -> Self {
        let start_time = time::Instant::now();
        let vars = Variables::new(&problem, lp);
        lp.update().unwrap(); // update to access variable names
        let varnames = vars
            .variables
            .iter()
            .map(|var| lp.get_obj_attr(grb::attr::VarName, &var).unwrap())
            .collect();
        SolverData {
            problem,
            vars,
            varnames,
            start_time,
            cpu,
            ncalls: 0,
        }
    }

    fn subtour_elimination<F>(&self, assignment: Vec<f64>, add: F) -> grb::Result<()>
    where
        F: Fn(grb::constr::IneqExpr) -> grb::Result<()>,
    {
        self.varnames
            .iter()
            .zip(assignment.iter())
            .filter(|(_, &value)| value > Self::EPSILON)
            .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));

        // build auxiliary LPs based on model
        for t in 0..self.problem.num_days {
            for v in 0..self.problem.num_vehicles {
                for k in 1..=self.problem.num_customers {
                    let mut lp = gurobi::Model::new(&format!("sep_{}_{}_{}", t, v, k))?;
                    lp.set_param(grb::param::Threads, 1)?;
                    lp.set_objective(0, gurobi::ModelSense::Maximize)?;

                    // generate variables
                    let mut z_vars =
                        Vec::<gurobi::Var>::with_capacity(self.problem.num_customers + 1);
                    let mut y_vars =
                        Vec::<Vec<gurobi::Var>>::with_capacity(self.problem.num_customers);
                    for i in 0..=self.problem.num_customers {
                        let name = format!("z_{}", i);
                        let coeff = if i != k {
                            -assignment[self.vars.visit_index(t, v, i)]
                        } else {
                            0.0
                        };
                        let var = lp.add_var(
                            &name,
                            gurobi::VarType::Continuous,
                            coeff,
                            0.0,
                            1.0,
                            std::iter::empty(),
                        )?;
                        z_vars.push(var);

                        let mut y_i_vars =
                            Vec::<gurobi::Var>::with_capacity(self.problem.num_customers - i);

                        for j in i + 1..=self.problem.num_customers {
                            let name = format!("y_{}_{}", i, j);
                            let coeff = assignment[self.vars.route_index(t, v, i, j)];
                            let var = lp.add_var(
                                &name,
                                gurobi::VarType::Continuous,
                                coeff,
                                0.0,
                                1.0,
                                std::iter::empty(),
                            )?;
                            y_i_vars.push(var);
                        }

                        y_vars.push(y_i_vars);
                    }

                    let z_var = |i: usize| z_vars[i];
                    let y_var = |i: usize, j: usize| y_vars[i][j - i - 1];

                    // add constraint: k is contained
                    lp.add_constr("K", grb::c!(z_var(k) == 1.0))?;

                    // add optional constraint: multiplication bound
                    for i in 0..=self.problem.num_customers {
                        for j in i + 1..=self.problem.num_customers {
                            lp.add_constr(
                                &format!("M_{}_{}", i, j),
                                grb::c!(z_var(i) + z_var(j) - y_var(i, j) <= 1),
                            )?;
                        }
                    }

                    // add constraint: i and j bounds
                    for i in 0..=self.problem.num_customers {
                        for j in i + 1..=self.problem.num_customers {
                            lp.add_constr(
                                &format!("I_{}_{}", i, j),
                                grb::c!(z_var(i) - y_var(i, j) >= 0),
                            )?;
                            lp.add_constr(
                                &format!("J_{}_{}", i, j),
                                grb::c!(z_var(j) - y_var(i, j) >= 0),
                            )?;
                        }
                    }

                    // solve separation problem
                    lp.optimize()?;

                    let status = lp.status()?;
                    eprintln!("# Sep({}, {}, {}) solution status: {:?}", t, v, k, status);

                    let objective = lp.get_attr(gurobi::attr::ObjVal)?;
                    eprintln!("# Sep({}, {}, {}) solution value: {}", t, v, k, objective);

                    if objective > SolverData::EPSILON {
                        let mut lhs = grb::expr::LinExpr::new();

                        for i in 0..=self.problem.num_customers {
                            for j in i + 1..=self.problem.num_customers {
                                let var = y_var(i, j);
                                let var_name = lp.get_obj_attr(grb::attr::VarName, &var)?;
                                let var_value = lp.get_obj_attr(grb::attr::X, &var)?;
                                if var_value > SolverData::EPSILON {
                                    eprintln!("#   - {}: {}", var_name, var_value);
                                    debug_assert!(var_value > 1.0 - SolverData::EPSILON);
                                    lhs.add_term(1.0, self.vars.route(t, v, i, j));
                                }
                            }

                            if i != k {
                                let var = z_var(i);
                                let var_name = lp.get_obj_attr(grb::attr::VarName, &var)?;
                                let var_value = lp.get_obj_attr(grb::attr::X, &var)?;
                                if var_value > SolverData::EPSILON {
                                    eprintln!("#   - {}: {}", var_name, var_value);
                                    debug_assert!(var_value > 1.0 - SolverData::EPSILON);
                                    lhs.add_term(-1.0, self.vars.visit(t, v, i));
                                }
                            }
                        }

                        add(grb::c!(lhs <= 0))?;
                    }
                }
            }
        }

        Ok(())
    }

    fn elapsed_time(&self) -> time::Duration {
        self.start_time.elapsed()
    }

    fn is_delivered(
        &self,
        solution: &[f64],
        t: usize,
        v: usize,
        source: usize,
        target: usize,
    ) -> bool {
        let var_route = self.vars.route_index(t, v, source, target);
        solution[var_route] > Self::EPSILON
    }

    fn get_delivery_amount(&self, solution: &[f64], t: usize, target: usize) -> usize {
        let mut result = 0.;
        for v in 0..self.problem.num_vehicles {
            let var_deliver = self.vars.deliver_index(t, v, target);
            result += solution[var_deliver];
        }
        result.round() as usize
    }

    /*
    fn get_routes(&self, solution: &[f64]) -> Routes {
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
                let mut i = 0; // last visited location
                loop {
                    let mut found = false;

                    #[allow(clippy::needless_range_loop)]
                    for j in 1..=self.problem.num_customers {
                        if i != j && !visited[j] && self.is_delivered(&solution, t, i, j) {
                            let quantity = self.get_delivery_amount(&solution, t, j);
                            visited[j] = true;
                            i = j;
                            found = true;

                            if quantity > 0 {
                                routes[t][route].push(Delivery {
                                    quantity,
                                    customer: j,
                                });
                            }
                            break;
                        }
                    }

                    if !found {
                        break;
                    }
                }

                if i == 0 {
                    // nothing found, no need to check further routes
                    break;
                }
            }
        }

        routes
    }
    */
}

impl<'a> grb::callback::Callback for SolverData<'a> {
    fn callback(&mut self, w: gurobi::Where) -> grb::callback::CbResult {
        match w {
            gurobi::Where::MIPSol(ctx) => {
                self.ncalls += 1;
                let assignment = ctx.get_solution(&self.vars.variables)?;
                /*
                let routes = self.get_routes(&assignment);
                let solution =
                    Solution::new(&self.problem, routes, self.elapsed_time(), self.cpu.clone());
                */
                eprintln!("# Incumbent {}!", self.ncalls);
                eprintln!("#    current obj: {}", ctx.obj()?);
                eprintln!("#       best obj: {}", ctx.obj_best()?);

                self.varnames
                    .iter()
                    .zip(assignment.iter())
                    .filter(|(_, &value)| value > Self::EPSILON)
                    .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));

                //eprintln!("{}", solution);

                self.subtour_elimination(assignment, |constr| ctx.add_lazy(constr))?;
            }
            gurobi::Where::MIPNode(ctx) => {
                let status = ctx.status()?;
                eprintln!(
                    "# MIPNode {} {} {:?}",
                    ctx.sol_cnt()?,
                    ctx.node_cnt()?,
                    status
                );
                eprintln!("#       best objective: {}", ctx.obj_best()?);
                eprintln!("#       best obj bound: {}", ctx.obj_bnd()?);

                if status == grb::Status::Optimal {
                    let assignment = ctx.get_solution(&self.vars.variables)?;
                    self.subtour_elimination(assignment, |constr| ctx.add_lazy(constr))?;
                }
            }
            _ => (),
        }

        Ok(())
    }
}

struct Solver {}

impl Solver {
    fn solve(problem: &Problem, cpu: String) -> grb::Result<()> {
        let mut lp = gurobi::Model::new("irp")?;
        lp.set_param(grb::param::LazyConstraints, 1)?;
        lp.set_param(grb::param::Threads, 1)?;
        lp.set_objective(0, gurobi::ModelSense::Minimize)?;

        let mut data = SolverData::new(&problem, &mut lp, cpu);

        // route degree constraints
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-2.0, data.vars.visit(t, v, i));
                    for j in 0..=problem.num_customers {
                        if i != j {
                            lhs.add_term(1.0, data.vars.route(t, v, i, j));
                        }
                    }

                    lp.add_constr(&format!("Rd_{}_{}_{}", t, v, i), grb::c!(lhs == 0.0))?;
                }
            }
        }

        // glue: visit depot if customer is visited
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 1..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(1.0, data.vars.visit(t, v, 0));
                    lhs.add_term(-1.0, data.vars.visit(t, v, i));

                    lp.add_constr(&format!("Gdc_{}_{}_{}", t, v, i), grb::c!(lhs >= 0.0))?;
                }
            }
        }

        // glue: disable delivery if we do not visit
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 1..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(problem.capacity as f64, data.vars.visit(t, v, i));
                    lhs.add_term(-1.0, data.vars.deliver(t, v, i));

                    lp.add_constr(&format!("Gdv_{}_{}_{}", t, v, i), grb::c!(lhs >= 0.0))?;
                }
            }
        }

        // inventory flow for depot
        for t in 0..problem.num_days {
            let mut lhs = grb::expr::LinExpr::new();
            for v in 0..problem.num_vehicles {
                for j in 1..=problem.num_customers {
                    lhs.add_term(-1.0, data.vars.deliver(t, v, j));
                }
            }
            lhs.add_term(-1.0, data.vars.inventory(t, 0)); // outgoing inventory
            let mut value = -problem.daily_level_change(0);

            if t == 0 {
                value -= problem.start_level(0);
            } else {
                lhs.add_term(1.0, data.vars.inventory(t - 1, 0)); // incoming inventory
            }

            lp.add_constr(&format!("Ifd_{}", t), grb::c!(lhs == value))?;
        }

        // inventory flow for customers
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                for v in 0..problem.num_vehicles {
                    lhs.add_term(1.0, data.vars.deliver(t, v, i)); // delivered
                }
                lhs.add_term(-1.0, data.vars.inventory(t, i)); // outgoing inventory
                let mut value = -problem.daily_level_change(i);

                if t == 0 {
                    value -= problem.start_level(i);
                } else {
                    lhs.add_term(1.0, data.vars.inventory(t - 1, i)); // incoming inventory
                }

                lp.add_constr(&format!("Ifc_{}_{}", t, i), grb::c!(lhs == value))?;
            }
        }

        //lp.write("/tmp/foo.lp")?;
        lp.optimize_with_callback(&mut data)?;

        Self::print_raw_solution(&data, &lp)?;

        /*
        let assignment = lp.get_obj_attr_batch(grb::attr::X, data.vars.variables.clone())?;
        let routes = data.get_routes(&assignment);
        let solution = Solution::new(&problem, routes, data.elapsed_time(), data.cpu);
        eprintln!("# Final solution");
        eprintln!("{}", solution);
        */

        Ok(())
    }

    fn print_raw_solution(data: &SolverData, lp: &gurobi::Model) -> grb::Result<()> {
        let status = lp.status()?;
        eprintln!("# MIP solution status: {:?}", status);

        let objective = lp.get_attr(gurobi::attr::ObjVal)?;
        eprintln!("# MIP solution value: {}", objective);

        // print raw solution:
        for var in data.vars.variables.iter() {
            let name = lp.get_obj_attr(grb::attr::VarName, &var)?;
            let value = lp.get_obj_attr(grb::attr::X, &var)?;
            if value > SolverData::EPSILON {
                eprintln!("#   - {}: {}", name, value);
            }
        }

        Ok(())
    }
}

pub fn solve(problem: Problem, cpu: String) {
    Solver::solve(&problem, cpu).unwrap();
}
