use super::*;
use grb::prelude as gurobi;
use std::time;

extern crate partitions;

struct Variables<'a> {
    problem: &'a Problem,
    variables: Vec<gurobi::Var>,
    route_range: (usize, usize),
    deliver_range: (usize, usize),
    visit_range: (usize, usize),
    inventory_range: (usize, usize),
}

#[allow(clippy::many_single_char_names)]
impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut gurobi::Model) -> Self {
        let num_variables_route = problem.num_days
            * problem.num_vehicles
            * (problem.num_customers + 1)
            * problem.num_customers;
        let num_variables_visit =
            problem.num_days * problem.num_vehicles * (problem.num_customers + 1);
        let num_variables_inventory = problem.num_days * (problem.num_customers + 1);
        let num_variables_deliver = problem.num_days * problem.num_vehicles * problem.num_customers;
        let num_variables = num_variables_route
            + num_variables_visit
            + num_variables_inventory
            + num_variables_deliver;
        let mut vars = Variables {
            problem,
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
                    for j in 0..=problem.num_customers {
                        if i != j {
                            let name = format!("r_{}_{}_{}_{}", t, v, i, j);
                            let coeff = 0.0;
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
                            gurobi::VarType::Continuous,
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
        let n = self.problem.num_customers;
        debug_assert!(i <= n);
        debug_assert!(j <= n);
        debug_assert!(i != j);

        let j_block_size = n;
        let i_j_block_size = (n + 1) * j_block_size;
        let v_i_j_block_size = self.problem.num_vehicles * i_j_block_size;
        let offset =
            self.route_range.0 + t * v_i_j_block_size + v * i_j_block_size + i * j_block_size;
        let result = offset + (if j > i { j - 1 } else { j });
        debug_assert!(result < self.route_range.1);

        result
    }

    fn route(&self, t: usize, v: usize, i: usize, j: usize) -> gurobi::Var {
        self.variables[self.route_index(t, v, i, j)]
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

enum SolverState {
    StartHeuristic,
    Normal,
}

struct SolverData<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    state: SolverState,
    varnames: Vec<String>,
    start_time: time::Instant,
    cpu: String,
    ncalls: usize,
}

impl<'a> SolverData<'a> {
    const EPSILON: f64 = 1e-7;

    fn new(problem: &'a Problem, lp: &mut gurobi::Model, cpu: String) -> Self {
        let start_time = time::Instant::now();
        let vars = Variables::new(problem, lp);
        lp.update().unwrap(); // update to access variable names
        let varnames = vars
            .variables
            .iter()
            .map(|var| lp.get_obj_attr(grb::attr::VarName, var).unwrap())
            .collect();

        SolverData {
            problem,
            vars,
            state: SolverState::StartHeuristic,
            varnames,
            start_time,
            cpu,
            ncalls: 0,
        }
    }

    fn integral_subtour_elimination<F>(&mut self, assignment: &[f64], add: F) -> grb::Result<bool>
    where
        F: Fn(grb::constr::IneqExpr) -> grb::Result<()>,
    {
        let mut added = false;

        self.varnames
            .iter()
            .zip(assignment.iter())
            .filter(|(_, &value)| value > Self::EPSILON)
            .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));

        // collect node sets of connected components for every day and every vehicle
        let mut sets: Vec<Vec<usize>> = Vec::new();
        for t in 0..self.problem.num_days {
            for v in 0..self.problem.num_vehicles {
                // find connected components with union-find data structure
                let mut uf = partitions::partition_vec![(); self.problem.num_customers + 1];
                for i in 0..=self.problem.num_customers {
                    for j in 0..=self.problem.num_customers {
                        if i != j {
                            let idx = self.vars.route_index(t, v, i, j);
                            if assignment[idx] > 0.5 {
                                uf.union(i, j);
                            }
                        }
                    }
                }

                // collect sets
                let mut component_count = 0;
                for set in uf.all_sets() {
                    let mut set_vec = Vec::new();
                    for (index, _) in set {
                        set_vec.push(index);
                    }
                    if set_vec.len() > 1 && !set_vec.contains(&0) {
                        sets.push(set_vec);
                        component_count += 1;
                    }
                }

                if component_count > 0 {
                    // we have a violated constraint
                    added = true;
                }
            }
        }

        if false {
            for set in sets.iter() {
                eprintln!("# Add node set for all days and vehicles:");
                for i in set.iter() {
                    eprintln!("#  * {}", i);
                }
            }
        }

        // now add all sets (whether violated or not) to all days and vehicles
        if added {
            for t in 0..self.problem.num_days {
                for v in 0..self.problem.num_vehicles {
                    // add subtour elimination constraints if necessary
                    for set in sets.iter() {
                        for k in set.iter() {
                            let mut lhs = grb::expr::LinExpr::new();

                            for i in set.iter() {
                                for j in set.iter() {
                                    if i != j {
                                        lhs.add_term(1.0, self.vars.route(t, v, *i, *j));
                                    }
                                }
                            }
                            for i in set.iter() {
                                if i != k {
                                    lhs.add_term(-1.0, self.vars.visit(t, v, *i));
                                }
                            }

                            add(grb::c!(lhs <= 0))?;
                        }
                    }
                }
            }
        }

        Ok(added)
    }

    fn elapsed_time(&self) -> time::Duration {
        self.start_time.elapsed()
    }

    fn find_next_site(&self, solution: &[f64], t: usize, v: usize, i: usize) -> Option<usize> {
        match self.state {
            SolverState::StartHeuristic => {
                // Simply find the next customer with a delivery
                for j in (i + 1)..=self.problem.num_customers {
                    let var_deliver = self.vars.deliver_index(t, v, j);
                    if solution[var_deliver] > 0.5 {
                        return Some(j);
                    }
                }
                None
            }
            SolverState::Normal => {
                // Find the next site by route variables
                for j in 0..=self.problem.num_customers {
                    if i != j {
                        let var_route = self.vars.route_index(t, v, i, j);
                        if solution[var_route] > 0.5 {
                            return Some(j);
                        }
                    }
                }
                None
            }
        }
    }

    fn get_delivery_amount(&self, solution: &[f64], t: usize, v: usize, target: usize) -> usize {
        if target == 0 {
            0
        } else {
            let var_deliver = self.vars.deliver_index(t, v, target);
            solution[var_deliver].round() as usize
        }
    }

    fn get_routes(&self, solution: &[f64]) -> Routes {
        let mut routes = Vec::with_capacity(self.problem.num_days);

        for t in 0..self.problem.num_days {
            routes.push(Vec::new());

            for v in 0..self.problem.num_vehicles {
                let mut initial_route = vec![Delivery {
                    quantity: 0,
                    customer: 0,
                }];

                let mut i = 0; // last visited site
                while let Some(j) = self.find_next_site(solution, t, v, i) {
                    if j == 0 {
                        break;
                    }

                    let quantity = self.get_delivery_amount(solution, t, v, j);
                    if quantity > 0 {
                        initial_route.push(Delivery {
                            quantity,
                            customer: j,
                        });
                    }
                    // note that results may deviate from intermediate MIP results
                    // because of the "if"

                    i = j
                }

                match self.state {
                    SolverState::StartHeuristic => {
                        let tsp_instance: Vec<(usize, f64, f64)> = initial_route
                            .iter()
                            .map(|delivery| self.problem.id_with_coords(delivery.customer))
                            .collect();
                        // println!("TSP INSTANCE: {:?}", tsp_instance);
                        let mut tsp_tour = lkh::run(&tsp_instance);
                        // println!("Resulting tour:");
                        // for site in tsp_tour.iter() {
                        //     println!(" * {}", site);
                        // }

                        let depot_position = tsp_tour
                            .iter()
                            .position(|site| *site == 0)
                            .expect("Depot is expected to be in TSP tour");
                        tsp_tour.rotate_left(depot_position);
                        let heuristic_route = tsp_tour
                            .iter()
                            .map(|site| Delivery {
                                quantity: self.get_delivery_amount(solution, t, v, *site),
                                customer: *site,
                            })
                            .collect();

                        routes[t].push(heuristic_route);
                    }
                    SolverState::Normal => {
                        routes[t].push(initial_route);
                    }
                }
            }
        }

        routes
    }
}

impl<'a> grb::callback::Callback for SolverData<'a> {
    fn callback(&mut self, w: gurobi::Where) -> grb::callback::CbResult {
        match w {
            gurobi::Where::MIPSol(ctx) => {
                self.ncalls += 1;
                let assignment = ctx.get_solution(&self.vars.variables)?;
                eprintln!("# Incumbent {}!", self.ncalls);
                eprintln!("#    current obj: {}", ctx.obj()?);
                eprintln!("#       best obj: {}", ctx.obj_best()?);

                let new_subtour_constraints =
                    self.integral_subtour_elimination(&assignment, |constr| ctx.add_lazy(constr))?;

                if !new_subtour_constraints {
                    let routes = self.get_routes(&assignment);
                    let solution =
                        Solution::new(self.problem, routes, self.elapsed_time(), self.cpu.clone());

                    eprintln!("{}", solution);
                }
            }
            gurobi::Where::MIPNode(ctx) => {
                let status = ctx.status()?;
                if status == grb::Status::Optimal {
                    eprintln!(
                        "# MIPNode {} {} {:?}",
                        ctx.sol_cnt()?,
                        ctx.node_cnt()?,
                        status
                    );
                    eprintln!("#       best objective: {}", ctx.obj_best()?);
                    eprintln!("#       best obj bound: {}", ctx.obj_bnd()?);
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
        let mut env = gurobi::Env::new("")?;
        env.set(grb::param::Threads, 1)?;

        let mut lp = gurobi::Model::with_env("irp", &env)?;

        lp.set_param(grb::param::LazyConstraints, 1)?;
        lp.set_param(grb::param::Method, 1)?; // use dual simplex
        lp.set_param(grb::param::Sifting, 0)?; // disable sifting

        // 2 = Devex ... partially very good ... needs more experiments on big instances
        //lp.set_param(grb::param::SimplexPricing, 2)?; // use Devex algorithm
        //lp.set_param(grb::param::NormAdjust, 1)?; // ... idk, this was good for devrun
        // maybe with NormAdjust setting?

        // The following parameters do not seem to have a notable effect:
        /*
        lp.set_param(grb::param::Quad, 0)?; // disable Quad precision for Simplex

        // precision used in DIMACS instances is 1e-2, so we use 9e-3 as absolute MIP gap
        lp.set_param(grb::param::MIPGapAbs, 9e-3)?;

        lp.set_param(grb::param::IntegralityFocus, 1)?;
        */

        lp.set_param(grb::param::Presolve, 2)?;

        lp.set_objective(0, gurobi::ModelSense::Minimize)?;

        let mut data = SolverData::new(problem, &mut lp, cpu);

        // route in-degree constraints
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-1.0, data.vars.visit(t, v, i));
                    for j in 0..=problem.num_customers {
                        if i != j {
                            lhs.add_term(1.0, data.vars.route(t, v, j, i));
                        }
                    }

                    lp.add_constr(&format!("Ri_{}_{}_{}", t, v, i), grb::c!(lhs == 0.0))?;
                }
            }
        }

        // route out-degree constraints
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-1.0, data.vars.visit(t, v, i));
                    for j in 0..=problem.num_customers {
                        if i != j {
                            lhs.add_term(1.0, data.vars.route(t, v, i, j));
                        }
                    }

                    lp.add_constr(&format!("Ro_{}_{}_{}", t, v, i), grb::c!(lhs == 0.0))?;
                }
            }
        }

        // visit depot if customer is visited
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 1..=problem.num_customers {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(1.0, data.vars.visit(t, v, 0));
                    lhs.add_term(-1.0, data.vars.visit(t, v, i));

                    lp.add_constr(&format!("DC_{}_{}_{}", t, v, i), grb::c!(lhs >= 0.0))?;
                }
            }
        }

        // at most one visit per day
        for t in 0..problem.num_days {
            for i in 1..=problem.num_customers {
                let mut lhs = grb::expr::LinExpr::new();
                for v in 0..problem.num_vehicles {
                    lhs.add_term(1.0, data.vars.visit(t, v, i));
                }
                lp.add_constr(&format!("V1d_{}_{}", t, i), grb::c!(lhs <= 1.0))?;
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

        // ensure vehicle capacity is not exceeded
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                let mut lhs = grb::expr::LinExpr::new();
                for i in 1..=problem.num_customers {
                    lhs.add_term(1.0, data.vars.deliver(t, v, i));
                }

                lp.add_constr(
                    &format!("VC_{}_{}", t, v),
                    grb::c!(lhs <= problem.capacity as f64),
                )?;
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
        // Get start solution
        lp.optimize()?;
        Self::print_raw_solution(&data, &lp)?;
        let assignment = lp.get_obj_attr_batch(grb::attr::X, data.vars.variables.clone())?;
        let routes = data.get_routes(&assignment);
        let solution = Solution::new(problem, routes, data.elapsed_time(), data.cpu.clone());
        eprintln!("# Start heuristic solution");
        eprintln!("{}", solution);

        data.state = SolverState::Normal;

        // Populate route costs
        for t in 0..problem.num_days {
            for v in 0..problem.num_vehicles {
                for i in 0..=problem.num_customers {
                    for j in 0..=problem.num_customers {
                        if i != j {
                            let var_route = data.vars.route(t, v, i, j);
                            let coeff = problem.distance(i, j).into();
                            lp.set_obj_attr(grb::attr::Obj, &var_route, coeff)?;
                            lp.set_obj_attr(grb::attr::Start, &var_route, 0.0)?;
                        }
                    }
                }
                for deliveries in solution.routes[t][v].windows(2) {
                    let var_route =
                        data.vars
                            .route(t, v, deliveries[0].customer, deliveries[1].customer);
                    lp.set_obj_attr(grb::attr::Start, &var_route, 1.0)?;
                }
            }
        }

        // Get final solution
        lp.optimize_with_callback(&mut data)?;
        Self::print_raw_solution(&data, &lp)?;
        let assignment = lp.get_obj_attr_batch(grb::attr::X, data.vars.variables.clone())?;
        let routes = data.get_routes(&assignment);
        let solution = Solution::new(problem, routes, data.elapsed_time(), data.cpu);
        eprintln!("# Final solution");
        eprintln!("{}", solution);

        Ok(())
    }

    fn print_raw_solution(data: &SolverData, lp: &gurobi::Model) -> grb::Result<()> {
        let status = lp.status()?;
        eprintln!("# MIP solution status: {:?}", status);

        let objective = lp.get_attr(gurobi::attr::ObjVal)?;
        eprintln!("# MIP solution value: {}", objective);

        // print raw solution:
        for var in data.vars.variables.iter() {
            let name = lp.get_obj_attr(grb::attr::VarName, var)?;
            let value = lp.get_obj_attr(grb::attr::X, var)?;
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
