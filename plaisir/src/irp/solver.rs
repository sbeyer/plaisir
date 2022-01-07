use super::*;
use grb::prelude as gurobi;
use std::cmp::Ordering;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time;

extern crate partitions;

const PRINT_VARIABLE_VALUES: bool = false;
const PRINT_ELIMINATED_SUBTOURS: bool = false;

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
    fn new(problem: &'a Problem, lp: &mut gurobi::Model) -> grb::Result<Self> {
        let num_variables_route =
            problem.num_days * problem.num_vehicles * problem.num_sites * problem.num_customers;
        let num_variables_visit = problem.num_days * problem.num_vehicles * problem.num_sites;
        let num_variables_inventory = problem.num_days * problem.num_sites;
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
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_sites() {
                    for j in problem.all_sites_except(i) {
                        let name = format!("r_{}_{}_{}_{}", t, v, i, j);
                        let coeff = problem.distance(i, j).into();
                        let bounds = (0.0, 1.0);
                        //let var = grb::add_binvar!(lp, name: &name, obj: coeff)?;
                        let var = lp.add_var(
                            &name,
                            gurobi::VarType::Binary,
                            coeff,
                            bounds.0,
                            bounds.1,
                            std::iter::empty(),
                        )?;
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
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_sites() {
                    let name = format!("v_{}_{}_{}", t, v, i);
                    let coeff = 0.0;
                    let bounds = (0.0, 1.0);
                    let var = lp.add_var(
                        &name,
                        gurobi::VarType::Continuous,
                        coeff,
                        bounds.0,
                        bounds.1,
                        std::iter::empty(),
                    )?;
                    debug_assert_eq!(vars.variables.len(), vars.visit_index(t, v, i));
                    lp.set_obj_attr(grb::attr::BranchPriority, &var, 1).unwrap();
                    vars.variables.push(var);
                }
            }
        }
        vars.visit_range.1 = vars.variables.len();
        debug_assert_eq!(num_variables_visit, vars.visit_range.1 - vars.visit_range.0);

        // inventory variables
        vars.inventory_range.0 = vars.visit_range.1;
        for t in problem.all_days() {
            for i in problem.all_sites() {
                let site = problem.site(i);
                let name = format!("i_{}_{}", t, i);
                let coeff = site.cost();
                let bounds = site.level_bounds();
                let var = lp.add_var(
                    &name,
                    gurobi::VarType::Continuous,
                    coeff,
                    bounds.0,
                    bounds.1 + site.level_change(),
                    std::iter::empty(),
                )?;
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
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_customers() {
                    let name = format!("d_{}_{}_{}", t, v, i);
                    let coeff = 0.0;
                    let bounds = (0.0, problem.capacity as f64);
                    let var = lp.add_var(
                        &name,
                        gurobi::VarType::Continuous,
                        coeff,
                        bounds.0,
                        bounds.1,
                        std::iter::empty(),
                    )?;
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

        Ok(vars)
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

        let i_block_size = self.problem.num_sites;
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

        let i_block_size = self.problem.num_sites;
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

#[derive(Eq)]
struct Delivery {
    quantity: usize,
    customer: usize,
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
    fn empty() -> Self {
        Self {
            routes: Vec::new(),
            cost_transportation: f64::INFINITY,
            cost_inventory_depot: f64::INFINITY,
            cost_inventory_customers: f64::INFINITY,
            cost_total: f64::INFINITY,
            processor: String::new(),
            time: time::Duration::default(),
        }
    }

    fn new(problem: &Problem, routes: Routes, time: time::Duration, cpu: &str) -> Self {
        let mut sol = Solution {
            routes,
            cost_transportation: 0.,
            cost_inventory_depot: 0.,
            cost_inventory_customers: 0.,
            cost_total: 0.,
            processor: cpu.to_string(),
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
        let mut inventory: Vec<f64> = problem
            .all_sites()
            .map(|i| problem.site(i).level_start())
            .collect();
        let mut cost_inventory = vec![0.; problem.num_sites];

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
            for i in problem.all_sites() {
                inventory[i] += problem.site(i).level_change();
            }

            // update inventory costs
            for i in problem.all_sites() {
                cost_inventory[i] += problem.site(i).cost() * inventory[i]
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

struct McfSubproblem {
    model: gurobi::Model,
    delivery_vars: Vec<Vec<Vec<gurobi::Var>>>,
}

impl McfSubproblem {
    fn new(env: &gurobi::Env, problem: &Problem) -> grb::Result<Self> {
        let mut model = gurobi::Model::with_env("deliveries", env)?;

        model.set_objective(0, gurobi::ModelSense::Minimize)?;

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
                                gurobi::VarType::Continuous,
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
                                        gurobi::VarType::Continuous,
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

    fn delivery_var(&self, t: usize, v: usize, i: usize) -> grb::Var {
        self.delivery_vars[t][v][i - 1]
    }
}

struct SolverData<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    varnames: Vec<String>,
    start_time: time::Instant,
    cpu: &'a str,
    ncalls: usize,
    mcf: McfSubproblem,
    best_solution: Solution,
    is_new_solution_just_set: bool,
}

impl<'a> SolverData<'a> {
    const EPSILON: f64 = 1e-7;

    fn new(
        problem: &'a Problem,
        lp: &mut gurobi::Model,
        env: &'a gurobi::Env,
        cpu: &'a str,
    ) -> grb::Result<Self> {
        let start_time = time::Instant::now();
        let vars = Variables::new(problem, lp)?;
        lp.update()?; // update to access variable names
        let varnames = vars
            .variables
            .iter()
            .map(|var| lp.get_obj_attr(grb::attr::VarName, var).unwrap())
            .collect();

        let mcf = McfSubproblem::new(env, problem)?;

        let best_solution = Solution::empty();

        Ok(SolverData {
            problem,
            vars,
            varnames,
            start_time,
            cpu,
            ncalls: 0,
            mcf,
            best_solution,
            is_new_solution_just_set: false,
        })
    }

    fn integral_subtour_elimination<F>(&mut self, assignment: &[f64], add: F) -> grb::Result<bool>
    where
        F: Fn(grb::constr::IneqExpr) -> grb::Result<()>,
    {
        let mut added = false;

        if PRINT_VARIABLE_VALUES {
            self.varnames
                .iter()
                .zip(assignment.iter())
                .filter(|(_, &value)| value > Self::EPSILON)
                .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));
        }

        // collect node sets of connected components for every day and every vehicle
        let mut sets: Vec<Vec<usize>> = Vec::new();
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                // find connected components with union-find data structure
                let mut uf = partitions::partition_vec![(); self.problem.num_sites];
                for i in self.problem.all_sites() {
                    for j in self.problem.all_sites_except(i) {
                        let idx = self.vars.route_index(t, v, i, j);
                        if assignment[idx] > 0.5 {
                            uf.union(i, j);
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

        if PRINT_ELIMINATED_SUBTOURS {
            for set in sets.iter() {
                eprintln!("# Add node set for all days and vehicles:");
                for i in set.iter() {
                    eprintln!("#  * {}", i);
                }
            }
        }

        // now add all sets (whether violated or not) to all days and vehicles
        if added {
            for t in self.problem.all_days() {
                for v in self.problem.all_vehicles() {
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

    fn get_best_solution_variable_assignment(&self) -> Vec<f64> {
        // Inventory levels, necessary for inventory variables
        let mut levels = self
            .problem
            .all_sites()
            .map(|i| self.problem.site(i).level_start() as isize)
            .collect::<Vec<_>>();

        let mut assignment = vec![0.0; self.vars.variables.len()];
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                let route = &self.best_solution.routes[t][v];
                if route.len() > 1 {
                    for delivery in route.iter() {
                        let i = delivery.customer;
                        let visit_index = self.vars.visit_index(t, v, i);
                        assignment[visit_index] = 1.0;
                        if i > 0 {
                            let deliver_index = self.vars.deliver_index(t, v, i);
                            assignment[deliver_index] = delivery.quantity as f64;

                            levels[0] -= delivery.quantity as isize;
                            levels[i] += delivery.quantity as isize;
                        }
                    }
                    for deliveries in route.windows(2) {
                        let i = deliveries[0].customer;
                        let j = deliveries[1].customer;
                        let index = self.vars.route_index(t, v, i, j);
                        assignment[index] = 1.0;
                    }
                    // Routes in our solutions do not end at the depot, but we want to
                    // go back to the depot at the end in the MIP solution.
                    {
                        let i = route[route.len() - 1].customer;
                        let j = 0;
                        let index = self.vars.route_index(t, v, i, j);
                        assignment[index] = 1.0;
                    }
                }
            }

            // Set inventory variable values
            #[allow(clippy::needless_range_loop)]
            for i in self.problem.all_sites() {
                levels[i] += self.problem.site(i).level_change() as isize;
                let index = self.vars.inventory_index(t, i);
                assignment[index] = levels[i] as f64;
            }
        }

        assignment
    }

    fn update_best_solution(&mut self, routes: Routes) {
        let solution = Solution::new(self.problem, routes, self.elapsed_time(), self.cpu);

        if solution.cost_total < self.best_solution.cost_total {
            eprintln!("{}", solution);
            self.best_solution = solution;
        } else {
            eprintln!(
                "# Found solution of objective value {} not better than {}",
                solution.cost_total, self.best_solution.cost_total
            );
        }
    }

    fn give_new_best_solution_to_solver(
        &mut self,
        ctx: grb::callback::MIPNodeCtx,
    ) -> grb::Result<()> {
        let best_objective = ctx.obj_best()?;
        if self.best_solution.cost_total < best_objective {
            let best_solution_assignment = self.get_best_solution_variable_assignment();
            let set_result =
                ctx.set_solution(self.vars.variables.iter().zip(best_solution_assignment))?;
            if let Some(value) = set_result {
                self.is_new_solution_just_set = true;
                eprintln!(
                    "# New best solution with objective value {} (old: {}) set successfully",
                    value, best_objective
                );
                if value != self.best_solution.cost_total {
                    eprintln!(
                        "# The new objective value deviates from expected value {}",
                        self.best_solution.cost_total
                    );
                }
            } else {
                eprintln!(
                    "# No new solution set, keeping best objective value  {}",
                    best_objective
                );
            }
        }

        Ok(())
    }

    fn elapsed_time(&self) -> time::Duration {
        self.start_time.elapsed()
    }

    fn find_next_site_in_route(
        &self,
        solution: &[f64],
        t: usize,
        v: usize,
        i: usize,
    ) -> Option<usize> {
        // Find the next site by route variables
        for j in self.problem.all_sites_except(i) {
            let var_route = self.vars.route_index(t, v, i, j);
            if solution[var_route] > 0.5 {
                return Some(j);
            }
        }
        None
    }

    fn get_visited_customers(&self, solution: &[f64], t: usize, v: usize) -> Vec<usize> {
        self.problem
            .all_customers()
            .filter(|i| {
                let var_deliver = self.vars.deliver_index(t, v, *i);
                solution[var_deliver] > 0.5
            })
            .collect()
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
        self.problem
            .all_days()
            .map(|t| {
                self.problem
                    .all_vehicles()
                    .map(|v| {
                        let mut route = vec![Delivery {
                            quantity: 0,
                            customer: 0,
                        }];

                        let mut i = 0; // last visited site
                        while let Some(j) = self.find_next_site_in_route(solution, t, v, i) {
                            if j == 0 {
                                break;
                            }

                            let quantity = self.get_delivery_amount(solution, t, v, j);
                            if quantity > 0 {
                                route.push(Delivery {
                                    quantity,
                                    customer: j,
                                });
                            }
                            // note that results may deviate from intermediate MIP results
                            // because of the "if"

                            i = j
                        }

                        route
                    })
                    .collect()
            })
            .collect()
    }

    /// Solve Minimum-Cost Flow LP to improve deliveries
    fn adjust_deliveries(&mut self, solution: &mut [f64]) -> grb::Result<bool> {
        // Update bounds of deliveries
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                for i in self.get_visited_customers(solution, t, v).into_iter() {
                    self.mcf.model.set_obj_attr(
                        grb::attr::UB,
                        &self.mcf.delivery_var(t, v, i),
                        self.problem.capacity as f64,
                    )?;
                }
            }
        }

        self.mcf.model.optimize()?;

        // Reset bounds of deliveries
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                for i in self.get_visited_customers(solution, t, v).into_iter() {
                    self.mcf.model.set_obj_attr(
                        grb::attr::UB,
                        &self.mcf.delivery_var(t, v, i),
                        0.0,
                    )?;
                }
            }
        }

        let status = self.mcf.model.status()?;
        if status == grb::Status::Optimal {
            if true {
                eprintln!("#### MIP solution status: {:?}", status);

                let objective = self.mcf.model.get_attr(gurobi::attr::ObjVal)?;
                eprintln!("#### MIP solution value: {}", objective);

                for delivery_vars_per_day in self.mcf.delivery_vars.iter() {
                    for delivery_vars_per_customer in delivery_vars_per_day.iter() {
                        for var in delivery_vars_per_customer.iter() {
                            let name = self.mcf.model.get_obj_attr(grb::attr::VarName, var)?;
                            let value = self.mcf.model.get_obj_attr(grb::attr::X, var)?;
                            if value > SolverData::EPSILON {
                                eprintln!("####   - {}: {}", name, value);
                            }
                        }
                    }
                }
            }

            if status == grb::Status::Optimal {
                for t in self.problem.all_days() {
                    for v in self.problem.all_vehicles() {
                        for i in self.problem.all_customers() {
                            let var = self.mcf.delivery_var(t, v, i);
                            let value = self.mcf.model.get_obj_attr(grb::attr::X, &var)?;
                            solution[self.vars.deliver_index(t, v, i)] = value;
                        }
                    }
                }
            }

            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// A really stupid heuristic that does not take anything into consideration that would be sane
    ///
    /// Note that this method just changes the delivery values although it would also have
    /// to change other values to get a feasible solution. We are assuming that all following
    /// methods ignore the parts that are not related to deliveries.
    fn fix_deliveries(&mut self, solution: &mut [f64]) -> grb::Result<bool> {
        for t in self.problem.all_days() {
            eprintln!("# Fixing day {}", t);

            // Compute load for each vehicle
            let load = self
                .problem
                .all_vehicles()
                .map(|v| {
                    self.problem
                        .all_customers()
                        .map(|i| self.get_delivery_amount(solution, t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!(
                "## Loads {}: {:?} capacity {}",
                t, load, self.problem.capacity
            );

            // Move delivieries to the same customer on different routes to the first route
            // with such a delivery. (Note that this does not change anything for customers
            // that are visited at most once.)
            for i in self.problem.all_customers() {
                let mut total_delivery: usize = self
                    .problem
                    .all_vehicles()
                    .map(|v| self.get_delivery_amount(solution, t, v, i))
                    .sum();

                for v in self.problem.all_vehicles() {
                    let value = self.get_delivery_amount(solution, t, v, i);
                    if value > 0 {
                        solution[self.vars.deliver_index(t, v, i)] = total_delivery as f64;
                        total_delivery = 0;
                    }
                }
            }

            // Compute load for each vehicle
            let mut load = self
                .problem
                .all_vehicles()
                .map(|v| {
                    self.problem
                        .all_customers()
                        .map(|i| self.get_delivery_amount(solution, t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!(
                "## Loads {}: {:?} capacity {}",
                t, load, self.problem.capacity
            );

            // Sort vehicles by descending load
            let mut sorted_vehicles = self.problem.all_vehicles().collect::<Vec<_>>();
            sorted_vehicles.sort_by_cached_key(|a| Reverse(load[*a]));

            // Move the smallest deliveries that overload the capacity to the next route
            for (i, v_ref) in sorted_vehicles[..sorted_vehicles.len() - 1]
                .iter()
                .enumerate()
            {
                let v = *v_ref;
                let v_next = sorted_vehicles[i + 1];
                if load[v] > self.problem.capacity {
                    let mut heap = self
                        .get_visited_customers(solution, t, v)
                        .into_iter()
                        .map(|i| {
                            Reverse(Delivery {
                                quantity: self.get_delivery_amount(solution, t, v, i),
                                customer: i,
                            })
                        })
                        .collect::<BinaryHeap<_>>();

                    while load[v] > self.problem.capacity {
                        if let Some(delivery) = heap.pop() {
                            let amount = delivery.0.quantity;
                            load[v] -= amount;
                            load[v_next] += amount;

                            let i = delivery.0.customer;
                            solution[self.vars.deliver_index(t, v, i)] = 0.0;
                            solution[self.vars.deliver_index(t, v_next, i)] = amount as f64;
                        } else {
                            panic!("Logic error: as long as the load is positive, we must have elements on the heap")
                        };
                    }
                }
            }

            // Compute load for each vehicle
            let load = self
                .problem
                .all_vehicles()
                .map(|v| {
                    self.problem
                        .all_customers()
                        .map(|i| self.get_delivery_amount(solution, t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!(
                "## Loads {}: {:?} capacity {}",
                t, load, self.problem.capacity
            );

            if load[sorted_vehicles[sorted_vehicles.len() - 1]] > self.problem.capacity {
                return Ok(false);
            }
        }

        Ok(true)
    }

    /// Runs LKH heuristic on visited sites to get a feasible route
    fn get_routes_heuristically(&self, solution: &[f64]) -> Routes {
        self.problem
            .all_days()
            .map(|t| {
                self.problem
                    .all_vehicles()
                    .map(|v| {
                        let mut visited_sites = self.get_visited_customers(solution, t, v);
                        visited_sites.push(0); // add depot

                        let tsp_instance: Vec<(usize, f64, f64)> = visited_sites
                            .iter()
                            .map(|site| {
                                let site = self.problem.site(*site);
                                let pos = &site.position();
                                (site.id(), pos.x, pos.y)
                            })
                            .collect();
                        let mut tsp_tour = lkh::run(&tsp_instance);

                        let depot_position = tsp_tour
                            .iter()
                            .position(|site| *site == 0)
                            .expect("Depot is expected to be in TSP tour");
                        tsp_tour.rotate_left(depot_position);

                        tsp_tour
                            .iter()
                            .map(|site| Delivery {
                                quantity: self.get_delivery_amount(solution, t, v, *site),
                                customer: *site,
                            })
                            .collect()
                    })
                    .collect()
            })
            .collect()
    }
}

impl<'a> grb::callback::Callback for SolverData<'a> {
    fn callback(&mut self, w: gurobi::Where) -> grb::callback::CbResult {
        match w {
            gurobi::Where::MIPSol(ctx) => {
                self.ncalls += 1;
                let mut assignment = ctx.get_solution(&self.vars.variables)?;
                eprintln!("# Incumbent {}!", self.ncalls);
                eprintln!("#    current obj: {}", ctx.obj()?);
                eprintln!("#       best obj: {}", ctx.obj_best()?);
                eprintln!(
                    "#           time: {}",
                    self.elapsed_time().as_millis() as f64 * 1e-3
                );

                if self.is_new_solution_just_set {
                    self.is_new_solution_just_set = false;
                } else {
                    {
                        if PRINT_VARIABLE_VALUES {
                            self.varnames
                                .iter()
                                .zip(assignment.iter())
                                .filter(|(_, &value)| value > Self::EPSILON)
                                .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));
                        }

                        let adjusted = self.adjust_deliveries(&mut assignment)?;
                        if adjusted {
                            let routes = self.get_routes_heuristically(&assignment);
                            self.update_best_solution(routes);
                        } else {
                            eprintln!(
                                "# Failed to find a feasible adjusted solution. \
                                   This is weird, because we are in MIPSol."
                            );
                        }
                    }

                    let new_subtour_constraints = self
                        .integral_subtour_elimination(&assignment, |constr| ctx.add_lazy(constr))?;

                    if !new_subtour_constraints {
                        let routes = self.get_routes(&assignment);
                        self.update_best_solution(routes);
                    }

                    eprintln!(
                        "# Callback finish time: {}",
                        self.elapsed_time().as_millis() as f64 * 1e-3
                    );
                }
            }
            gurobi::Where::MIPNode(ctx) => {
                let status = ctx.status()?;
                if status == grb::Status::Optimal {
                    eprintln!(
                        "# MIPNode #sols {} #nodes {} status {:?}",
                        ctx.sol_cnt()?,
                        ctx.node_cnt()?,
                        status,
                    );
                    eprintln!("#       best objective: {}", ctx.obj_best()?);
                    eprintln!("#       best obj bound: {}", ctx.obj_bnd()?);
                    eprintln!(
                        "#                 time: {}",
                        self.elapsed_time().as_millis() as f64 * 1e-3
                    );

                    if true {
                        let mut assignment = ctx.get_solution(&self.vars.variables)?;
                        if PRINT_VARIABLE_VALUES {
                            self.varnames
                                .iter()
                                .zip(assignment.iter())
                                .filter(|(_, &value)| value > Self::EPSILON)
                                .for_each(|(var, value)| eprintln!("#   - {}: {}", var, value));
                        }
                        let fixed = self.adjust_deliveries(&mut assignment)?
                            && self.fix_deliveries(&mut assignment)?;
                        if fixed {
                            let routes = self.get_routes_heuristically(&assignment);
                            self.update_best_solution(routes);
                        } else {
                            eprintln!("# Failed to find a feasible solution.");
                        }
                    }

                    self.give_new_best_solution_to_solver(ctx)?;

                    eprintln!(
                        "# MIPNode finish time: {}",
                        self.elapsed_time().as_millis() as f64 * 1e-3
                    );
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

        let mut data = SolverData::new(problem, &mut lp, &env, &cpu)?;

        // route in-degree constraints
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_sites() {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-1.0, data.vars.visit(t, v, i));
                    for j in problem.all_sites_except(i) {
                        lhs.add_term(1.0, data.vars.route(t, v, j, i));
                    }

                    lp.add_constr(&format!("Ri_{}_{}_{}", t, v, i), grb::c!(lhs == 0.0))?;
                }
            }
        }

        // route out-degree constraints
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_sites() {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-1.0, data.vars.visit(t, v, i));
                    for j in problem.all_sites_except(i) {
                        lhs.add_term(1.0, data.vars.route(t, v, i, j));
                    }

                    lp.add_constr(&format!("Ro_{}_{}_{}", t, v, i), grb::c!(lhs == 0.0))?;
                }
            }
        }

        // visit depot if customer is visited
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_customers() {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(1.0, data.vars.visit(t, v, 0));
                    lhs.add_term(-1.0, data.vars.visit(t, v, i));

                    lp.add_constr(&format!("DC_{}_{}_{}", t, v, i), grb::c!(lhs >= 0.0))?;
                }
            }
        }

        // at most one visit per day
        for t in problem.all_days() {
            for i in problem.all_customers() {
                let mut lhs = grb::expr::LinExpr::new();
                for v in problem.all_vehicles() {
                    lhs.add_term(1.0, data.vars.visit(t, v, i));
                }
                lp.add_constr(&format!("V1d_{}_{}", t, i), grb::c!(lhs <= 1.0))?;
            }
        }

        // glue: disable delivery if we do not visit
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_customers() {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(problem.capacity as f64, data.vars.visit(t, v, i));
                    lhs.add_term(-1.0, data.vars.deliver(t, v, i));

                    lp.add_constr(&format!("Gdv_{}_{}_{}", t, v, i), grb::c!(lhs >= 0.0))?;
                }
            }
        }

        // ensure vehicle capacity is not exceeded
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                let mut lhs = grb::expr::LinExpr::new();
                for i in problem.all_customers() {
                    lhs.add_term(1.0, data.vars.deliver(t, v, i));
                }

                lp.add_constr(
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
                    lhs.add_term(-1.0, data.vars.deliver(t, v, j));
                }
            }
            lhs.add_term(-1.0, data.vars.inventory(t, 0)); // outgoing inventory
            let depot = problem.site(0);
            let mut value = -depot.level_change();

            if t == 0 {
                value -= depot.level_start();
            } else {
                lhs.add_term(1.0, data.vars.inventory(t - 1, 0)); // incoming inventory
            }

            lp.add_constr(&format!("Ifd_{}", t), grb::c!(lhs == value))?;
        }

        // inventory flow for customers
        for t in problem.all_days() {
            for i in problem.all_customers() {
                let mut lhs = grb::expr::LinExpr::new();
                for v in problem.all_vehicles() {
                    lhs.add_term(1.0, data.vars.deliver(t, v, i)); // delivered
                }
                lhs.add_term(-1.0, data.vars.inventory(t, i)); // outgoing inventory
                let customer = problem.site(i);
                let mut value = -customer.level_change();

                if t == 0 {
                    value -= customer.level_start();
                } else {
                    lhs.add_term(1.0, data.vars.inventory(t - 1, i)); // incoming inventory
                }

                lp.add_constr(&format!("Ifc_{}_{}", t, i), grb::c!(lhs == value))?;
            }
        }

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

        if PRINT_VARIABLE_VALUES {
            for var in data.vars.variables.iter() {
                let name = lp.get_obj_attr(grb::attr::VarName, var)?;
                let value = lp.get_obj_attr(grb::attr::X, var)?;
                if value > SolverData::EPSILON {
                    eprintln!("#   - {}: {}", name, value);
                }
            }
        }

        Ok(())
    }
}

pub fn solve(problem: Problem, cpu: String) {
    Solver::solve(&problem, cpu).unwrap();
}
