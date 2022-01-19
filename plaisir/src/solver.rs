use crate::delivery::Solver as DeliverySolver;
use crate::delivery::{Deliveries, Delivery};
use crate::heuristic::*;
use crate::problem::*;
use crate::solution::*;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::time;

extern crate partitions;

const PRINT_VARIABLE_VALUES: bool = false;
const PRINT_ELIMINATED_SUBTOURS: bool = false;

struct Variables<'a> {
    problem: &'a Problem,
    variables: Vec<grb::Var>,
    route_range: (usize, usize),
    deliver_range: (usize, usize),
    visit_range: (usize, usize),
    inventory_range: (usize, usize),
}

#[allow(clippy::many_single_char_names)]
impl<'a> Variables<'a> {
    fn new(problem: &'a Problem, lp: &mut grb::Model) -> grb::Result<Self> {
        let num_variables_route = problem.num_days
            * problem.num_vehicles
            * (problem.num_sites * problem.num_customers / 2);
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
                    for j in problem.all_sites_after(i) {
                        let name = format!("r_{t}_{v}_{i}_{j}");
                        let coeff = problem.distance(i, j).into();
                        let bounds = (0.0, if i == 0 { 2.0 } else { 1.0 });
                        //let var = grb::add_binvar!(lp, name: &name, obj: coeff)?;
                        let var = lp.add_var(
                            &name,
                            if i == 0 {
                                grb::VarType::Integer
                            } else {
                                grb::VarType::Binary
                            },
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
                    let name = format!("v_{t}_{v}_{i}");
                    let coeff = 0.0;
                    let bounds = (0.0, 1.0);
                    let var = lp.add_var(
                        &name,
                        grb::VarType::Binary,
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
                let name = format!("i_{t}_{i}");
                let coeff = site.cost();
                let bounds = site.level_bounds();
                let var = lp.add_var(
                    &name,
                    grb::VarType::Continuous,
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
                    let name = format!("d_{t}_{v}_{i}");
                    let coeff = 0.0;
                    let bounds = (0.0, problem.capacity as f64);
                    let var = lp.add_var(
                        &name,
                        grb::VarType::Continuous,
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

    fn route_index(&self, t: DayId, v: VehicleId, i: SiteId, j: SiteId) -> usize {
        let t = t as usize;
        let v = v as usize;
        let i = i as usize;
        let j = j as usize;

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

    fn route_index_undirected(&self, t: DayId, v: VehicleId, i: SiteId, j: SiteId) -> usize {
        if i < j {
            self.route_index(t, v, i, j)
        } else {
            self.route_index(t, v, j, i)
        }
    }

    fn route(&self, t: DayId, v: VehicleId, i: SiteId, j: SiteId) -> grb::Var {
        let index = self.route_index_undirected(t, v, i, j);
        self.variables[index]
    }

    fn visit_index(&self, t: DayId, v: VehicleId, i: SiteId) -> usize {
        let t = t as usize;
        let v = v as usize;
        let i = i as usize;

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

    fn visit(&self, t: DayId, v: VehicleId, i: SiteId) -> grb::Var {
        self.variables[self.visit_index(t, v, i)]
    }

    fn inventory_index(&self, t: DayId, i: SiteId) -> usize {
        let t = t as usize;
        let i = i as usize;

        debug_assert!(t < self.problem.num_days);
        debug_assert!(i <= self.problem.num_customers);

        let i_block_size = self.problem.num_sites;
        let offset = self.inventory_range.0 + t * i_block_size;
        let result = offset + i;
        debug_assert!(result < self.inventory_range.1);

        result
    }

    fn inventory(&self, t: DayId, i: SiteId) -> grb::Var {
        self.variables[self.inventory_index(t, i)]
    }

    fn deliver_index(&self, t: DayId, v: VehicleId, i: SiteId) -> usize {
        let t = t as usize;
        let v = v as usize;
        let i = i as usize;

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

    fn deliver(&self, t: DayId, v: VehicleId, i: SiteId) -> grb::Var {
        self.variables[self.deliver_index(t, v, i)]
    }
}

struct SolverData<'a> {
    problem: &'a Problem,
    vars: Variables<'a>,
    varnames: Vec<String>,
    start_time: time::Instant,
    cpu: &'a str,
    ncalls: usize,
    deliveries: DeliverySolver<'a>,
    best_solution: Solution,
    is_new_solution_just_set: bool,
    heuristic: RandomHeuristic<'a>,
}

impl<'a> SolverData<'a> {
    const EPSILON: f64 = 1e-7;

    fn new(
        problem: &'a Problem,
        lp: &mut grb::Model,
        env: &'a grb::Env,
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

        let deliveries = DeliverySolver::new(env, problem)?;

        let best_solution = Solution::empty();

        let heuristic = RandomHeuristic::new(problem);

        Ok(SolverData {
            problem,
            vars,
            varnames,
            start_time,
            cpu,
            ncalls: 0,
            deliveries,
            best_solution,
            is_new_solution_just_set: false,
            heuristic,
        })
    }

    fn integral_subtour_elimination<F>(&mut self, assignment: &[f64], add: F) -> grb::Result<bool>
    where
        F: Fn(grb::constr::IneqExpr) -> grb::Result<()>,
    {
        let mut added = false;

        eprintln!("# Subtour elimination run, time {}", self.elapsed_seconds());

        if PRINT_VARIABLE_VALUES {
            self.varnames
                .iter()
                .zip(assignment.iter())
                .filter(|(_, &value)| value > Self::EPSILON)
                .for_each(|(var, value)| eprintln!("#   - {var}: {value}"));
        }

        // collect node sets of connected components for every day and every vehicle
        let mut sets: Vec<Vec<SiteId>> = Vec::new();
        for t in self.problem.all_days() {
            for v in self.problem.all_vehicles() {
                // find connected components with union-find data structure
                let mut uf = partitions::partition_vec![(); self.problem.num_sites];
                for i in self.problem.all_sites() {
                    for j in self.problem.all_sites_after(i) {
                        let idx = self.vars.route_index(t, v, i, j);
                        if assignment[idx] > 0.5 {
                            uf.union(i as usize, j as usize);
                        }
                    }
                }

                // collect sets
                let mut component_count = 0;
                for set in uf.all_sets() {
                    let mut set_vec = Vec::new();
                    for (index, _) in set {
                        set_vec.push(index as SiteId);
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
                    eprintln!("#  * {i}");
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
                                    if i < j {
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

        eprintln!(
            "# Subtour elimination -> constraints added? {}, time {}",
            added,
            self.elapsed_seconds()
        );

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
                let route = &self.best_solution.route(t, v);
                if route.len() > 1 {
                    for delivery in route.iter() {
                        let i = delivery.customer;
                        let visit_index = self.vars.visit_index(t, v, i);
                        assignment[visit_index] = 1.0;
                        if i > 0 {
                            let deliver_index = self.vars.deliver_index(t, v, i);
                            assignment[deliver_index] = delivery.quantity as f64;

                            levels[0] -= delivery.quantity as isize;
                            levels[i as usize] += delivery.quantity as isize;
                        }
                    }
                    for deliveries in route.windows(2) {
                        let i = deliveries[0].customer;
                        let j = deliveries[1].customer;
                        let index = self.vars.route_index_undirected(t, v, i, j);
                        assignment[index] = 1.0;
                    }
                    // Routes in our solutions do not end at the depot, but we want to
                    // go back to the depot at the end in the MIP solution.
                    {
                        let i = route[route.len() - 1].customer;
                        let j = 0;
                        let index = self.vars.route_index_undirected(t, v, i, j);
                        assignment[index] += 1.0;
                    }
                }
            }

            // Set inventory variable values
            for i in self.problem.all_sites() {
                levels[i as usize] += self.problem.site(i).level_change() as isize;
                let index = self.vars.inventory_index(t, i);
                assignment[index] = levels[i as usize] as f64;
            }
        }

        assignment
    }

    fn update_best_solution(&mut self, schedule: Schedule) {
        let solution = Solution::new(self.problem, schedule, self.elapsed_seconds(), self.cpu);

        if solution.value() < self.best_solution.value() {
            eprintln!("{}", solution);
            self.best_solution = solution;
        } else {
            eprintln!(
                "# Found solution of objective value {} not better than {}",
                solution.value(),
                self.best_solution.value()
            );
        }
    }

    fn give_new_best_solution_to_solver(
        &mut self,
        ctx: grb::callback::MIPNodeCtx,
    ) -> grb::Result<()> {
        let best_objective = ctx.obj_best()?;
        if self.best_solution.value() < best_objective {
            let best_solution_assignment = self.get_best_solution_variable_assignment();
            let set_result =
                ctx.set_solution(self.vars.variables.iter().zip(best_solution_assignment))?;
            if set_result.is_some() {
                self.is_new_solution_just_set = true;
                eprintln!(
                    "# New best solution with objective value {} (old: {}) set successfully",
                    self.best_solution.value(),
                    best_objective
                );
            } else {
                eprintln!("# No new solution set, keeping best objective value {best_objective}");
            }
        }

        Ok(())
    }

    fn elapsed_seconds(&self) -> f64 {
        self.start_time.elapsed().as_millis() as f64 * 1e-3
    }

    fn is_edge_in_route(
        &self,
        assignment: &[f64],
        t: DayId,
        v: VehicleId,
        i: SiteId,
        j: SiteId,
    ) -> bool {
        let var_route = self.vars.route_index_undirected(t, v, i, j);
        assignment[var_route].round() > Self::EPSILON
    }

    fn get_delivery_amount(
        &self,
        assignment: &[f64],
        t: DayId,
        v: VehicleId,
        target: SiteId,
    ) -> usize {
        if target == 0 {
            0
        } else {
            let var_deliver = self.vars.deliver_index(t, v, target);
            assignment[var_deliver].round() as usize
        }
    }

    fn get_schedule(&self, assignment: &[f64]) -> Schedule {
        eprintln!("# Get schedule, time {}", self.elapsed_seconds());

        Schedule(
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

                            let mut adjacencies =
                                vec![Vec::with_capacity(2); self.problem.num_sites];
                            for i in self.problem.all_sites() {
                                for j in self.problem.all_sites_after(i) {
                                    if self.is_edge_in_route(assignment, t, v, i, j) {
                                        adjacencies[i as usize].push(j);
                                        adjacencies[j as usize].push(i);
                                    }
                                }
                            }

                            let mut visited = vec![false; self.problem.num_sites];
                            visited[0] = true;

                            let mut i = 0; // last visited site
                            loop {
                                let mut found = false;

                                for j in adjacencies[i].iter() {
                                    if !visited[*j as usize] {
                                        let quantity =
                                            self.get_delivery_amount(assignment, t, v, *j);

                                        if quantity > 0 && *j != 0 {
                                            route.push(Delivery {
                                                quantity,
                                                customer: *j,
                                            });
                                        }

                                        visited[*j as usize] = true;
                                        found = true;
                                        i = *j as usize;

                                        break;
                                    }
                                }

                                if !found || i == 0 {
                                    break;
                                }
                            }

                            route
                        })
                        .collect()
                })
                .collect(),
        )
    }

    /// Solve Minimum-Cost Flow LP to improve deliveries (based on currently visited customers)
    fn adjust_deliveries(&mut self, assignment: &[f64]) -> grb::Result<Option<Deliveries>> {
        eprintln!("# Adjust deliveries, time {}", self.elapsed_seconds());

        self.deliveries.set_all_statuses(|t, v, i| {
            let var_deliver = self.vars.deliver_index(t, v, i);
            assignment[var_deliver] > 0.5
        })?;
        let result = self.deliveries.solve()?;

        Ok(result)
    }

    fn fractional_delivery_heuristic(
        &mut self,
        assignment: &[f64],
    ) -> grb::Result<Option<Schedule>> {
        let vehicle_choices = self
            .problem
            .all_days()
            .map(|t| {
                self.problem
                    .all_customers()
                    .map(|i| {
                        let mut vehicle_delivery = self
                            .problem
                            .all_vehicles()
                            .map(|v| (v, assignment[self.vars.deliver_index(t, v, i)]))
                            .filter(|(_, delivery)| *delivery > Self::EPSILON)
                            .collect::<Vec<_>>();

                        vehicle_delivery.sort_unstable_by(|(_, delivery_a), (_, delivery_b)| {
                            delivery_b.partial_cmp(delivery_a).unwrap()
                        });

                        vehicle_delivery
                            .into_iter()
                            .map(|(vehicle, _)| vehicle)
                            .collect::<Vec<_>>()
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        eprintln!("# Compute new deliveries, time {}", self.elapsed_seconds());
        self.deliveries.set_all_statuses(|t, v, i| {
            let customer_vehicle_choices = &vehicle_choices[t as usize][i as usize - 1];
            if customer_vehicle_choices.is_empty() {
                false
            } else {
                customer_vehicle_choices[0] == v
            }
        })?;
        let opt_deliveries = match self.deliveries.solve()? {
            Some(deliveries) => Some(deliveries),
            None => match self.adjust_deliveries(assignment)? {
                Some(mut deliveries) => {
                    eprintln!("# Attempting fallback heuristic...");
                    let fixed = self.fix_deliveries_fallback(&mut deliveries)?;
                    if fixed {
                        Some(deliveries)
                    } else {
                        None
                    }
                }
                None => None,
            },
        };

        Ok(opt_deliveries.map(|deliveries| self.get_schedule_heuristically(&deliveries)))
    }

    /// A really stupid heuristic that does not take anything into consideration that would be sane
    fn fix_deliveries_fallback(&mut self, deliveries: &mut Deliveries) -> grb::Result<bool> {
        eprintln!("# Fix deliveries, time {}", self.elapsed_seconds());
        for t in self.problem.all_days() {
            eprintln!("# Fixing day {t}");

            // Compute load for each vehicle
            let load = self
                .problem
                .all_vehicles()
                .map(|v| {
                    self.problem
                        .all_customers()
                        .map(|i| deliveries.get(t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!("## Loads {t}: {load:?} capacity {}", self.problem.capacity);

            // Move delivieries to the same customer on different routes to the first route
            // with such a delivery. (Note that this does not change anything for customers
            // that are visited at most once.)
            for i in self.problem.all_customers() {
                let mut total_delivery: usize = self
                    .problem
                    .all_vehicles()
                    .map(|v| deliveries.get(t, v, i))
                    .sum();

                for v in self.problem.all_vehicles() {
                    let value = deliveries.get(t, v, i);
                    if value > 0 {
                        deliveries.set(t, v, i, total_delivery);
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
                        .map(|i| deliveries.get(t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!("## Loads {t}: {load:?} capacity {}", self.problem.capacity);

            // Sort vehicles by descending load
            let mut sorted_vehicles = self.problem.all_vehicles().collect::<Vec<_>>();
            sorted_vehicles.sort_by_cached_key(|a| Reverse(load[*a as usize]));

            // Move the smallest deliveries that overload the capacity to the next route
            for (i, v_ref) in sorted_vehicles[..sorted_vehicles.len() - 1]
                .iter()
                .enumerate()
            {
                let v = *v_ref;
                let v_next = sorted_vehicles[i + 1];
                if load[v as usize] > self.problem.capacity {
                    let mut heap = deliveries
                        .get_all_delivered_customers(t, v)
                        .into_iter()
                        .map(|i| {
                            Reverse(Delivery {
                                quantity: deliveries.get(t, v, i),
                                customer: i,
                            })
                        })
                        .collect::<BinaryHeap<_>>();

                    while load[v as usize] > self.problem.capacity {
                        if let Some(delivery) = heap.pop() {
                            let amount = delivery.0.quantity;
                            load[v as usize] -= amount;
                            load[v_next as usize] += amount;

                            deliveries.change_vehicle(t, v, v_next, delivery.0.customer);
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
                        .map(|i| deliveries.get(t, v, i))
                        .sum()
                })
                .collect::<Vec<usize>>();
            eprintln!("## Loads {t}: {load:?} capacity {}", self.problem.capacity);

            if load[sorted_vehicles[sorted_vehicles.len() - 1] as usize] > self.problem.capacity {
                return Ok(false);
            }
        }

        Ok(true)
        // Ok(None)
    }

    /// Runs LKH heuristic on visited sites to get a feasible route
    fn get_schedule_heuristically(&self, deliveries: &Deliveries) -> Schedule {
        eprintln!(
            "# Get schedule heuristically, time {}",
            self.elapsed_seconds()
        );

        Schedule::new_via_heuristic(self.problem, deliveries)
    }

    fn run_heuristic(&mut self) -> grb::Result<()> {
        self.heuristic.solve(&mut self.deliveries)?;

        Ok(())
    }
}

impl<'a> grb::callback::Callback for SolverData<'a> {
    fn callback(&mut self, w: grb::callback::Where) -> grb::callback::CbResult {
        match w {
            grb::callback::Where::MIPSol(ctx) => {
                self.ncalls += 1;
                let assignment = ctx.get_solution(&self.vars.variables)?;
                eprintln!("# Incumbent {}!", self.ncalls);
                eprintln!("#    current obj: {}", ctx.obj()?);
                eprintln!("#       best obj: {}", ctx.obj_best()?);
                eprintln!("#           time: {}", self.elapsed_seconds());

                if self.is_new_solution_just_set {
                    self.is_new_solution_just_set = false;
                } else {
                    {
                        if PRINT_VARIABLE_VALUES {
                            self.varnames
                                .iter()
                                .zip(assignment.iter())
                                .filter(|(_, &value)| value > Self::EPSILON)
                                .for_each(|(var, value)| eprintln!("#   - {var}: {value}"));
                        }

                        match self.adjust_deliveries(&assignment)? {
                            Some(deliveries) => {
                                let schedule = self.get_schedule_heuristically(&deliveries);
                                self.update_best_solution(schedule);
                            }
                            None => {
                                eprintln!(
                                    "# Failed to find a feasible adjusted solution. \
                                   This is weird, because we are in MIPSol."
                                );
                            }
                        }
                    }

                    let new_subtour_constraints = self
                        .integral_subtour_elimination(&assignment, |constr| ctx.add_lazy(constr))?;

                    if !new_subtour_constraints {
                        let schedule = self.get_schedule(&assignment);
                        self.update_best_solution(schedule);
                    }

                    eprintln!("# Callback finish time: {}", self.elapsed_seconds());
                }
            }
            grb::callback::Where::MIPNode(ctx) => {
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
                    eprintln!("#                 time: {}", self.elapsed_seconds());

                    if true {
                        let assignment = ctx.get_solution(&self.vars.variables)?;
                        if PRINT_VARIABLE_VALUES {
                            self.varnames
                                .iter()
                                .zip(assignment.iter())
                                .filter(|(_, &value)| value > Self::EPSILON)
                                .for_each(|(var, value)| eprintln!("#   - {var}: {value}"));
                        }
                        match self.fractional_delivery_heuristic(&assignment)? {
                            Some(schedule) => {
                                self.update_best_solution(schedule);
                            }
                            None => {
                                eprintln!("# Failed to find a feasible solution.");
                            }
                        }
                    }

                    self.give_new_best_solution_to_solver(ctx)?;

                    eprintln!("# MIPNode finish time: {}", self.elapsed_seconds());
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
        let mut env = grb::Env::new("")?;
        env.set(grb::param::Threads, 1)?;

        let mut lp = grb::Model::with_env("irp", &env)?;

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

        lp.set_objective(0, grb::ModelSense::Minimize)?;

        let mut data = SolverData::new(problem, &mut lp, &env, &cpu)?;

        // route degree constraints
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_sites() {
                    let mut lhs = grb::expr::LinExpr::new();

                    lhs.add_term(-2.0, data.vars.visit(t, v, i));
                    for j in problem.all_sites_except(i) {
                        lhs.add_term(1.0, data.vars.route(t, v, i, j));
                    }

                    lp.add_constr(&format!("Rd_{t}_{v}_{i}"), grb::c!(lhs == 0.0))?;
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

                    lp.add_constr(&format!("DC_{t}_{v}_{i}"), grb::c!(lhs >= 0.0))?;
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
                lp.add_constr(&format!("V1d_{t}_{i}"), grb::c!(lhs <= 1.0))?;
            }
        }

        // glue: disable delivery if we do not visit
        for t in problem.all_days() {
            for v in problem.all_vehicles() {
                for i in problem.all_customers() {
                    let mut lhs = grb::expr::LinExpr::new();
                    lhs.add_term(problem.capacity as f64, data.vars.visit(t, v, i));
                    lhs.add_term(-1.0, data.vars.deliver(t, v, i));

                    lp.add_constr(&format!("Gdv_{t}_{v}_{i}"), grb::c!(lhs >= 0.0))?;
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

            lp.add_constr(&format!("Ifd_{t}"), grb::c!(lhs == value))?;
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

                lp.add_constr(&format!("Ifc_{t}_{i}"), grb::c!(lhs == value))?;
            }
        }

        data.run_heuristic()?;

        lp.optimize_with_callback(&mut data)?;
        Self::print_raw_solution(&data, &lp)?;
        let assignment = lp.get_obj_attr_batch(grb::attr::X, data.vars.variables.clone())?;
        let schedule = data.get_schedule(&assignment);
        let solution = Solution::new(problem, schedule, data.elapsed_seconds(), data.cpu);
        eprintln!("# Final solution");
        eprintln!("{solution}");

        Ok(())
    }

    fn print_raw_solution(data: &SolverData, lp: &grb::Model) -> grb::Result<()> {
        let status = lp.status()?;
        eprintln!("# MIP solution status: {status:?}");

        let objective = lp.get_attr(grb::attr::ObjVal)?;
        eprintln!("# MIP solution value: {objective}");

        if PRINT_VARIABLE_VALUES {
            for var in data.vars.variables.iter() {
                let name = lp.get_obj_attr(grb::attr::VarName, var)?;
                let value = lp.get_obj_attr(grb::attr::X, var)?;
                if value > SolverData::EPSILON {
                    eprintln!("#   - {name}: {value}");
                }
            }
        }

        Ok(())
    }
}

pub fn solve(problem: Problem, cpu: String) {
    Solver::solve(&problem, cpu).unwrap();
}
