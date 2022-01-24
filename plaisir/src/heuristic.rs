use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use crate::route::Solver as RouteSolver;
use crate::solution::*;
use rand::distributions::Uniform;
use rand::Rng;
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Clone, Debug)]
struct VehiclePlan(Vec<Vec<Option<VehicleId>>>);

/// Genetic algorithm for IRP
pub struct GeneticHeuristic<'a> {
    problem: &'a Problem,
    dist_day: Uniform<DayId>,
    dist_customer: Uniform<SiteId>,
    rng: rand_xoshiro::Xoshiro128StarStar,
}

impl<'a> GeneticHeuristic<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        const SEED: [u8; 16] = [
            42, 228, 59, 86, 175, 57, 79, 176, 13, 49, 245, 187, 66, 136, 74, 182,
        ];
        let rng = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED);

        let dist_day = Uniform::from(0..problem.num_days as DayId);
        let dist_customer = Uniform::from(0..problem.num_customers as SiteId);
        Self {
            problem,
            dist_day,
            dist_customer,
            rng,
        }
    }

    pub fn solve(
        &mut self,
        delivery_solver: &mut DeliverySolver,
        route_solver: &mut RouteSolver,
        solution_pool: &mut SolutionPool,
    ) -> grb::Result<()> {
        // This threshold should be understood as follows: consider the best and worst solution in
        // the current solution pool, and a new best solution. Scale these values such that the new
        // best solution is zero and the worst solution is one. The (previous) best solution lies
        // in between. If that solution is below this threshold, we have a neglectable improvement.
        const NEGLECTABLE_IMPROVEMENT_THRESHOLD: f64 = 0.1;

        /// Stop the algorithm after this number of iterations that led to either no or only small
        /// improvements... Note that the actual number of iterations can be smaller due to other
        /// factors, e.g. the number of infeasible solutions found.
        const MAX_NEGLECTABLE_IMPROVEMENT_ITERATIONS: usize = 3000;

        /// Guarantee at least this number of iterations that led to either no or only small
        /// improvements.
        const MIN_NEGLECTABLE_IMPROVEMENT_ITERATIONS: usize = 100;

        if solution_pool.solutions.len() < 5 {
            return Ok(());
        }

        let mut previous_best_value = solution_pool.get_best().unwrap().value();

        eprintln!("# GeneticHeuristic start");
        let mut count_iteration = 0;
        let mut count_infeasible = 0;
        let mut count_neglectable_improvement = 0;
        loop {
            let sol_idx1 = self.rng.gen_range(0..solution_pool.solutions.len());
            let mut sol_idx2 = self.rng.gen_range(0..solution_pool.solutions.len());
            while sol_idx2 == sol_idx1 {
                sol_idx2 = self.rng.gen_range(0..solution_pool.solutions.len());
            }

            let day_idx1 = self.rng.sample(self.dist_day);

            let customer_idx_range = (
                self.rng.sample(self.dist_customer),
                self.rng.sample(self.dist_customer),
            );
            let customer_idx_range = if customer_idx_range.0 < customer_idx_range.1 {
                (customer_idx_range.0, customer_idx_range.1)
            } else {
                (customer_idx_range.1, customer_idx_range.0)
            };

            let solution1 = &solution_pool.solutions[sol_idx1];
            let solution2 = &solution_pool.solutions[sol_idx2];

            let vehicle_plan1 = self.create_vehicle_plan_from_solution(solution1);

            let vehicle_plans = self
                .problem
                .all_days()
                .map(|source_day| {
                    let mut vehicle_plan = vehicle_plan1.clone();
                    self.crossover_vehicle_plan(
                        &mut vehicle_plan,
                        day_idx1 as DayId,
                        solution2,
                        source_day,
                        customer_idx_range,
                    );
                    vehicle_plan
                })
                .collect::<Vec<VehiclePlan>>();

            for vehicle_plan in vehicle_plans.iter() {
                count_iteration += 1;

                // Compute deliveries from vehicle plan
                delivery_solver.set_all_statuses(|t, v, i| {
                    if let Some(vehicle_choice) = vehicle_plan.0[t as usize][i as usize - 1] {
                        vehicle_choice == v
                    } else {
                        false
                    }
                })?;
                let opt_deliveries = delivery_solver.solve()?;

                if let Some(deliveries) = opt_deliveries {
                    let schedule =
                        Schedule::new_via_heuristic(self.problem, &deliveries, route_solver);
                    let worst_value = solution_pool.get_worst_value();
                    let (new_best, opt_solution) = solution_pool.add(self.problem, schedule);
                    if new_best {
                        let solution = opt_solution.unwrap();
                        eprintln!(
                            "# New best solution of value {} (old: {})",
                            solution.value(),
                            previous_best_value
                        );
                        let previous_scaled_for_threshold = (previous_best_value
                            - solution.value())
                            / (worst_value - solution.value());
                        /*
                        eprintln!(
                            "# NEW BEST = {}, PREVIOUS BEST = {}, WORST = {}, SCALED = {}, NEGLECTABLE = {:?}",
                            solution.value(),
                            previous_best_value,
                            worst_value,
                            previous_scaled_for_threshold,
                            previous_scaled_for_threshold < NEGLECTABLE_IMPROVEMENT_THRESHOLD
                        );
                        */
                        if previous_scaled_for_threshold < NEGLECTABLE_IMPROVEMENT_THRESHOLD {
                            count_neglectable_improvement += 1;
                        } else {
                            count_neglectable_improvement = 0;
                            count_infeasible = 0;
                        }
                        previous_best_value = solution.value();
                    } else {
                        count_neglectable_improvement += 1;
                    }
                } else {
                    count_infeasible += 1;
                    count_neglectable_improvement += 1;
                }

                if count_iteration % 100 == 50 {
                    eprintln!("# GeneticHeuristic Iteration {count_iteration} ({count_infeasible} infeasible of {count_neglectable_improvement} neglectable improvements)");
                }
            }

            if count_neglectable_improvement >= MIN_NEGLECTABLE_IMPROVEMENT_ITERATIONS
                && count_neglectable_improvement as f64
                    >= MAX_NEGLECTABLE_IMPROVEMENT_ITERATIONS as f64
                        * (1.0 - (count_infeasible as f64 / count_neglectable_improvement as f64))
            {
                break;
            }
        }

        eprintln!("# GeneticHeuristic end after {count_iteration} iterations ({count_infeasible} infeasible of {count_neglectable_improvement} neglectable improvements)");
        Ok(())
    }

    fn create_vehicle_plan_from_solution(&self, solution: &Solution) -> VehiclePlan {
        VehiclePlan(
            self.problem
                .all_days()
                .map(|t| {
                    let mut day_plan = vec![None; self.problem.num_customers];

                    for v in self.problem.all_vehicles() {
                        let route = solution.route(t, v);
                        for delivery in route.iter() {
                            debug_assert_ne!(delivery.customer, 0);

                            day_plan[delivery.customer as usize - 1] = Some(v);
                        }
                    }

                    day_plan
                })
                .collect(),
        )
    }

    fn crossover_vehicle_plan(
        &self,
        vehicle_plan: &mut VehiclePlan,
        target_idx: DayId,
        crossover_solution: &Solution,
        crossover_day: DayId,
        customer_idx_range: (SiteId, SiteId),
    ) {
        // Prepare crossover by cleaning the range
        #[allow(clippy::needless_range_loop)]
        for i in customer_idx_range.0..=customer_idx_range.1 {
            vehicle_plan.0[target_idx as usize][i as usize] = None;
        }

        // Apply crossover
        for v in self.problem.all_vehicles() {
            let route = crossover_solution.route(crossover_day, v);
            for delivery in route.iter() {
                debug_assert_ne!(delivery.customer, 0);

                let i = delivery.customer;
                if i >= customer_idx_range.0 && i <= customer_idx_range.1 {
                    vehicle_plan.0[target_idx as usize][i as usize] = Some(v);
                }
            }
        }
    }
}
