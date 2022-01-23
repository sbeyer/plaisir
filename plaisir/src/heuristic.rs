use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use crate::solution::*;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Clone, Debug)]
struct VehiclePlan(Vec<Vec<Option<VehicleId>>>);

#[allow(dead_code)]
pub struct RandomHeuristic<'a> {
    problem: &'a Problem,
    dist_vehicle: Uniform<VehicleId>,
    dist_visit: Uniform<f64>,
    rng_vehicle: rand_xoshiro::Xoshiro128StarStar,
    rng_visit_threshold: rand_xoshiro::Xoshiro128StarStar,
    rng_visit: rand_xoshiro::Xoshiro128StarStar,
}

#[allow(dead_code)]
impl<'a> RandomHeuristic<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        const SEED_VEHICLE: [u8; 16] = [
            42, 228, 59, 86, 175, 57, 79, 176, 13, 49, 245, 187, 66, 136, 74, 182,
        ];
        const SEED_VISIT_THRESHOLD: [u8; 16] = [
            220, 169, 125, 15, 9, 75, 254, 75, 143, 241, 88, 78, 76, 61, 234, 233,
        ];
        const SEED_VISIT: [u8; 16] = [
            107, 106, 176, 127, 47, 108, 48, 99, 115, 98, 130, 167, 146, 167, 44, 165,
        ];
        let rng_vehicle = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED_VEHICLE);
        let rng_visit_threshold = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED_VISIT_THRESHOLD);
        let rng_visit = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED_VISIT);
        let dist_vehicle = Uniform::from(0..(problem.num_vehicles as VehicleId));
        let dist_visit = Uniform::from(0f64..1f64);
        Self {
            problem,
            dist_vehicle,
            dist_visit,
            rng_vehicle,
            rng_visit_threshold,
            rng_visit,
        }
    }

    pub fn solve(
        &mut self,
        delivery_solver: &mut DeliverySolver,
        solution_pool: &mut SolutionPool,
    ) -> grb::Result<()> {
        const LOWER_BOUND_SEARCH_DESCENT: f64 = 0.67;
        const LOWER_BOUND_SEARCH_ABSDIFF: f64 = 0.01;

        eprintln!(
            "# RandomHeuristic Step 1: Find good lower bounds for thresholds (increasing by day)"
        );
        let mut increasing_threshold_bounds = vec![1.0; self.problem.num_days];
        for i in 0..increasing_threshold_bounds.len() {
            eprintln!("# RandomHeuristic Step 1.{i}: Find good lower bound for threshold {i}");
            loop {
                let infeasible_freq = self.threshold_loop(
                    &increasing_threshold_bounds,
                    delivery_solver,
                    solution_pool,
                    5000,
                    1000,
                )?;

                let threshold_old = increasing_threshold_bounds[i];
                increasing_threshold_bounds[i] *=
                    LOWER_BOUND_SEARCH_DESCENT * (1.0 - infeasible_freq) + infeasible_freq;

                // Finish up if threshold didn't change too much
                if approx::abs_diff_eq!(
                    threshold_old,
                    increasing_threshold_bounds[i],
                    epsilon = LOWER_BOUND_SEARCH_ABSDIFF
                ) {
                    break;
                }
            }
        }
        eprintln!("# RandomHeuristic Finished Step 1 with {increasing_threshold_bounds:?}");

        eprintln!("#");
        eprintln!("# RandomHeuristic Step 2: Sample thresholds and try to find better solutions");
        loop {
            let thresholds = increasing_threshold_bounds
                .iter()
                .map(|lb| self.rng_visit_threshold.gen_range(*lb..=1.0))
                .collect::<Vec<_>>();

            self.threshold_loop(&thresholds, delivery_solver, solution_pool, 10000, 2000)?;
        }
    }

    /// Make a loop for a specific threshold configurations
    fn threshold_loop(
        &mut self,
        thresholds: &[f64],
        delivery_solver: &mut DeliverySolver,
        solution_pool: &mut SolutionPool,
        max_count_no_improvement: usize,
        max_count_infeasibles: usize,
    ) -> grb::Result<f64> {
        eprintln!("# Finding random solutions with visit thresholds {thresholds:?}",);
        let mut counter_no_improvement = 0usize;
        let mut counter_infeasible = 0usize;
        loop {
            let vehicle_plan = self.make_random_vehicle_plan(thresholds);
            //eprintln!("# plan {:?}", vehicle_plan);

            delivery_solver.set_all_statuses(|t, v, i| {
                if let Some(vehicle_choice) = vehicle_plan.0[t as usize][i as usize - 1] {
                    vehicle_choice == v
                } else {
                    false
                }
            })?;
            let opt_deliveries = delivery_solver.solve()?;
            if let Some(deliveries) = opt_deliveries {
                //eprintln!("# -> deliveries {deliveries:?}");
                let schedule = Schedule::new_via_heuristic(self.problem, &deliveries);
                //eprintln!("# -> schedule {schedule:?}");
                // TODO: struct SolutionPool or something like that
                let (new_best, opt_solution) = solution_pool.add(self.problem, schedule);

                if new_best {
                    let solution = opt_solution.unwrap();
                    eprintln!(
                        "# New best solution of value {} with visit probabilities {thresholds:?} after {counter_no_improvement} iterations since last best value",
                        solution.value(),
                    );
                    counter_no_improvement = 0;
                    counter_infeasible = 0;
                } else {
                    counter_no_improvement += 1;
                }
            } else {
                counter_infeasible += 1;
                counter_no_improvement += 1;
            }

            if counter_no_improvement > 0 && counter_no_improvement % 1000 == 0 {
                eprintln!("# No new best solution found after {counter_no_improvement} iterations of which {counter_infeasible} were infeasible");
            }
            if counter_no_improvement >= max_count_no_improvement {
                eprintln!("# Stop attempts using these thresholds due to too many iterations ({counter_no_improvement}) without improvements, having {counter_infeasible} infeasibles");

                break;
            }
            if counter_infeasible >= max_count_infeasibles {
                eprintln!("# Stop attempts using these thresholds due to too many infeasibles ({counter_infeasible}) during {counter_no_improvement} iterations without improvement");
                break;
            }
        }

        Ok(counter_infeasible as f64 / counter_no_improvement as f64)
    }

    /// Make a random vehicle plan based on visit probabilities `thresholds` (one for each day)
    fn make_random_vehicle_plan(&mut self, thresholds: &[f64]) -> VehiclePlan {
        debug_assert_eq!(self.problem.all_days().count(), thresholds.len());

        // Notes:
        //  - The first used vehicles per day are always 0,1,...,#vehicles
        //    in ascending order.
        VehiclePlan(
            thresholds
                .iter()
                .map(|threshold| {
                    let mut vehicles_used = 0;
                    self.problem
                        .all_customers()
                        .map(|_| {
                            let value = self.dist_visit.sample(&mut self.rng_visit);
                            if value < *threshold {
                                let mut vehicle = self.dist_vehicle.sample(&mut self.rng_vehicle);
                                if vehicles_used < vehicle {
                                    vehicle = vehicles_used;
                                    vehicles_used += 1;
                                }
                                Some(vehicle)
                            } else {
                                None
                            }
                        })
                        .collect()
                })
                .collect(),
        )
    }
}

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
        solution_pool: &mut SolutionPool,
    ) -> grb::Result<()> {
        if solution_pool.solutions.len() < 5 {
            return Ok(());
        }

        eprintln!("# GeneticHeuristic start");
        let mut count_iteration = 0;
        let mut count_infeasible = 0;
        let mut count_no_improvement = 0;
        loop {
            count_iteration += 1;

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
                    let schedule = Schedule::new_via_heuristic(self.problem, &deliveries);
                    let (new_best, opt_solution) = solution_pool.add(self.problem, schedule);
                    if new_best {
                        let solution = opt_solution.unwrap();
                        eprintln!("# New best solution of value {}", solution.value(),);
                        count_no_improvement = 0;
                        count_infeasible = 0;
                    } else {
                        count_no_improvement += 1;
                    }
                } else {
                    count_infeasible += 1;
                    count_no_improvement += 1;
                }

                if count_iteration % 100 == 50 {
                    eprintln!("# GeneticHeuristic Iteration {count_iteration} (#{count_infeasible} infeasible of #{count_no_improvement} no improvement)");
                }
            }

            if count_infeasible >= 1000 || count_no_improvement >= 2000 {
                break;
            }
        }

        eprintln!("# GeneticHeuristic end after {count_iteration} iterations (#{count_infeasible} infeasible of #{count_no_improvement} no improvement)");
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
                        for delivery in route.iter().skip(1) {
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
            for delivery in route.iter().skip(1) {
                debug_assert_ne!(delivery.customer, 0);

                let i = delivery.customer;
                if i >= customer_idx_range.0 && i <= customer_idx_range.1 {
                    vehicle_plan.0[target_idx as usize][i as usize] = Some(v);
                }
            }
        }
    }
}
