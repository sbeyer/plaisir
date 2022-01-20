use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use crate::solution::*;
use rand::distributions::{Distribution, Uniform};
use rand::Rng;
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Debug)]
struct VehiclePlan(Vec<Vec<Option<VehicleId>>>);

pub struct RandomHeuristic<'a> {
    problem: &'a Problem,
    dist_vehicle: Uniform<VehicleId>,
    dist_visit: Uniform<f64>,
    rng_vehicle: rand_xoshiro::Xoshiro128StarStar,
    rng_visit_threshold: rand_xoshiro::Xoshiro128StarStar,
    rng_visit: rand_xoshiro::Xoshiro128StarStar,
}

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

    pub fn solve(&mut self, delivery_solver: &mut DeliverySolver) -> grb::Result<()> {
        let mut best_solution = Solution::empty();

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
                    &mut best_solution,
                    5000,
                    1000,
                )?;

                let threshold_old = increasing_threshold_bounds[i];
                increasing_threshold_bounds[i] *= 0.8 * (1.0 - infeasible_freq) + infeasible_freq;

                // Finish up if threshold didn't change too much
                if approx::abs_diff_eq!(
                    threshold_old,
                    increasing_threshold_bounds[i],
                    epsilon = 0.01
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
                .map(|lb| self.rng_visit_threshold.gen_range(*lb..1.0))
                .collect::<Vec<_>>();

            self.threshold_loop(
                &thresholds,
                delivery_solver,
                &mut best_solution,
                10000,
                2000,
            )?;
        }
    }

    /// Make a loop for a specific threshold configurations
    fn threshold_loop(
        &mut self,
        thresholds: &[f64],
        delivery_solver: &mut DeliverySolver,
        best_solution: &mut Solution,
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
                let solution = Solution::new(
                    self.problem,
                    schedule,
                    1337.0,                           // TODO
                    "Intel Core i5-10210U @ 1.60GHz", // TODO
                );

                if solution.value() < best_solution.value() {
                    eprintln!(
                        "# New best solution of value {} (improving {}) with visit probabilities {thresholds:?} after {counter_no_improvement} iterations since last best value",
                        solution.value(),
                        best_solution.value(),
                    );
                    eprintln!("{}", solution);
                    *best_solution = solution;
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
