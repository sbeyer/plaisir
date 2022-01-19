use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use crate::solution::*;
use rand_distr::{Distribution, SkewNormal, Uniform};
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Debug)]
struct VehiclePlan(Vec<Vec<Option<VehicleId>>>);

pub struct RandomHeuristic<'a> {
    problem: &'a Problem,
    dist_vehicle: Uniform<VehicleId>,
    dist_visit_threshold: SkewNormal<f64>,
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
        let dist_visit_threshold = SkewNormal::<f64>::new(0.95, 0.3, -5.0).unwrap();
        let dist_visit = Uniform::from(0f64..1f64);
        Self {
            problem,
            dist_vehicle,
            dist_visit_threshold,
            dist_visit,
            rng_vehicle,
            rng_visit_threshold,
            rng_visit,
        }
    }

    pub fn solve(&mut self, delivery_solver: &mut DeliverySolver) -> grb::Result<()> {
        let start_time = std::time::Instant::now();
        let mut best_solution = Solution::empty();
        let mut counter_no_improvement = 0usize;
        let mut counter_infeasible = 0usize;
        loop {
            let threshold = self
                .dist_visit_threshold
                .sample(&mut self.rng_visit_threshold);
            let threshold = if threshold < 0.1 {
                1.0 - threshold
            } else {
                threshold
            };
            let vehicle_plan = self.make_random_vehicle_plan(threshold);
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
                    start_time.elapsed().as_millis() as f64 * 1e-3,
                    "Intel Core i5-10210U @ 1.60GHz", // TODO
                );

                if solution.value() < best_solution.value() {
                    eprintln!(
                        "# New best solution of value {} (improving {}) with visit probability {} after {} iterations since last best value",
                        solution.value(),
                        best_solution.value(),
                        threshold * 100.0,
                        counter_no_improvement
                    );
                    eprintln!("{}", solution);
                    best_solution = solution;
                    counter_no_improvement = 0;
                    counter_infeasible = 0;
                }
            } else {
                counter_infeasible += 1;
            }

            if counter_no_improvement > 0 && counter_no_improvement % 1000 == 0 {
                eprintln!("# No new best solution found after {counter_no_improvement} iterations of which {counter_infeasible} were infeasible");
            }

            counter_no_improvement += 1;
        }
    }

    fn make_random_vehicle_plan(&mut self, threshold: f64) -> VehiclePlan {
        VehiclePlan(
            self.problem
                .all_days()
                .map(|_| {
                    self.problem
                        .all_customers()
                        .map(|_| {
                            let value = self.dist_visit.sample(&mut self.rng_visit);
                            if value < threshold {
                                Some(self.dist_vehicle.sample(&mut self.rng_vehicle))
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
