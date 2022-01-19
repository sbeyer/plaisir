use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use rand::distributions::{Distribution, Uniform};
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Debug)]
struct VehiclePlan(Vec<Vec<Option<VehicleId>>>);

pub struct RandomHeuristic<'a> {
    problem: &'a Problem,
    dist: Uniform<VehicleId>,
    rng: rand_xoshiro::Xoshiro128StarStar,
}

impl<'a> RandomHeuristic<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        const SEED: [u8; 16] = [
            42, 228, 59, 86, 175, 57, 79, 176, 13, 49, 245, 187, 66, 136, 74, 182,
        ];
        let rng = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED);
        let dist = Uniform::from(0..(problem.num_vehicles as VehicleId));
        Self { problem, dist, rng }
    }

    pub fn solve(&mut self, delivery_solver: &mut DeliverySolver) -> grb::Result<()> {
        loop {
            let vehicle_plan = self.make_random_vehicle_plan();
            eprintln!("# plan {:?}", vehicle_plan);

            delivery_solver.set_all_statuses(|t, v, i| {
                if let Some(vehicle_choice) = vehicle_plan.0[t as usize][i as usize - 1] {
                    vehicle_choice == v
                } else {
                    false
                }
            })?;
            let opt_deliveries = delivery_solver.solve()?;
            if let Some(deliveries) = opt_deliveries {
                eprintln!("# -> deliveries {:?}", deliveries);
            } else {
                eprintln!("# -> infeasible deliveries");
            }
        }
    }

    fn make_random_vehicle_plan(&mut self) -> VehiclePlan {
        VehiclePlan(
            self.problem
                .all_days()
                .map(|_| {
                    self.problem
                        .all_customers()
                        .map(|_| Some(self.dist.sample(&mut self.rng)))
                        .collect()
                })
                .collect(),
        )
    }
}
