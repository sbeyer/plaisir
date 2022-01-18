use crate::problem::*;
use rand::distributions::{Distribution, Uniform};
use rand_xoshiro::rand_core::SeedableRng;

/// Assigns (day, customer) to vehicle
#[derive(Debug)]
struct VehiclePlan(Vec<Vec<Option<usize>>>);

pub struct RandomHeuristic<'a> {
    problem: &'a Problem,
    rng: rand_xoshiro::Xoshiro128StarStar,
}

impl<'a> RandomHeuristic<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        const SEED: [u8; 16] = [
            42, 228, 59, 86, 175, 57, 79, 176, 13, 49, 245, 187, 66, 136, 74, 182,
        ];
        Self {
            problem,
            rng: rand_xoshiro::Xoshiro128StarStar::from_seed(SEED),
        }
    }

    pub fn solve(&mut self) {
        loop {
            let vehicle_plan = self.make_random_vehicle_plan();
            eprintln!("# plan {:?}", vehicle_plan);
        }
    }

    fn make_random_vehicle_plan(&mut self) -> VehiclePlan {
        let dist = Uniform::from(0..self.problem.num_vehicles);
        VehiclePlan(
            self.problem
                .all_days()
                .map(|_| {
                    self.problem
                        .all_customers()
                        .map(|_| Some(dist.sample(&mut self.rng)))
                        .collect()
                })
                .collect(),
        )
    }
}
