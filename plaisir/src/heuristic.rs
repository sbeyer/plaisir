use crate::delivery::Solver as DeliverySolver;
use crate::problem::*;
use crate::route::Solver as RouteSolver;
use crate::solution::*;
use rand::distributions::Uniform;
use rand::Rng;
use rand_xoshiro::rand_core::SeedableRng;

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

/// The probability that the best solution is chosen for crossover.
const ALPHA_MALE_THRESHOLD: f64 = 0.666667;

/// Assigns (day, customer) to vehicle
#[derive(Clone, Debug)]
struct VehiclePlan(Vec<VehicleDayPlan>);

#[derive(Clone, Debug)]
struct VehicleDayPlan(Vec<Option<VehicleId>>);

impl VehicleDayPlan {
    fn new(problem: &Problem) -> Self {
        Self(vec![None; problem.num_customers])
    }

    fn set(&mut self, customer: SiteId, vehicle: VehicleId) {
        self.0[customer as usize - 1] = Some(vehicle);
    }

    fn get(&self, i: SiteId) -> Option<VehicleId> {
        self.0[i as usize - 1]
    }

    fn clear(&mut self, (start, end): (SiteId, SiteId)) {
        for site in start..=end {
            let idx = site as usize - 1;
            self.0[idx] = None;
        }
    }
}

impl VehiclePlan {
    fn from_solution(problem: &Problem, solution: &Solution) -> Self {
        Self(
            problem
                .all_days()
                .map(|t| {
                    let mut day_plan = VehicleDayPlan::new(problem);

                    for v in problem.all_vehicles() {
                        let route = solution.route(t, v);
                        for delivery in route.iter() {
                            debug_assert_ne!(delivery.customer, 0);

                            day_plan.set(delivery.customer, v);
                        }
                    }

                    day_plan
                })
                .collect(),
        )
    }

    fn get_vehicle(&self, t: DayId, i: SiteId) -> Option<VehicleId> {
        self.0[t as usize].get(i)
    }

    fn set(&mut self, t: DayId, i: SiteId, vehicle: VehicleId) {
        self.0[t as usize].set(i, vehicle);
    }

    fn clear(&mut self, t: DayId, customer_range: (SiteId, SiteId)) {
        self.0[t as usize].clear(customer_range);
    }
}

/// Genetic algorithm for IRP
pub struct GeneticHeuristic<'a> {
    problem: &'a Problem,
    dist_day: Uniform<DayId>,
    dist_vehicle: Uniform<VehicleId>,
    dist_customer: Uniform<SiteId>,
    dist_01: Uniform<f64>,
    rng: rand_xoshiro::Xoshiro128StarStar,
}

impl<'a> GeneticHeuristic<'a> {
    pub fn new(problem: &'a Problem) -> Self {
        const SEED: [u8; 16] = [
            42, 228, 59, 86, 175, 57, 79, 176, 13, 49, 245, 187, 66, 136, 74, 182,
        ];
        let rng = rand_xoshiro::Xoshiro128StarStar::from_seed(SEED);

        let dist_day = Uniform::from(0..problem.num_days as DayId);
        let dist_vehicle = Uniform::from(0..problem.num_vehicles as VehicleId);
        let dist_customer = Uniform::from(1..problem.num_sites as SiteId);
        let dist_01 = Uniform::<f64>::from(0.0..1.0);
        Self {
            problem,
            dist_day,
            dist_vehicle,
            dist_customer,
            dist_01,
            rng,
        }
    }

    pub fn solve(
        &mut self,
        delivery_solver: &mut DeliverySolver,
        route_solver: &mut RouteSolver,
        solution_pool: &mut SolutionPool,
    ) -> grb::Result<()> {
        if solution_pool.solutions.len() < 3 {
            return Ok(());
        }

        let mut previous_best_value = solution_pool.get_best().unwrap().value();

        eprintln!("# GeneticHeuristic start");
        let mut count_iteration = 0;
        let mut count_infeasible = 0;
        let mut count_neglectable_improvement = 0;
        loop {
            let alpha_male_crossover = self.rng.sample(self.dist_01) < ALPHA_MALE_THRESHOLD;
            let (solution1, solution2) = if alpha_male_crossover {
                let sol_idx = self.rng.gen_range(0..solution_pool.solutions.len());

                (
                    solution_pool.get_best().unwrap(),
                    &solution_pool.solutions[sol_idx],
                )
            } else {
                let sol_idx1 = self.rng.gen_range(0..solution_pool.solutions.len());
                let mut sol_idx2 = self.rng.gen_range(0..solution_pool.solutions.len());
                while sol_idx2 == sol_idx1 {
                    sol_idx2 = self.rng.gen_range(0..solution_pool.solutions.len());
                }

                (
                    &solution_pool.solutions[sol_idx1],
                    &solution_pool.solutions[sol_idx2],
                )
            };

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

            let vehicle_plan1 = VehiclePlan::from_solution(self.problem, solution1);

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
                self.try_solution_from_vehicle_plan(
                    vehicle_plan,
                    delivery_solver,
                    route_solver,
                    solution_pool,
                    &mut previous_best_value,
                    &mut count_iteration,
                    &mut count_neglectable_improvement,
                    &mut count_infeasible,
                )?;
            }

            if alpha_male_crossover
                && count_iteration == count_neglectable_improvement
                && count_iteration >= MIN_NEGLECTABLE_IMPROVEMENT_ITERATIONS
            {
                eprintln!("# GeneticHeuristic Mutation Intermezzo");
                for i in customer_idx_range.0..=customer_idx_range.1 {
                    let mut vehicle_plan = vehicle_plan1.clone();
                    let opt_vehicle = vehicle_plan.get_vehicle(day_idx1, i);
                    if opt_vehicle.is_some() {
                        vehicle_plan.clear(day_idx1, (i, i));
                    } else {
                        let new_vehicle = self.rng.sample(self.dist_vehicle);
                        vehicle_plan.set(day_idx1, i, new_vehicle);
                    }

                    self.try_solution_from_vehicle_plan(
                        &vehicle_plan,
                        delivery_solver,
                        route_solver,
                        solution_pool,
                        &mut previous_best_value,
                        &mut count_iteration,
                        &mut count_neglectable_improvement,
                        &mut count_infeasible,
                    )?;
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

    fn crossover_vehicle_plan(
        &self,
        vehicle_plan: &mut VehiclePlan,
        target_idx: DayId,
        crossover_solution: &Solution,
        crossover_day: DayId,
        customer_idx_range: (SiteId, SiteId),
    ) {
        // Prepare crossover by cleaning the range
        vehicle_plan.clear(target_idx, customer_idx_range);

        // Apply crossover
        for v in self.problem.all_vehicles() {
            let route = crossover_solution.route(crossover_day, v);
            for delivery in route.iter() {
                debug_assert_ne!(delivery.customer, 0);

                let i = delivery.customer;
                if i >= customer_idx_range.0 && i <= customer_idx_range.1 {
                    vehicle_plan.set(target_idx, i, v);
                }
            }
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn try_solution_from_vehicle_plan(
        &mut self,
        vehicle_plan: &VehiclePlan,
        delivery_solver: &mut DeliverySolver,
        route_solver: &mut RouteSolver,
        solution_pool: &mut SolutionPool,
        previous_best_value: &mut f64,
        count_iteration: &mut usize,
        count_neglectable_improvement: &mut usize,
        count_infeasible: &mut usize,
    ) -> grb::Result<()> {
        *count_iteration += 1;

        // Compute deliveries from vehicle plan
        delivery_solver.set_all_statuses(|t, v, i| {
            if let Some(vehicle_choice) = vehicle_plan.get_vehicle(t, i) {
                vehicle_choice == v
            } else {
                false
            }
        })?;
        let opt_deliveries = delivery_solver.solve()?;

        if let Some(deliveries) = opt_deliveries {
            let schedule = Schedule::new_via_heuristic(self.problem, &deliveries, route_solver);
            let worst_value = solution_pool.get_worst_value();
            let (new_best, opt_solution) = solution_pool.add(self.problem, schedule);
            if new_best {
                let solution = opt_solution.unwrap();
                eprintln!(
                    "# New best solution of value {} (old: {})",
                    solution.value(),
                    *previous_best_value
                );
                let previous_scaled_for_threshold =
                    (*previous_best_value - solution.value()) / (worst_value - solution.value());
                if previous_scaled_for_threshold < NEGLECTABLE_IMPROVEMENT_THRESHOLD {
                    *count_neglectable_improvement += 1;
                } else {
                    *count_neglectable_improvement = 0;
                    *count_infeasible = 0;
                }
                *previous_best_value = solution.value();
            } else {
                *count_neglectable_improvement += 1;
            }
        } else {
            *count_infeasible += 1;
            *count_neglectable_improvement += 1;
        }

        if *count_iteration % 100 == 50 {
            eprintln!("# GeneticHeuristic Iteration {count_iteration} ({count_infeasible} infeasible of {count_neglectable_improvement} neglectable improvements)");
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod vehicle_plan {
        use super::*;

        #[test]
        fn set_and_get() {
            let mut vp = VehicleDayPlan(vec![None, None, None, None]);

            vp.set(3, 2);
            assert_eq!(vp.get(3), Some(2));

            vp.set(4, 2);
            assert_eq!(vp.get(4), Some(2));

            vp.set(2, 1);
            assert_eq!(vp.get(2), Some(1));

            vp.set(1, 0);
            assert_eq!(vp.get(1), Some(0));

            vp.clear((3, 3));
            assert_eq!(vp.get(1), Some(0));
            assert_eq!(vp.get(2), Some(1));
            assert_eq!(vp.get(3), None);
            assert_eq!(vp.get(4), Some(2));

            vp.clear((2, 4));
            assert_eq!(vp.get(1), Some(0));
            assert_eq!(vp.get(2), None);
            assert_eq!(vp.get(3), None);
            assert_eq!(vp.get(4), None);

            vp.set(2, 0);
            vp.set(3, 0);
            vp.set(4, 0);

            vp.clear((2, 3));
            assert_eq!(vp.get(1), Some(0));
            assert_eq!(vp.get(2), None);
            assert_eq!(vp.get(3), None);
            assert_eq!(vp.get(4), Some(0));
        }
    }
}
