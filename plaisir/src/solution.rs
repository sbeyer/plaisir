use crate::delivery::{Deliveries, Delivery};
use crate::problem::{DayId, Problem, SiteId, VehicleId};
use crate::route::Solver as RouteSolver;
use std::fmt;
use std::time;

pub type Route = Vec<Delivery>;

#[derive(Debug)]
pub struct Schedule(pub Vec<Vec<Route>>);

impl Schedule {
    pub fn new_via_heuristic(problem: &Problem, deliveries: &Deliveries) -> Self {
        Schedule(
            problem
                .all_days()
                .map(|t| {
                    problem
                        .all_vehicles()
                        .map(|v| {
                            let tsp_tour = RouteSolver::solve(problem, deliveries, t, v);

                            tsp_tour
                                .iter()
                                .map(|site| Delivery {
                                    quantity: if *site == 0 {
                                        0
                                    } else {
                                        deliveries.get(t, v, *site as SiteId)
                                    },
                                    customer: *site as SiteId,
                                })
                                .collect()
                        })
                        .collect()
                })
                .collect(),
        )
    }
}

pub struct Solution {
    schedule: Schedule,
    cost_transportation: f64,
    cost_inventory_depot: f64,
    cost_inventory_customers: f64,
    cost_total: f64,
    processor: String,
    time: f64,
}

impl Solution {
    pub fn empty() -> Self {
        Self {
            schedule: Schedule(Vec::new()),
            cost_transportation: f64::INFINITY,
            cost_inventory_depot: f64::INFINITY,
            cost_inventory_customers: f64::INFINITY,
            cost_total: f64::INFINITY,
            processor: String::new(),
            time: 0.0,
        }
    }

    pub fn new(problem: &Problem, schedule: Schedule, time: f64, cpu: &str) -> Self {
        let mut sol = Solution {
            schedule,
            cost_transportation: 0.,
            cost_inventory_depot: 0.,
            cost_inventory_customers: 0.,
            cost_total: 0.,
            processor: cpu.to_string(),
            time,
        };

        // transportation cost
        for day_schedule in sol.schedule.0.iter() {
            for route in day_schedule.iter() {
                let mut tour: Vec<SiteId> = route.iter().map(|x| x.customer).collect();
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

        for day_schedule in sol.schedule.0.iter() {
            // step one: deliveries
            for route in day_schedule.iter() {
                for Delivery { quantity, customer } in route.iter() {
                    inventory[*customer as usize] += *quantity as f64;
                    inventory[0] -= *quantity as f64;
                }
            }

            // step two: daily change (production at depot, consumption at customers)
            for i in problem.all_sites() {
                inventory[i as usize] += problem.site(i).level_change();
            }

            // update inventory costs
            for i in problem.all_sites() {
                cost_inventory[i as usize] += problem.site(i).cost() * inventory[i as usize]
            }
        }

        sol.cost_inventory_depot = cost_inventory[0];
        sol.cost_inventory_customers = cost_inventory[1..].iter().sum();

        // total cost
        sol.cost_total =
            sol.cost_transportation + sol.cost_inventory_depot + sol.cost_inventory_customers;

        sol
    }

    pub fn route(&self, t: DayId, v: VehicleId) -> &Route {
        &self.schedule.0[t as usize][v as usize]
    }

    pub fn value(&self) -> f64 {
        self.cost_total
    }
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (t, day_schedule) in self.schedule.0.iter().enumerate() {
            writeln!(f, "Day {}", t + 1)?;
            for (route_idx, route) in day_schedule.iter().enumerate() {
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
        writeln!(f, "{}", self.time)?;

        Ok(())
    }
}

pub struct SolutionPool<'a> {
    solutions: Vec<Solution>,
    idx_best: usize,
    idx_worst: usize,
    time_init: time::Instant,
    cpu: &'a str,
}

impl<'a> SolutionPool<'a> {
    pub fn new(capacity: usize, cpu: &'a str) -> Self {
        SolutionPool {
            solutions: Vec::with_capacity(capacity + 1),
            idx_best: 0,
            idx_worst: 0,
            time_init: time::Instant::now(),
            cpu,
        }
    }

    /// Attempts to add a new solution to the pool and return it if it really has been added
    pub fn add(&mut self, problem: &Problem, schedule: Schedule) -> Option<&Solution> {
        let solution = Solution::new(problem, schedule, self.elapsed_seconds(), self.cpu);

        let len = self.solutions.len();
        if len == self.solutions.capacity() - 1 {
            if solution.value() < self.solutions[self.idx_worst].value() {
                // overwrite worst
                let idx_new = self.idx_worst;
                self.solutions[idx_new] = solution;

                // find new worst index
                self.idx_worst = self
                    .solutions
                    .iter()
                    .enumerate()
                    .max_by(|(_, candidate_a), (_, candidate_b)| {
                        candidate_a
                            .value()
                            .partial_cmp(&candidate_b.value())
                            .unwrap()
                    })
                    .unwrap()
                    .0;

                self.update_best(idx_new);

                Some(&self.solutions[idx_new])
            } else {
                None
            }
        } else if len == 0 {
            self.solutions.push(solution);
            Some(&self.solutions[0])
        } else {
            self.solutions.push(solution);
            self.update_best(len);

            if self.solutions[len].value() > self.solutions[self.idx_worst].value() {
                self.idx_worst = len;
            }
            Some(&self.solutions[len])
        }
    }

    fn update_best(&mut self, idx: usize) {
        if self.solutions[idx].value() < self.solutions[self.idx_best].value() {
            eprintln!("{}", self.solutions[idx]);
            self.idx_best = idx;
        }
    }

    pub fn get_best(&self) -> &Solution {
        debug_assert!(!self.solutions.is_empty());

        &self.solutions[self.idx_best]
    }

    fn elapsed_seconds(&self) -> f64 {
        self.time_init.elapsed().as_millis() as f64 * 1e-3
    }
}
