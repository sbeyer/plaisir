use crate::delivery::Delivery;
use crate::problem::{DayId, Problem, SiteId, VehicleId};
use std::fmt;

pub type Route = Vec<Delivery>;
pub type Schedule = Vec<Vec<Route>>;

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
            schedule: Vec::new(),
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
        for day_schedule in sol.schedule.iter() {
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

        for day_schedule in sol.schedule.iter() {
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
        &self.schedule[t as usize][v as usize]
    }

    pub fn value(&self) -> f64 {
        self.cost_total
    }
}

impl fmt::Display for Solution {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (t, day_schedule) in self.schedule.iter().enumerate() {
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
