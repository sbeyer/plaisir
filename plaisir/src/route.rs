use crate::delivery::Deliveries;
use crate::problem::*;

pub struct Solver {}

impl Solver {
    pub fn solve(problem: &Problem, deliveries: &Deliveries, t: DayId, v: VehicleId) -> Vec<usize> {
        let mut visited_sites = deliveries.get_all_delivered_customers(t, v);

        match visited_sites.len() {
            0 => vec![],
            1 => vec![visited_sites[0] as usize],
            2 => vec![visited_sites[0] as usize, visited_sites[1] as usize],
            _ => {
                visited_sites.push(0); // add depot

                let tsp_instance = visited_sites
                    .iter()
                    .map(|site| {
                        let site = problem.site(*site);
                        let pos = &site.position();
                        (site.id() as usize, pos.x, pos.y)
                    })
                    .collect::<Vec<_>>();
                let mut tsp_tour = lkh::run(&tsp_instance);

                let depot_position = tsp_tour
                    .iter()
                    .position(|site| *site == 0)
                    .expect("Depot is expected to be in TSP tour");
                tsp_tour.rotate_left(depot_position + 1);
                tsp_tour.pop();

                tsp_tour
            }
        }
    }
}
