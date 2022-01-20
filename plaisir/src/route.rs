use crate::delivery::Deliveries;
use crate::problem::*;

pub struct Solver {}

impl Solver {
    pub fn solve(problem: &Problem, deliveries: &Deliveries, t: DayId, v: VehicleId) -> Vec<usize> {
        let mut visited_sites = deliveries.get_all_delivered_customers(t, v);

        if visited_sites.is_empty() {
            vec![0usize]
        } else {
            visited_sites.push(0); // add depot

            if visited_sites.len() > 3 {
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
                tsp_tour.rotate_left(depot_position);

                tsp_tour
            } else {
                let mut tsp_tour = visited_sites
                    .iter()
                    .map(|site| *site as usize)
                    .collect::<Vec<_>>();
                tsp_tour.sort_unstable();

                tsp_tour
            }
        }
    }
}
