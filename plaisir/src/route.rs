use crate::delivery::{Deliveries, Delivery};
use crate::problem::*;

type Route = Vec<SiteId>;

struct RouteTrie {
    root: RouteTrieNode,
    number: usize,
}

impl RouteTrie {
    fn new() -> Self {
        Self {
            root: RouteTrieNode::new(0),
            number: 0,
        }
    }

    fn set(&mut self, customer_set: &[SiteId], tour: Vec<usize>) {
        //-> &Route {
        self.number += 1;
        self.root.set(customer_set, tour);
    }

    fn get(&self, customer_set: &[SiteId]) -> Option<&Route> {
        self.root.get(customer_set)
    }
}

struct RouteTrieNode {
    id: SiteId,
    next: Vec<RouteTrieNode>,
    route: Option<Route>,
}

impl RouteTrieNode {
    fn new(id: SiteId) -> Self {
        Self {
            id,
            next: Vec::new(),
            route: None,
        }
    }

    fn set(&mut self, customer_set: &[SiteId], tour: Vec<usize>) {
        // -> &Route { <- TODO, but Rust's borrow checker doesn't like it...
        if customer_set.is_empty() {
            self.route = Some(tour.into_iter().map(|x| x as SiteId).collect());
            //self.route.as_ref().unwrap() // TODO
        } else {
            let opt_next = self
                .next
                .iter_mut()
                .find(|candidate| candidate.id == customer_set[0]);

            if let Some(next) = opt_next {
                next.set(&customer_set[1..], tour); // TODO without semicolon
            } else {
                let next = Self::new(customer_set[0]);
                self.next.push(next);
                let last = self.next.len() - 1;
                self.next[last].set(&customer_set[1..], tour); // TODO without semicolon
            }
        }
    }

    fn get(&self, customer_set: &[SiteId]) -> Option<&Route> {
        if customer_set.is_empty() {
            if let Some(route) = &self.route {
                Some(route)
            } else {
                None
            }
        } else {
            let opt_next = self
                .next
                .iter()
                .find(|candidate| candidate.id == customer_set[0]);

            if let Some(next) = opt_next {
                next.get(&customer_set[1..])
            } else {
                None
            }
        }
    }
}

pub struct Solver {
    saved: RouteTrie,

    // Number of nontrivial calls
    ncalls: usize,
}

impl Solver {
    pub fn new() -> Self {
        Self {
            saved: RouteTrie::new(),
            ncalls: 0,
        }
    }

    pub fn solve(
        &mut self,
        problem: &Problem,
        deliveries: &Deliveries,
        t: DayId,
        v: VehicleId,
    ) -> Vec<Delivery> {
        let visited_customers = deliveries.get_all_delivered_customers(t, v);

        if visited_customers.len() < 3 {
            &visited_customers
        } else {
            self.ncalls += 1;

            if self.ncalls % 2500 == 0 {
                eprintln!("# Route solver called (with non-trivial routes) {} times and {} routes are saved", self.ncalls, self.saved.number)
            }

            if let Some(saved_route) = self.saved.get(&visited_customers) {
                saved_route
            } else {
                let mut tsp_instance = visited_customers
                    .iter()
                    .map(|site| {
                        let site = problem.site(*site);
                        let pos = &site.position();
                        (site.id() as usize, pos.x, pos.y)
                    })
                    .collect::<Vec<_>>();
                let depot_pos = problem.site(0).position();
                tsp_instance.push((0, depot_pos.x, depot_pos.y));

                let mut tsp_tour = lkh::run(&tsp_instance);

                let depot_position = tsp_tour
                    .iter()
                    .position(|site| *site == 0)
                    .expect("Depot is expected to be in TSP tour");
                tsp_tour.rotate_left(depot_position + 1);
                tsp_tour.pop();

                self.saved.set(&visited_customers, tsp_tour);
                self.saved.get(&visited_customers).unwrap()
            }
        }
            .iter()
            .map(|site| Delivery {
                quantity: deliveries.get(t, v, *site),
                customer: *site,
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod route_tries {
        use super::*;

        #[test]
        fn set_and_get() {
            let mut trie = RouteTrie::new();
            assert_eq!(trie.number, 0);

            let k_2_3_5 = &[2, 3, 5];
            let t_2_3_5 = vec![0, 5, 2, 3];
            let e_2_3_5 = t_2_3_5.iter().map(|x| *x as SiteId).collect::<Vec<_>>();
            trie.set(k_2_3_5, t_2_3_5);
            assert_eq!(trie.number, 1);

            let k_1 = &[1];
            let t_1 = vec![0, 1];
            let e_1 = t_1.iter().map(|x| *x as SiteId).collect::<Vec<_>>();
            trie.set(k_1, t_1);
            assert_eq!(trie.number, 2);

            let k_2 = &[2];
            let t_2 = vec![0, 2];
            let e_2 = t_2.iter().map(|x| *x as SiteId).collect::<Vec<_>>();
            trie.set(k_2, t_2);
            assert_eq!(trie.number, 3);

            let k_1_3_5 = &[1, 3, 5];
            let t_1_3_5 = vec![0, 1, 5, 3];
            let e_1_3_5 = t_1_3_5.iter().map(|x| *x as SiteId).collect::<Vec<_>>();
            trie.set(k_1_3_5, t_1_3_5);
            assert_eq!(trie.number, 4);

            let a_1 = trie.get(k_1);
            assert!(a_1.is_some());
            assert_eq!(*a_1.unwrap(), e_1);

            let a_2 = trie.get(k_2);
            assert!(a_2.is_some());
            assert_eq!(*a_2.unwrap(), e_2);

            let a_1_3_5 = trie.get(k_1_3_5);
            assert!(a_1_3_5.is_some());
            assert_eq!(*a_1_3_5.unwrap(), e_1_3_5);

            let a_2_3_5 = trie.get(k_2_3_5);
            assert!(a_2_3_5.is_some());
            assert_eq!(*a_2_3_5.unwrap(), e_2_3_5);
        }
    }
}
