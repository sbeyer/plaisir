use lkh_sys::NodeCoords as Coords;
use std::os::raw::c_int;

pub fn run(coords: &[(usize, f64, f64)]) -> Vec<usize> {
    if coords.len() <= 3 {
        // handle trivial cases that lkh_sys will not handle
        coords.iter().map(|(id, _, _)| *id).collect()
    } else {
        // convert data for lkh_sys::run()
        let data: Vec<Coords> = coords
            .iter()
            .map(|(_, x, y)| Coords { x: *x, y: *y })
            .collect();

        // invoke lkh_sys::run()
        let tour: *const c_int = unsafe { lkh_sys::run(coords.len() as c_int, data.as_ptr()) };
        // Note that ownership of tour is kept by LKH (it's a global variable there)

        // prepare return vector
        (0..coords.len())
            .map(|i| coords[unsafe { *tour.add(i) } as usize - 1].0)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn verify(coords: &[(usize, f64, f64)]) {
        let tour = run(&coords);
        assert_eq!(tour.len(), coords.len());

        for (site, _, _) in coords.iter() {
            assert!(
                tour.contains(site),
                "tour {:?} does not contain site {}",
                tour,
                site
            );
        }
    }

    #[test]
    fn run_on_empty_input() {
        let coords = vec![];
        verify(&coords);
    }

    #[test]
    fn run_on_single_site() {
        let coords = vec![(512, 1.34, 3.14)];
        verify(&coords);
    }

    #[test]
    fn run_on_two_sites() {
        let coords = vec![(512, 1.34, 3.14), (256, 4.13, 3.41)];
        verify(&coords);
    }

    #[test]
    fn run_on_three_sites() {
        let coords = vec![(512, 1.34, 3.14), (256, 4.13, 3.41), (17, 4.31, 1.43)];
        verify(&coords);
    }

    #[test]
    fn run_on_many_sites() {
        let coords = vec![
            (12, 0.0, 0.0),
            (42, 1.0, 0.0),
            (7, 1.0, 1.0),
            (113, 0.5, 1.5),
            (2, 0.0, 2.0),
            (13, 0.0, 1.0),
            (29, 0.0, 0.5),
        ];
        verify(&coords);
    }
}
