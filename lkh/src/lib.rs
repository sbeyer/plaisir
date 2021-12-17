use lkh_sys::NodeCoords as Coords;
use std::os::raw::c_int;

pub fn run(coords: &[(usize, f64, f64)]) -> Vec<usize> {
    // dummy code that will be replaced
    let data: Vec<Coords> = coords
        .iter()
        .map(|(_, x, y)| Coords { x: *x, y: *y })
        .collect();

    let tour: *const c_int = unsafe { lkh_sys::run(coords.len() as c_int, data.as_ptr()) };
    // Note that ownership of tour is kept by LKH (it's a global variable there)

    (0..coords.len())
        .map(|i| coords[unsafe { *tour.add(i) } as usize - 1].0)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
