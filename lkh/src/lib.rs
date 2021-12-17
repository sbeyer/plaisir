use lkh_sys::NodeCoords as Coords;
use std::os::raw::c_int;

pub fn run(dimension: usize) -> Vec<usize> {
    // dummy code that will be replaced
    let mut data = Vec::<Coords>::new();
    for i in 0..dimension {
        data.push(Coords {
            x: 0.0,
            y: i as u64 as f64,
        });
    }

    let tour: *const c_int = unsafe { lkh_sys::run(dimension as c_int, data.as_ptr()) };
    // Note that ownership of tour is kept by LKH (it's a global variable there)

    (0..dimension)
        .map(|i| unsafe { *tour.add(i) } as usize - 1)
        .collect()
}

// use cargo test -- --nocapture to see output
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_on_ten_sites() {
        let tour = run(10);
        println!("Resulting tour:");
        for site in tour.iter() {
            println!(" * {}", site);
        }
        assert_eq!(tour.len(), 10);
    }
}
