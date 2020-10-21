use std::error::Error;
use std::fmt;
use std::fs;

pub mod solver;

#[derive(Debug)]
pub struct Position {
    pub x: f64,
    pub y: f64,
}

impl Position {
    pub fn new(x: f64, y: f64) -> Position {
        Position { x: x, y: y }
    }

    fn real_distance(&self, other: &Self) -> f64 {
        let xdist = self.x - other.x;
        let ydist = self.y - other.y;
        f64::sqrt(xdist * xdist + ydist * ydist)
    }

    pub fn distance(&self, other: &Self) -> i32 {
        self.real_distance(&other).round() as i32
    }
}

impl fmt::Display for Position {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", self.x, self.y)
    }
}

#[derive(Debug)]
pub struct Depot {
    pub position: Position,
    pub start_level: i32,
    pub daily_production: i32,
    pub daily_cost: f64,
}

impl fmt::Display for Depot {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Depot at {} w/ initial level {}, daily production {}, daily cost {}",
            self.position, self.start_level, self.daily_production, self.daily_cost
        )
    }
}

#[derive(Debug)]
pub struct Customer {
    pub id: usize,
    pub position: Position,
    pub start_level: i32,
    pub max_level: i32,
    pub min_level: i32,
    pub daily_consumption: i32,
    pub daily_cost: f64,
}

impl fmt::Display for Customer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Customer {} at {} w/ levels {}:{}:{}, daily consumption {}, daily cost {}",
            self.id,
            self.position,
            self.min_level,
            self.start_level,
            self.max_level,
            self.daily_consumption,
            self.daily_cost
        )
    }
}

#[derive(Debug)]
pub struct Problem {
    pub num_customers: usize,
    pub num_days: usize,
    pub num_vehicles: usize,
    pub capacity: i32,
    pub depot: Depot,
    pub customers: Vec<Customer>,
}

impl Problem {
    pub fn read_from_file(filename: String) -> Result<Problem, Box<dyn Error>> {
        let contents = fs::read_to_string(filename)?;
        let mut iter = contents.split_whitespace();

        let num_nodes = iter.next().unwrap().parse::<usize>()?;
        let num_days = iter.next().unwrap().parse::<usize>()?;
        let capacity = iter.next().unwrap().parse::<i32>()?;
        let num_vehicles = iter.next().unwrap().parse::<usize>()?;
        let depot_id = iter.next().unwrap().parse::<i32>()?;
        assert!(depot_id == 0, "Depot identifier is not 0");
        let depot_x = iter.next().unwrap().parse::<f64>()?;
        let depot_y = iter.next().unwrap().parse::<f64>()?;
        let depot_level = iter.next().unwrap().parse::<i32>()?;
        let depot_production = iter.next().unwrap().parse::<i32>()?;
        let depot_cost = iter.next().unwrap().parse::<f64>()?;

        let mut customers = Vec::new();
        let mut index: usize = 0;
        while index < num_nodes - 1 {
            index += 1;
            let c_id = iter.next().unwrap().parse::<usize>()?;
            assert!(c_id == index, "Customer does not have expected index");
            let c_x = iter.next().unwrap().parse::<f64>()?;
            let c_y = iter.next().unwrap().parse::<f64>()?;
            let c_start_level = iter.next().unwrap().parse::<i32>()?;
            let c_max_level = iter.next().unwrap().parse::<i32>()?;
            let c_min_level = iter.next().unwrap().parse::<i32>()?;
            let c_consumption = iter.next().unwrap().parse::<i32>()?;
            let c_cost = iter.next().unwrap().parse::<f64>()?;

            customers.push(Customer {
                id: c_id,
                position: Position::new(c_x, c_y),
                start_level: c_start_level,
                max_level: c_max_level,
                min_level: c_min_level,
                daily_consumption: c_consumption,
                daily_cost: c_cost,
            });
        }
        assert!(iter.next() == None, "There is junk at the end of the file");

        Ok(Problem {
            num_customers: num_nodes - 1,
            num_days: num_days,
            num_vehicles: num_vehicles,
            capacity: capacity,
            depot: Depot {
                position: Position::new(depot_x, depot_y),
                start_level: depot_level,
                daily_production: depot_production,
                daily_cost: depot_cost,
            },
            customers: customers,
        })
    }
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Problem with {} customers for {} days with {} vehicles of capacity {}:\n",
            self.num_customers, self.num_days, self.num_vehicles, self.capacity
        )?;
        write!(f, "    {}\n", self.depot)?;
        Ok(for customer in &self.customers {
            write!(f, "    {}\n", customer)?
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    mod position {
        use super::*;

        #[test]
        fn computes_zero_distance() {
            let pos1 = Position::new(12.5, 23.0);
            let pos2 = Position::new(12.5, 23.0);
            assert_eq!(pos1.distance(&pos2), 0);
            assert_eq!(pos2.distance(&pos1), 0);
        }

        #[test]
        fn computes_horizontal_distance() {
            let pos1 = Position::new(12.5, 23.0);
            let pos2 = Position::new(22.5, 23.0);
            assert_eq!(pos1.distance(&pos2), 10);
            assert_eq!(pos2.distance(&pos1), 10);
        }

        #[test]
        fn computes_vertical_distance() {
            let pos1 = Position::new(12.5, 23.0);
            let pos2 = Position::new(12.5, 13.0);
            assert_eq!(pos1.distance(&pos2), 10);
            assert_eq!(pos2.distance(&pos1), 10);
        }

        #[test]
        fn computes_diagonal_distance_rounded_down() {
            let pos1 = Position::new(12.5, 23.0);
            let pos2 = Position::new(13.5, 24.0);
            assert_eq!(pos1.distance(&pos2), 1);
            assert_eq!(pos2.distance(&pos1), 1);
        }

        #[test]
        fn computes_diagonal_distance_rounded_up() {
            let pos1 = Position::new(12.5, 23.0);
            let pos2 = Position::new(14.0, 24.5);
            assert_eq!(pos1.distance(&pos2), 2);
            assert_eq!(pos2.distance(&pos1), 2);
        }
    }
}
