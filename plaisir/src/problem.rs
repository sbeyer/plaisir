use std::error::Error;
use std::fmt;
use std::fs;

#[derive(Debug)]
pub struct Position {
    pub x: f64,
    pub y: f64,
}

impl Position {
    pub fn new(x: f64, y: f64) -> Position {
        Position { x, y }
    }

    fn real_distance(&self, other: &Self) -> f64 {
        let xdist = self.x - other.x;
        let ydist = self.y - other.y;
        f64::sqrt(xdist * xdist + ydist * ydist)
    }

    pub fn distance(&self, other: &Self) -> i32 {
        self.real_distance(other).round() as i32
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

pub type SiteId = u8;
pub type DayId = u8;
pub type VehicleId = u8;
pub type Load = usize;

#[derive(Debug)]
pub struct Customer {
    pub id: SiteId,
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

pub trait Site {
    /// Returns the numeric identifer
    fn id(&self) -> SiteId;

    /// Returns the start inventory level
    fn level_start(&self) -> f64;

    /// Returns the daily inventory level change
    fn level_change(&self) -> f64;

    /// Returns minimum and maximum inventory levels
    fn level_bounds(&self) -> (f64, f64);

    /// Returns the daily cost for each inventory unit
    fn cost(&self) -> f64;

    /// Returns the position
    fn position(&self) -> &Position;
}

impl Site for Depot {
    fn id(&self) -> SiteId {
        0
    }

    fn level_start(&self) -> f64 {
        self.start_level.into()
    }

    fn level_change(&self) -> f64 {
        self.daily_production.into()
    }

    fn level_bounds(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }

    fn cost(&self) -> f64 {
        self.daily_cost
    }

    fn position(&self) -> &Position {
        &self.position
    }
}

impl Site for Customer {
    fn id(&self) -> SiteId {
        self.id
    }

    fn level_start(&self) -> f64 {
        self.start_level.into()
    }

    fn level_change(&self) -> f64 {
        (-self.daily_consumption).into()
    }

    fn level_bounds(&self) -> (f64, f64) {
        (self.min_level.into(), self.max_level.into())
    }

    fn cost(&self) -> f64 {
        self.daily_cost
    }

    fn position(&self) -> &Position {
        &self.position
    }
}

#[derive(Debug)]
pub struct Problem {
    pub num_sites: usize,
    pub num_customers: usize,
    pub num_days: usize,
    pub num_vehicles: usize,
    pub capacity: Load,
    pub depot: Depot,
    pub customers: Vec<Customer>,
}

impl Problem {
    pub fn read_from_file(filename: String) -> Result<Problem, Box<dyn Error>> {
        let contents = fs::read_to_string(filename)?;
        let mut iter = contents.split_whitespace();

        let num_sites = iter.next().unwrap().parse::<usize>()?;
        let num_days = iter.next().unwrap().parse::<usize>()?;
        let capacity = iter.next().unwrap().parse::<Load>()?;
        let num_vehicles = iter.next().unwrap().parse::<usize>()?;
        let depot_id = iter.next().unwrap().parse::<usize>()?;
        assert!(depot_id == 0, "Depot identifier is not 0");
        let depot_x = iter.next().unwrap().parse::<f64>()?;
        let depot_y = iter.next().unwrap().parse::<f64>()?;
        let depot_level = iter.next().unwrap().parse::<i32>()?;
        let depot_production = iter.next().unwrap().parse::<i32>()?;
        let depot_cost = iter.next().unwrap().parse::<f64>()?;

        let mut customers = Vec::new();
        let mut index: usize = 0;
        while index < num_sites - 1 {
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
                id: c_id as SiteId,
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
            num_sites,
            num_customers: num_sites - 1,
            num_days,
            num_vehicles,
            capacity,
            depot: Depot {
                position: Position::new(depot_x, depot_y),
                start_level: depot_level,
                daily_production: depot_production,
                daily_cost: depot_cost,
            },
            customers,
        })
    }

    pub fn site(&self, index: SiteId) -> &dyn Site {
        if index == 0 {
            &self.depot
        } else {
            &self.customers[index as usize - 1]
        }
    }

    pub fn distance(&self, i: SiteId, j: SiteId) -> i32 {
        if i == j {
            0
        } else {
            let pos1 = &self.site(i).position();
            let pos2 = &self.site(j).position();
            pos1.distance(pos2)
        }
    }

    pub fn all_days(&self) -> impl Iterator<Item = DayId> {
        0..(self.num_days as DayId)
    }

    pub fn all_vehicles(&self) -> impl Iterator<Item = VehicleId> {
        0..(self.num_vehicles as VehicleId)
    }

    pub fn all_sites(&self) -> impl Iterator<Item = SiteId> {
        0..(self.num_sites as SiteId)
    }

    pub fn all_sites_except(&self, exception: SiteId) -> impl Iterator<Item = SiteId> {
        self.all_sites().filter(move |i| *i != exception)
    }

    pub fn all_sites_after(&self, after: SiteId) -> impl Iterator<Item = SiteId> {
        after + 1..(self.num_sites as SiteId)
    }

    pub fn all_customers(&self) -> impl Iterator<Item = SiteId> {
        1..(self.num_sites as SiteId)
    }
}

impl fmt::Display for Problem {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Problem with {} customers for {} days with {} vehicles of capacity {}:",
            self.num_customers, self.num_days, self.num_vehicles, self.capacity
        )?;
        writeln!(f, "    {}", self.depot)?;
        for customer in &self.customers {
            writeln!(f, "    {}", customer)?
        }
        Ok(())
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
