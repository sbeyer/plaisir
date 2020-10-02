use std::error::Error;
use std::fmt;
use std::fs;

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    println!("Reading problem from file {:?}", config.filename);

    let problem = Problem::read_from_file(config.filename)?;
    println!("{}", problem);

    Ok(())
}

#[derive(Debug)]
pub struct Position {
    pub x: f64,
    pub y: f64,
}

impl Position {
    pub fn new(x: f64, y: f64) -> Position {
        Position { x: x, y: y }
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
    pub id: i32,
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
    pub num_nodes: i32,
    pub num_days: i32,
    pub num_vehicles: i32,
    pub capacity: i32,
    pub depot: Depot,
    pub customers: Vec<Customer>,
}

impl Problem {
    pub fn read_from_file(filename: String) -> Result<Problem, Box<dyn Error>> {
        let contents = fs::read_to_string(filename)?;
        let mut iter = contents.split_whitespace();

        let num_nodes = iter.next().unwrap().parse::<i32>()?;
        let num_days = iter.next().unwrap().parse::<i32>()?;
        let capacity = iter.next().unwrap().parse::<i32>()?;
        let num_vehicles = iter.next().unwrap().parse::<i32>()?;
        let depot_id = iter.next().unwrap().parse::<i32>()?;
        assert!(depot_id == 0, "Depot identifier is not 0");
        let depot_x = iter.next().unwrap().parse::<f64>()?;
        let depot_y = iter.next().unwrap().parse::<f64>()?;
        let depot_level = iter.next().unwrap().parse::<i32>()?;
        let depot_production = iter.next().unwrap().parse::<i32>()?;
        let depot_cost = iter.next().unwrap().parse::<f64>()?;

        let mut customers = Vec::new();
        let mut index = 0;
        while index < num_nodes - 1 {
            index += 1;
            let c_id = iter.next().unwrap().parse::<i32>()?;
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
            num_nodes: num_nodes,
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
            "Problem with {} locations for {} days with {} vehicles of capacity {}:\n",
            self.num_nodes, self.num_days, self.num_vehicles, self.capacity
        )?;
        write!(f, "    {}\n", self.depot)?;
        Ok(for customer in &self.customers {
            write!(f, "    {}\n", customer)?
        })
    }
}

pub struct Config {
    pub filename: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() != 2 {
            return Err("need filename argument");
        }

        let filename = args[1].clone();

        Ok(Config { filename })
    }
}
