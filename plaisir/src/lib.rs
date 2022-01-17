use std::error::Error;

mod delivery;
mod heuristic;
mod problem;
mod solution;
mod solver;

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    println!("Reading problem from file {:?}", config.filename);

    let problem = problem::Problem::read_from_file(config.filename)?;
    println!("{}", problem);

    solver::solve(problem, config.cpu);

    Ok(())
}

pub struct Config {
    pub filename: String,
    pub cpu: String,
}

impl Config {
    pub fn new(args: &[String]) -> Result<Config, &'static str> {
        if args.len() != 2 {
            return Err("need filename argument");
        }

        let filename = args[1].clone();
        let cpu = "Intel Core i5-10210U @ 1.60GHz".to_string(); // hardcoded for now

        Ok(Config { filename, cpu })
    }
}
