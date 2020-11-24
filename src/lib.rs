use std::error::Error;

mod irp;

pub fn run(config: Config) -> Result<(), Box<dyn Error>> {
    println!("Reading problem from file {:?}", config.filename);

    let problem = irp::Problem::read_from_file(config.filename)?;
    println!("{}", problem);

    irp::solver::solve(problem);

    Ok(())
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
