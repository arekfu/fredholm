use structopt::StructOpt;

/// Parameters for the Monte Carlo simulation
#[derive(StructOpt)]
pub struct Params {
    /// the coefficient for the linear term in the kernel
    #[structopt(long, default_value = "1.0")]
    pub alpha: f64,

    /// the cross section
    #[structopt(long, default_value = "1.0")]
    pub sigma: f64,

    /// the absorption probability
    #[structopt(long, default_value = "0.9")]
    pub lambda: f64,
}

/// Solve a Fredholm equation of the second kind with Monte Carlo
#[derive(StructOpt)]
pub struct FredholmConfig {
    /// the number of replicas
    #[structopt(short, long, default_value = "1000")]
    pub replicas: usize,

    /// the maximum packet size
    #[structopt(short, long, default_value = "1000")]
    pub packet_size: usize,

    /// the cutoff for the coordinate
    #[structopt(long, default_value = "10.0")]
    pub cutoff: f64,

    /// the Russian roulette threshold
    #[structopt(long, default_value = "0.5")]
    pub roulette_threshold: f64,

    /// the number of bins for the output histogram
    #[structopt(long, default_value = "100")]
    pub bins: usize,

    /// the output file
    #[structopt(short, long, default_value = "fredholm.dat")]
    pub output: String,

    /// debug flag
    #[structopt(short, long)]
    pub debug: bool,

    /// number of threads
    #[structopt(short, long, default_value = "1")]
    pub threads: usize,

    #[structopt(flatten)]
    pub params: Params,
}

