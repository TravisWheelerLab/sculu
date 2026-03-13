use anyhow::{anyhow, bail, Result};
use clap::Parser;
use sculu::{
    common::format_seconds,
    types::{Cli, Command},
};
use std::{
    fs::{self, File},
    io::BufWriter,
    path::PathBuf,
    time::Instant,
};

// --------------------------------------------------
fn main() {
    if let Err(e) = run(Cli::parse()) {
        eprintln!("Error: {e}");
        std::process::exit(1);
    }
}

// --------------------------------------------------
fn run(args: Cli) -> Result<()> {
    let log_file = args.logfile.unwrap_or(PathBuf::from("debug.log"));
    let num_threads = args.num_threads.unwrap_or(num_cpus::get());

    if let Some(dir) = log_file.parent() {
        if !dir.exists() {
            fs::create_dir_all(dir)?;
        }
    }

    env_logger::Builder::new()
        .filter_level(log::LevelFilter::Debug)
        .target(match log_file.to_str() {
            Some("-") => env_logger::Target::Stdout,
            Some(path) => env_logger::Target::Pipe(Box::new(BufWriter::new(
                File::create(&log_file).map_err(|e| anyhow!("{path}: {e}"))?,
            ))),
            _ => bail!(r#"Cannot parse --log-file "{}""#, log_file.display()),
        })
        .init();

    let start = Instant::now();
    match &args.command.expect("Missing command") {
        Command::Config(args) => {
            sculu::write_config(args)?;
            println!(r#"Wrote config file to "{}""#, args.outfile.display());
            Ok(())
        }
        Command::Components(args) => {
            sculu::build_components::build(args, num_threads)?;
            println!(
                r#"Wrote output to "{}" in {}."#,
                args.outdir.display(),
                format_seconds(start.elapsed().as_secs()),
            );
            Ok(())
        }
        Command::Cluster(args) => {
            sculu::cluster::cluster_component(args, num_threads)?;
            println!(
                r#"Wrote output to "{}" in {}."#,
                args.outdir.display(),
                format_seconds(start.elapsed().as_secs()),
            );
            Ok(())
        }
        Command::Concat(args) => {
            sculu::concat_files(args)?;
            println!(
                r#"Wrote output to "{}" in {}."#,
                args.outfile.display(),
                format_seconds(start.elapsed().as_secs()),
            );
            Ok(())
        }
        Command::Run(args) => {
            // Do all the things!
            sculu::run_sculu(args, num_threads)?;
            println!(
                r#"Wrote output to "{}" in {}."#,
                args.outdir.display(),
                format_seconds(start.elapsed().as_secs()),
            );
            Ok(())
        }
    }
}
