use anyhow::Result;
use assert_cmd::Command;
use predicates::prelude::*;
use tempfile::tempdir;

const PRG: &str = "sculu";

struct RunArgs<'a> {
    consensi: &'a str,
    instances: &'a str,
}

// --------------------------------------------------
#[test]
fn usage() -> Result<()> {
    for flag in &["-h", "--help"] {
        Command::cargo_bin(PRG)?
            .arg(flag)
            .assert()
            .stdout(predicate::str::contains("Usage"));
    }
    Ok(())
}

// --------------------------------------------------
// TODO: get MSA/STK to replace instances
//#[test]
//fn run1() -> Result<()> {
//    run(RunArgs {
//        consensus: "tests/inputs/consensi.fa",
//        instances: "tests/inputs/instances",
//    })
//}

// --------------------------------------------------
fn run(args: RunArgs) -> Result<()> {
    let outdir = tempdir()?;
    let outdir_name = outdir.path().to_string_lossy();
    let outfile = outdir.path().join("final.fa");
    let outname = outfile.to_string_lossy().to_string();

    let args = vec![
        "run",
        "--alphabet",
        "dna",
        "--consensi",
        args.consensi,
        "--instances",
        args.instances,
        "--config",
        "tests/inputs/sculu.toml",
        "--outfile",
        &outname,
        "--outdir",
        &outdir_name,
    ];

    dbg!(Command::cargo_bin(PRG)?.args(args.clone()));

    let output = Command::cargo_bin(PRG)?.args(args).output().unwrap();
    dbg!(&output);
    assert!(output.status.success());

    // Ensure the output directory was created
    assert!(outdir.path().is_dir());

    // Ensure the final consensi file was created
    assert!(outfile.exists());

    Ok(())
}
