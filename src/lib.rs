pub mod build_components;
pub mod cluster;
pub mod common;
pub mod concat;
pub mod graph;
pub mod types;

use crate::{
    common::{default_config, open_for_write},
    types::{
        ClusterArgs, ClusterResult, ComponentsArgs, ConcatArgs, ConfigArgs,
        RmBlastOutput, RunArgs,
    },
};
use anyhow::{bail, Result};
use log::debug;

// --------------------------------------------------
pub fn run_sculu(args: &RunArgs, num_threads: usize) -> Result<()> {
    let intermediate_dir = args.outdir.join("intermediate");
    let built_components = build_components::build(
        &ComponentsArgs {
            alphabet: args.alphabet.clone(),
            consensus: args.consensus.clone(),
            alignments: args.alignments.clone(),
            outdir: intermediate_dir.clone(),
            config: args.config.clone(),
        },
        num_threads,
    )?;
    debug!("built_components =\n{built_components:?}");

    let mut clusters: Vec<ClusterResult> = vec![];
    for component in built_components.components {
        let res = cluster::cluster_component(
            &ClusterArgs {
                alphabet: args.alphabet.clone(),
                consensus: built_components.consensus_path.clone(),
                alignments: built_components.alignments.clone(),
                instances: built_components.instances_dir.clone(),
                outdir: intermediate_dir.clone(),
                config: args.config.clone(),
                component,
            },
            num_threads,
        )?;

        clusters.push(res);
    }
    debug!("clusters =\n{clusters:#?}");

    concat::concat_files(&ConcatArgs {
        consensus_path: built_components.consensus_path,
        singletons: built_components.singletons.clone(),
        components: clusters
            .iter()
            .map(|v| v.consensus_path.clone())
            .collect::<Vec<_>>(),
        outdir: args.outdir.clone(),
    })?;

    Ok(())
}

// --------------------------------------------------
pub fn write_config(args: &ConfigArgs) -> Result<()> {
    if args.outfile.exists() {
        bail!(r#"Will not overwrite "{}""#, args.outfile.display());
    }

    let mut outfile = open_for_write(&args.outfile)?;
    writeln!(&mut outfile, "{}", toml::to_string(&default_config())?)?;
    Ok(())
}
