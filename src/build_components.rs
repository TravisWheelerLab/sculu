use crate::{
    common::{
        copy_fasta, get_config, open, open_for_write, read_instances_dir, run_blastp,
        run_rmblastn,
    },
    graph,
    types::{
        AlignedConsensusPair, BuiltComponents, CheckFamilyInstancesArgs,
        CheckedFamilyResult, Components, ComponentsArgs, Config, Direction, Partition,
        RmBlastOutput, SelectInstancesArgs, SequenceAlphabet, StringPair,
    },
};
use anyhow::{anyhow, bail, Result};
use itertools::Itertools;
use kseq::parse_reader;
use log::debug;
use noodles_fasta::{self, io::Reader as FastaReader, io::Writer as FastaWriter};
use rayon::prelude::*;
use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    fs,
    io::{BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};

// --------------------------------------------------
pub fn build(args: &ComponentsArgs, num_threads: usize) -> Result<BuiltComponents> {
    let config = get_config(&args.config)?;

    if !&args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    let taken_instances_dir = args.outdir.join("instances");
    fs::create_dir_all(&taken_instances_dir)?;

    let checked = check_family_instances(CheckFamilyInstancesArgs {
        consensus_path: &args.consensus,
        instances_dir: &args.instances,
        out_dir: &args.outdir,
        taken_instances_dir: &taken_instances_dir,
        config: &config,
        alphabet: &args.alphabet,
        num_threads,
    })?;

    // This will only BLAST if there are no components from a previous run
    let components = align_consensus_to_self(
        &checked.consensus_path,
        &args.outdir,
        &args.alphabet,
        &config,
        num_threads,
    )?;
    let num_components = components.components.len();
    debug!(
        "Built {num_components} component{}",
        if num_components == 1 { "" } else { "s" }
    );

    Ok(BuiltComponents {
        singletons: components.singletons,
        components: components.components,
        alignments: components.alignments,
        consensus_path: checked.consensus_path,
        instances_dir: taken_instances_dir,
    })
}

// --------------------------------------------------
// This function will first filter the instances for those that
// are most representative of the family, and then it will filter
// the consensus for those that have sufficient supporting instances.
fn check_family_instances(
    args: CheckFamilyInstancesArgs,
) -> Result<CheckedFamilyResult> {
    debug!(
        r#"Checking consensus file "{}""#,
        args.consensus_path.display()
    );

    let instances = read_instances_dir(&args.instances_dir)?;

    debug!(
        "Found {} instance files in '{}'",
        instances.len(),
        args.instances_dir.display(),
    );

    let working_dir = args.out_dir.join("select");
    fs::create_dir_all(&working_dir)?;

    let min_num_instances = match args.alphabet {
        SequenceAlphabet::Dna => args.config.general.min_num_instances_dna,
        _ => args.config.general.min_num_instances_prot,
    };

    let now = Instant::now();
    instances.par_iter().try_for_each(
        |(family_name, instance_path)| -> Result<()> {
            let taken_path = args.taken_instances_dir.join(format!("{family_name}.fa"));

            match select_instances(SelectInstancesArgs {
                consensus_path: args.consensus_path,
                family_name: &family_name,
                from_path: &instance_path,
                to_path: &taken_path,
                working_dir: &working_dir,
                alphabet: &args.alphabet,
                config: &args.config,
                num_threads: args.num_threads,
            }) {
                Err(e) => eprintln!("Warning: {e}"),
                Ok(num_taken) => {
                    debug!("Took {num_taken} instances for {family_name}");
                }
            }
            Ok(())
        },
    )?;

    // Remove families with too few instances
    let too_few_path = args.out_dir.join("too_few.fa");
    let mut too_few_writer =
        FastaWriter::new(BufWriter::new(open_for_write(&too_few_path)?));
    let taken_instances = read_instances_dir(&args.taken_instances_dir)?;
    for (family_name, instance_path) in taken_instances {
        let mut reader = FastaReader::new(open(&instance_path)?);
        let num_taken = reader.records().count();
        if num_taken < min_num_instances {
            debug!(
                r#"Removing family "{family_name}", too few instances ({num_taken})"#,
            );
            let num_taken = copy_fasta(
                &[family_name.clone()],
                args.consensus_path,
                &mut too_few_writer,
            )?;
            if num_taken != 1 {
                bail!(
                    r#"Failed to copy family "{family_name}" to {}"#,
                    too_few_path.display()
                );
            }
            fs::remove_file(instance_path)?;
        }
    }

    let family_to_instance = read_instances_dir(args.taken_instances_dir)?;
    debug!(
        "Finished instance selection for {} consensus in {:?}",
        family_to_instance.len(),
        now.elapsed()
    );

    // Create a new consensus file containing only the families that have instances
    let taken_consensus_path = args.out_dir.join("consensus.fa");
    let mut out_consensus = open_for_write(&taken_consensus_path)?;
    debug!(
        r#"Creating new consensus file "{}""#,
        taken_consensus_path.display()
    );

    let mut reader = parse_reader(open(args.consensus_path)?)?;
    let mut consensus_names: HashMap<String, u32> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        let family = rec.head().to_string();
        if family_to_instance.contains_key(&family) {
            debug!(r#"... adding "{family}""#);
            consensus_names
                .entry(rec.head().to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
            writeln!(out_consensus, ">{family}\n{}", rec.seq())?;
        } else {
            debug!(r#"... skipping "{family}""#);
        }
    }

    // This doesn't seem possible
    //let mut dups: Vec<_> = consensus_names
    //    .iter()
    //    .flat_map(|(name, &count)| (count > 1).then_some(name))
    //    .collect();

    //if !dups.is_empty() {
    //    // Have to sort for tests
    //    dups.sort();
    //    bail!(
    //        "The following consensus IDs are duplicated: {}",
    //        dups.iter().join(", ")
    //    );
    //}

    Ok(CheckedFamilyResult {
        consensus_path: taken_consensus_path,
        family_to_instance,
    })
}

// --------------------------------------------------
// Align the consensus sequences to themselves to identify clusters, e.g.:
// [
//     [ "Charlie13a" ],
//     [ "Tigger3c" ],
//     [ "AluSc", "AluYm1", "AluYh3", "AluYh9", "AluYb8", "AluYb9" ],
//     [ "Charlie1a", "Charlie1", "Charlie2a", "Charlie2b" ],
// ]
// These will be written into a "components" output directory
// A "singletons" file will contain the components with only one member
// The other 1..N will be written to files "component-N"
// Returns the component file paths
//
pub fn align_consensus_to_self(
    consensus_path: &PathBuf,
    out_dir: &Path,
    alphabet: &SequenceAlphabet,
    config: &Config,
    num_threads: usize,
) -> Result<Components> {
    debug!("Aligning consensus sequences to self");

    let components_dir = out_dir.join("components");
    fs::create_dir_all(&components_dir)?;

    let num_components = fs::read_dir(&components_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, components_dir.display()))?
        .filter_map(|entry| entry.ok())
        .collect::<Vec<_>>()
        .len();

    if num_components > 0 {
        debug!("Reusing existing component files");
    } else {
        let blast_dir = out_dir.join("consensus_cluster");
        let blast_out = match alphabet {
            SequenceAlphabet::Dna => run_rmblastn(
                &blast_dir,
                consensus_path,
                consensus_path,
                config,
                num_threads,
            )?,
            _ => run_blastp(
                &blast_dir,
                consensus_path,
                consensus_path,
                config,
                num_threads,
            )?,
        };

        // There may be multiple hits per pair, take highest
        let best_alignments = blast_dir.join("best.tsv");
        take_best_alignments(&blast_out, &best_alignments, config)?;

        // Filter out: Hits to self, insufficient coverage
        let alignment_file = components_dir.join("alignment.tsv");
        let mut alignment_wtr = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(&alignment_file)
            .map_err(|e| {
                anyhow!("Failed to write '{}': {e}", &alignment_file.display())
            })?;

        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&best_alignments)
            .map_err(|e| {
                anyhow!("Failed to read '{}': {e}", &best_alignments.display())
            })?;

        //let mut saw_consensus: HashSet<String> = HashSet::new();
        let mut records = vec![];
        for res in reader.records() {
            let record: RmBlastOutput = res?.deserialize(None)?;
            let query_dir = if record.query_start < record.query_end {
                Direction::Forward
            } else {
                Direction::Reverse
            };
            let subject_dir = if record.subject_start < record.subject_end {
                Direction::Forward
            } else {
                Direction::Reverse
            };
            //let _ = saw_consensus.insert(record.target.to_string());
            //let _ = saw_consensus.insert(record.query.to_string());
            alignment_wtr.serialize(AlignedConsensusPair {
                target: record.target.to_string(),
                query: record.query.clone(),
                is_flipped: query_dir != subject_dir,
            })?;
            records.push(record);
        }

        //debug!("saw_consensus = {saw_consensus:?}");
        //let consensus_names = get_sequence_ids(consensus_path)?;
        //debug!("consensus_names = {consensus_names:?}");
        //let missing_consensus: Vec<_> = consensus_names
        //    .into_iter()
        //    .filter(|name| !saw_consensus.contains(name))
        //    .collect();
        //debug!("missing_consensus = {missing_consensus:?}");

        let mut components = graph::connected_components(records);

        // Sort by size of component group ascending
        components.sort_by_key(|v| v.len());

        debug!("components = {components:#?}");

        let (singles, multis) =
            components.split_at(components.partition_point(|v| v.len() == 1));

        if !singles.is_empty() {
            let singletons_file = components_dir.join("singletons");
            let mut singletons = open_for_write(&singletons_file)?;
            writeln!(singletons, "{}", singles.iter().flatten().join("\n"))?;
        }

        let width = multis.len().to_string().len();
        for (num, component) in multis.iter().enumerate() {
            let component_path =
                components_dir.join(format!("component-{num:0width$}"));
            let mut file = open_for_write(&component_path)?;
            writeln!(file, "{}", component.join("\n"))?;
        }
    }

    let mut singletons: Option<PathBuf> = None;
    let mut alignments: Option<PathBuf> = None;
    let mut components = vec![];
    let entries = fs::read_dir(&components_dir).map_err(|e| {
        anyhow!(r#"Failed to read "{}": {e}"#, components_dir.display())
    })?;

    for entry in entries {
        let entry = entry?;
        let file_name = entry.file_name().to_string_lossy().to_string();
        let path = entry.path().to_path_buf();

        if file_name == "singletons" {
            singletons = Some(path)
        } else if file_name == "alignment.tsv" {
            alignments = Some(path)
        } else if file_name.starts_with("component-") {
            components.push((file_name, path))
        }
    }

    let alignments = alignments.expect("Failed to find 'alignments.tsv'");

    // Sort by name, and the files are named in order of increasing size
    components.sort_by_key(|t| t.0.clone());

    Ok(Components {
        singletons,
        components: components.into_iter().map(|(_name, path)| path).collect(),
        alignments,
    })
}

// --------------------------------------------------
fn select_instances(args: SelectInstancesArgs) -> Result<usize> {
    if args.to_path.is_file() {
        let mut reader = FastaReader::new(open(args.to_path)?);
        let num = reader.records().count();
        debug!("Reusing existing instance file: {}", args.to_path.display());
        return Ok(num);
    }

    let blast_dir = args.working_dir.join(args.family_name);
    fs::create_dir_all(&blast_dir)?;

    // Extract the family's consensus sequence
    let db_path = blast_dir.join("db.fa");
    if !db_path.exists() {
        let mut writer = FastaWriter::new(BufWriter::new(open_for_write(&db_path)?));
        let num_taken = copy_fasta(
            &[args.family_name.to_string()],
            args.consensus_path,
            &mut writer,
        )?;

        if num_taken != 1 {
            fs::remove_file(db_path)?;
            bail!(
                r#"Failed to find family "{}" in consensus "{}""#,
                args.family_name,
                args.consensus_path.display()
            );
        }
    }

    // Get the length of the consensus
    let mut reader = FastaReader::new(open(&db_path)?);
    let consensus_len = reader
        .records()
        .next()
        .unwrap_or_else(|| panic!("failed to read {db_path:?}"))
        .map(|rec| rec.sequence().len())?;

    // BLAST the instances to the consensus
    let blast_out = match args.alphabet {
        SequenceAlphabet::Dna => run_rmblastn(
            &blast_dir,
            &db_path,
            args.from_path,
            args.config,
            args.num_threads,
        )?,
        _ => run_blastp(
            &blast_dir,
            &db_path,
            args.from_path,
            args.config,
            args.num_threads,
        )?,
    };
    let alignments = parse_alignment(&blast_out)?;

    // Remove short alignments, find the highest score for each hit
    let min_len = match args.alphabet {
        SequenceAlphabet::Dna => args.config.general.min_instance_sequence_length_dna,
        _ => args.config.general.min_instance_sequence_length_prot,
    };
    let mut targets: HashMap<String, usize> = HashMap::new();
    for aln in alignments.iter().filter(|aln| aln.query_len >= min_len) {
        if let Some(val) = targets.get_mut(&aln.target) {
            *val = max(aln.score, *val);
        } else {
            targets.insert(aln.target.clone(), aln.score);
        }
    }

    // Create a lookup of the (target, high score)
    let wanted: HashSet<(String, usize)> = targets.into_iter().collect();

    // Use the "wanted" hash to filter the alignments to the single best hit
    let mut filtered: Vec<_> = alignments
        .into_iter()
        .filter(|aln| wanted.contains(&(aln.target.clone(), aln.score)))
        .collect();

    // Sort the hits by CpG-adjusted Kimura divergence ascending
    filtered.sort_by(|a, b| a.cpg_kdiv.partial_cmp(&b.cpg_kdiv).unwrap());

    // Split into top 75%, bottom 25%
    let num = filtered.len();
    let three_quarters = num / 2 + num / 4;
    let (best, worst) = filtered.split_at_mut(three_quarters);

    // Sort by query_len descending
    best.sort_by(|a, b| b.query_len.cmp(&a.query_len));
    worst.sort_by(|a, b| b.query_len.cmp(&a.query_len));

    // Create an array representing 10-bp chunks of the consensus to measure
    // coverage by the instances
    let mut consensus_cov = vec![0; consensus_len.div_ceil(10)];
    let mut i = 0; // Index into "best"
    let mut j = 0; // Index into "worst"
    let mut wanted = HashSet::new(); // Target names of the instances we want
    let mut coverage_reached = false; // Flag for sufficient coverage has been met

    loop {
        // Exit if target number of instances (100?)
        // AND the coverage depth target ("coverageReached") has been met.
        // OR if we've exhausted the best/worse arrays
        if wanted.len() >= args.config.general.max_num_instances && coverage_reached
            || ((i == best.len()) && (j == worst.len()))
        {
            break;
        }

        let (aln, partition) = if i < best.len() {
            let val = best[i].clone();
            i += 1;
            (val, Partition::Top)
        } else {
            let val = worst[j].clone();
            j += 1;
            (val, Partition::Bottom)
        };

        let bins: Vec<_> =
            ((aln.subject_start / 10)..(aln.subject_end.div_ceil(10))).collect();
        let new_cov: Vec<_> = bins.iter().map(|i| consensus_cov[*i] + 1).collect();
        let supports_cov = new_cov.iter().any(|&val| val <= 10);

        if partition == Partition::Top || supports_cov {
            // Add the instance
            wanted.insert(aln.target);

            // Increment the consensus coverage
            for bin in bins {
                consensus_cov[bin] += 1;
            }
        }

        // Iterate over all positions in the coverage array.
        // If all positions in the consensus are above the coverage
        // depth target (10), set flag ("coverageReached")
        coverage_reached = consensus_cov
            .iter()
            .all(|val| *val >= args.config.general.min_consensus_coverage);
    }

    let mut num_taken = 0;
    if !wanted.is_empty() {
        let mut reader = FastaReader::new(BufReader::new(open(args.from_path)?));
        let mut fasta_writer =
            FastaWriter::new(BufWriter::new(open_for_write(args.to_path)?));
        for record in reader.records().map_while(Result::ok) {
            let name = String::from_utf8(record.name().to_vec())?;
            if wanted.contains(&name) {
                fasta_writer.write_record(&record)?;
                num_taken += 1;
            }
        }
    }

    Ok(num_taken)
}

// --------------------------------------------------
fn take_best_alignments(
    blast_out: &PathBuf,
    output: &PathBuf,
    config: &Config,
) -> Result<()> {
    let now = Instant::now();

    if fs::metadata(output).map_or(0, |meta| meta.len()) > 0 {
        debug!("Reusing best alignments '{}'", output.display());
    } else {
        debug!("Taking best alignments from '{}'", blast_out.display());
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(blast_out)
            .map_err(|e| anyhow!("Failed to read '{}': {e}", &blast_out.display()))?;

        let mut writer = csv::WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(output)
            .map_err(|e| anyhow!("Failed to write '{}': {e}", &output.display()))?;

        let mut taken = HashMap::<StringPair, RmBlastOutput>::new();
        for res in reader.records() {
            let record: RmBlastOutput = res?.deserialize(None)?;

            // Be sure to allow record.query == record.target
            // to identify singletons that only cluster to self
            if record.pident >= config.general.percent_id_for_components * 100.0 {
                taken
                    .entry(StringPair(record.query.clone(), record.target.clone()))
                    .and_modify(|val| {
                        if val.score > record.score {
                            *val = record.clone();
                        }
                    })
                    .or_insert(record);
            }
        }

        for record in taken.values() {
            writer.serialize(record)?;
        }

        debug!("Wrote to '{}' in {:?}", output.display(), now.elapsed());
    }

    Ok(())
}

// --------------------------------------------------
fn parse_alignment(blast_out: &PathBuf) -> Result<Vec<RmBlastOutput>> {
    // BLAST output fails to include headers
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(blast_out)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &blast_out.display()))?;

    let mut records = vec![];
    for res in reader.records() {
        let rec = res?;
        let blast_rec: RmBlastOutput = rec.deserialize(None)?;
        records.push(blast_rec);
    }

    Ok(records)
}
