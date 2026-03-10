mod graph;

use anyhow::{anyhow, bail, Result};
use bio::alphabets::dna::revcomp;
use chrono::Duration;
use clap::{builder::PossibleValue, Parser, ValueEnum};
use csv::{ReaderBuilder, WriterBuilder};
use itertools::Itertools;
use kseq::parse_reader;
use log::debug;
use newick::Newick;
use noodles_fasta::{
    self,
    io::Reader as FastaReader,
    io::Writer as FastaWriter,
    record::{Definition as FastaDefinition, Record as FastaRecord},
};
use rayon::prelude::*;
use regex::Regex;
use serde::{Deserialize, Serialize};
use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    fmt,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use which::which;

#[derive(Debug, Clone, PartialEq)]
pub enum SequenceAlphabet {
    Dna,
    Protein,
}

impl ValueEnum for SequenceAlphabet {
    fn value_variants<'a>() -> &'a [Self] {
        &[SequenceAlphabet::Dna, SequenceAlphabet::Protein]
    }

    fn to_possible_value<'a>(&self) -> Option<PossibleValue> {
        Some(match self {
            SequenceAlphabet::Dna => PossibleValue::new("dna"),
            SequenceAlphabet::Protein => PossibleValue::new("protein"),
        })
    }
}

#[derive(Parser, Debug)]
pub enum Command {
    /// Generate config TOML
    Config(ConfigArgs),

    /// Align consensi to self and create "components"
    Components(ComponentsArgs),

    /// Cluster output from "components"
    Cluster(ClusterArgs),

    /// Concatenate singletons from "components" and output from "cluster"
    Concat(ConcatArgs),

    /// Run all steps (components, cluster, concat)
    Run(RunArgs),
}

/// SCULU subfamily clustering tool
#[derive(Parser, Debug)]
#[command(arg_required_else_help = true, version, about)]
pub struct Cli {
    /// SCULU command
    #[command(subcommand)]
    pub command: Option<Command>,

    /// Log output
    #[arg(long, value_name = "LOGFILE")]
    pub logfile: Option<PathBuf>,

    /// Number of threads for rmblastn/Refiner
    #[arg(long, value_name = "THREADS")]
    pub num_threads: Option<usize>,
}

#[derive(Debug, Parser)]
#[command()]
pub struct ConfigArgs {
    /// Output file
    #[arg(value_name = "OUTFILE", default_value = "sculu.toml")]
    pub outfile: PathBuf,
}

#[derive(Debug, Parser)]
#[command()]
pub struct ConcatArgs {
    /// Output file
    #[arg(short, long, value_name = "OUTFILE", required = true)]
    pub outfile: PathBuf,

    /// The filtered FASTA consensi from "components"
    #[arg(long, value_name = "CONSENSI", required = true)]
    pub consensi: PathBuf,

    /// Singletons file from "components"
    #[arg(long, value_name = "SINGLETONS")]
    pub singletons: Option<PathBuf>,

    /// Merged components from "cluster"
    #[arg(long, value_name = "COMPONENTS", num_args = 0..)]
    pub components: Vec<PathBuf>,
}

#[derive(Debug, Parser, Clone)]
#[command()]
pub struct RunArgs {
    /// Sequence alphabet
    #[arg(short, long, value_name = "ALPHABET", required = true)]
    pub alphabet: SequenceAlphabet,

    /// FASTA file of subfamily consensi
    #[arg(long, value_name = "CONSENSI", required = true)]
    pub consensi: PathBuf,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES", required = true)]
    pub instances: PathBuf,

    /// One file from the "components" action
    #[arg(long, value_name = "COMPONENT")]
    pub component: Option<PathBuf>,

    /// Output directory
    #[arg(long, value_name = "OUTDIR")]
    pub outdir: PathBuf,

    /// Output file
    #[arg(long, value_name = "OUTFILE", default_value = "families.fa")]
    pub outfile: PathBuf,

    /// Config file
    #[arg(long, value_name = "CONFIG")]
    pub config: Option<PathBuf>,
}

#[derive(Debug, Parser, Clone)]
#[command()]
pub struct ComponentsArgs {
    /// Sequence alphabet
    #[arg(short, long, value_name = "ALPHABET", required = true)]
    pub alphabet: SequenceAlphabet,

    /// FASTA file of subfamily consensi
    #[arg(long, value_name = "CONSENSI", required = true)]
    pub consensi: PathBuf,

    /// Directory of instance files for each subfamily
    #[arg(long, value_name = "INSTANCES", required = true)]
    pub instances: PathBuf,

    /// Output directory
    #[arg(long, value_name = "OUTDIR", default_value = "sculu-out")]
    pub outdir: PathBuf,

    /// Config file
    #[arg(long, value_name = "CONFIG")]
    pub config: Option<PathBuf>,
}

#[derive(Debug, Parser, Clone)]
#[command()]
pub struct ClusterArgs {
    /// Sequence alphabet
    #[arg(short, long, value_name = "ALPHABET", required = true)]
    pub alphabet: SequenceAlphabet,

    /// The alignment.tsv file from "components"
    #[arg(long, value_name = "ALIGNMENTS", required = true)]
    pub alignments: PathBuf,

    /// The filtered FASTA consensi from "components"
    #[arg(long, value_name = "CONSENSI", required = true)]
    pub consensi: PathBuf,

    /// Directory of filtered instance files from "components"
    #[arg(long, value_name = "INSTANCES", required = true)]
    pub instances: PathBuf,

    /// A file from "components" action
    #[arg(long, value_name = "COMPONENT", required = true)]
    pub component: PathBuf,

    /// Output directory
    #[arg(long, value_name = "OUTDIR")]
    pub outdir: PathBuf,

    /// Config file
    #[arg(long, value_name = "CONFIG")]
    pub config: Option<PathBuf>,
}

/// Struct of sculu.toml file
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Config {
    general: GeneralConfig,
    rmblastn: RmblastnConfig,
    blastp: BlastpConfig,
}

/// "general" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
struct GeneralConfig {
    confidence_margin: isize,
    independence_threshold: f64,
    lambda: f64,
    percent_id_for_components: f64,
    max_num_instances: usize,
    min_align_cover: f64,
    min_consensus_coverage: usize,
    min_instance_sequence_length_dna: usize,
    min_instance_sequence_length_prot: usize,
    min_len_similarity: f64,
    min_num_instances_dna: usize,
    min_num_instances_prot: usize,
}

/// "blastp" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
struct BlastpConfig {
    matrix: Option<PathBuf>,
    gap_open: usize,
    gap_extend: usize,
    word_size: usize,
    mask_level: usize,
    dust: bool,
    complexity_adjust: bool,
    min_raw_gapped_score: usize,
    xdrop_gap: usize,
    xdrop_ungap: usize,
    xdrop_gap_final: usize,
}

/// "rmblast" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
struct RmblastnConfig {
    matrix: Option<PathBuf>,
    gap_open: usize,
    gap_extend: usize,
    word_size: usize,
    mask_level: usize,
    dust: bool,
    complexity_adjust: bool,
    min_raw_gapped_score: usize,
    xdrop_gap: usize,
    xdrop_ungap: usize,
    xdrop_gap_final: usize,
}

#[derive(Debug, PartialEq)]
enum Partition {
    Top,
    Bottom,
}

#[derive(Debug, Deserialize, Serialize, Clone, PartialEq)]
struct RmBlastOutput {
    score: usize,
    target: String,
    query: String,
    query_len: usize,
    query_start: usize,
    query_end: usize,
    subject_len: usize,
    subject_start: usize,
    subject_end: usize,
    cpg_kdiv: f64,
    pident: f64,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
struct AlignedConsensusPair {
    target: String,
    query: String,
    is_flipped: bool,
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
enum Direction {
    #[serde(rename = "forward")]
    Forward,
    #[serde(rename = "reverse")]
    Reverse,
}

#[derive(Debug, Serialize, Deserialize)]
struct AlignmentScore {
    score: usize,
    target: String,
    query: String,
}

#[derive(Debug, PartialEq)]
struct Independence {
    f1: String,
    f2: String,
    val: f64,
}

#[derive(Debug)]
pub struct Components {
    pub singletons: Option<PathBuf>,
    pub alignments: PathBuf,
    pub components: Vec<PathBuf>,
}

#[derive(Debug)]
pub struct BuiltComponents {
    pub singletons: Option<PathBuf>,
    pub alignments: PathBuf,
    pub components: Vec<PathBuf>,
    pub consensi: PathBuf,
    pub instances_dir: PathBuf,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
struct StringPair(String, String);

impl StringPair {
    pub fn new(s1: String, s2: String) -> StringPair {
        // Store the strings in ascending order
        if s1 < s2 {
            StringPair(s1, s2)
        } else {
            StringPair(s2, s1)
        }
    }
}

impl fmt::Display for StringPair {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}/{}", self.0, self.1)
    }
}

#[derive(Debug)]
struct Winners {
    clear_winners: HashMap<String, u32>,
    winning_sets: HashMap<StringPair, u32>,
}

#[derive(Debug)]
struct MergeFamilies<'a> {
    family1: String,
    family2: String,
    outdir: PathBuf,
    taken_instances_dir: &'a PathBuf,
    num_threads: usize,
    alphabet: SequenceAlphabet,
    flipped: bool,
}

// --------------------------------------------------
fn default_config() -> Config {
    Config {
        general: GeneralConfig {
            confidence_margin: 3,
            independence_threshold: 0.5,
            lambda: 0.1227,
            percent_id_for_components: 0.70,
            max_num_instances: 100,
            min_align_cover: 0.9,
            min_consensus_coverage: 5,
            min_instance_sequence_length_dna: 30,
            min_instance_sequence_length_prot: 12,
            min_len_similarity: 0.9,
            min_num_instances_dna: 10,
            min_num_instances_prot: 1,
        },
        blastp: BlastpConfig {
            matrix: Some(PathBuf::from("/path/to/matrix.txt")),
            gap_open: 20,
            gap_extend: 5,
            word_size: 7,
            mask_level: 101,
            min_raw_gapped_score: 400,
            xdrop_gap: 800,
            xdrop_ungap: 200,
            xdrop_gap_final: 400,
            dust: false,
            complexity_adjust: false,
        },
        rmblastn: RmblastnConfig {
            matrix: Some(PathBuf::from("/path/to/matrix.txt")),
            gap_open: 20,
            gap_extend: 5,
            word_size: 7,
            mask_level: 101,
            min_raw_gapped_score: 400,
            xdrop_gap: 800,
            xdrop_ungap: 200,
            xdrop_gap_final: 400,
            dust: false,
            complexity_adjust: false,
        },
    }
}

// --------------------------------------------------
pub fn write_config(args: &ConfigArgs) -> Result<()> {
    if args.outfile.exists() {
        bail!(r#"Will not overwrite "{}""#, args.outfile.display());
    }

    let mut outfile = open_for_write(&args.outfile)?;
    writeln!(&mut outfile, "{}", toml::to_string(&default_config())?)?;
    println!(r#"Wrote "{}""#, args.outfile.display());
    Ok(())
}

// --------------------------------------------------
pub fn build_components(
    args: &ComponentsArgs,
    num_threads: usize,
) -> Result<BuiltComponents> {
    let config = match args.config.as_ref() {
        Some(file) => get_config(file)?,
        _ => default_config(),
    };

    if !&args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    let taken_instances_dir = args.outdir.join("instances");
    fs::create_dir_all(&taken_instances_dir)?;

    let (consensi_file, _family_to_instance) = check_family_instances(
        &args.consensi,
        &args.instances,
        &args.outdir,
        &taken_instances_dir,
        &config,
        &args.alphabet,
        num_threads,
    )?;

    // This will only BLAST if there are no components from a previous run
    let components = align_consensi_to_self(
        &consensi_file,
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
        consensi: consensi_file,
        instances_dir: taken_instances_dir,
    })
}

// --------------------------------------------------
pub fn process_component(args: &ClusterArgs, num_threads: usize) -> Result<PathBuf> {
    let component_file = &args.component;
    let config = match args.config.as_ref() {
        Some(file) => get_config(file)?,
        _ => default_config(),
    };

    if !&args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    // If given an component file to process, assume the
    // consensi and instances have been properly filtered
    // and use given args as-is.
    let family_to_instance = read_instances_dir(&args.instances)?;

    let merged = merge_component(
        &args.outdir,
        &args.alphabet,
        component_file,
        &args.consensi,
        family_to_instance,
        &args.instances,
        &args.alignments,
        &config,
        num_threads,
    )?;

    debug!(
        "'{}' merged into '{}'",
        component_file.display(),
        merged.display()
    );

    Ok(merged)
}

// --------------------------------------------------
pub fn concat_files(
    consensi: &PathBuf,
    singletons: &Option<PathBuf>,
    merged: &Vec<PathBuf>,
    outfile: &PathBuf,
) -> Result<()> {
    let mut fasta_writer = FastaWriter::new(BufWriter::new(open_for_write(outfile)?));

    if let Some(file) = singletons {
        let singletons = read_lines(file)?;
        debug!(
            "Copying {} sequences from singletons file",
            singletons.len()
        );
        copy_fasta(&singletons, consensi, &mut fasta_writer)?;
    }

    for component_file in merged {
        debug!(r#"Copying from "{}""#, component_file.display());
        let mut reader = FastaReader::new(BufReader::new(open(component_file)?));
        for record in reader.records().map_while(Result::ok) {
            fasta_writer.write_record(&record)?;
        }
    }

    debug!(r#"Final output written to "{}""#, outfile.display());

    Ok(())
}

// --------------------------------------------------
/// Read the (optional) "sculu.toml" file
fn get_config(config_file: &PathBuf) -> Result<Config> {
    let mut file = open(config_file)?;
    let mut contents = String::new();
    let _bytes = file.read_to_string(&mut contents);
    let config: Config = toml::from_str(&contents)?;
    Ok(config)
}

// --------------------------------------------------
// Align the consensi to itself to identify clusters, e.g.:
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
pub fn align_consensi_to_self(
    consensi: &Path,
    out_dir: &PathBuf,
    alphabet: &SequenceAlphabet,
    config: &Config,
    num_threads: usize,
) -> Result<Components> {
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
        let blast_dir = out_dir.join("consensi_cluster");
        let blast_out = match alphabet {
            SequenceAlphabet::Dna => {
                run_rmblastn(&blast_dir, consensi, consensi, config, num_threads)?
            }
            _ => run_blastp(&blast_dir, consensi, consensi, config, num_threads)?,
        };

        // There may be multiple hits per pair, take highest
        let best_alignments = blast_dir.join("best.tsv");
        take_best_alignments(&blast_out, &best_alignments, config)?;

        // Filter out: Hits to self, insufficient coverage
        let alignment_file = components_dir.join("alignment.tsv");
        let mut alignment_wtr = WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(&alignment_file)
            .map_err(|e| {
                anyhow!("Failed to write '{}': {e}", &alignment_file.display())
            })?;

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(&best_alignments)
            .map_err(|e| {
                anyhow!("Failed to read '{}': {e}", &best_alignments.display())
            })?;

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
            alignment_wtr.serialize(AlignedConsensusPair {
                target: record.target.to_string(),
                query: record.query.clone(),
                is_flipped: query_dir != subject_dir,
            })?;
            records.push(record);
        }

        let mut components = graph::connected_components(records);
        // Sort by size of component group ascending
        components.sort_by_key(|v| v.len());

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
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_path(blast_out)
            .map_err(|e| anyhow!("Failed to read '{}': {e}", &blast_out.display()))?;

        let mut writer = WriterBuilder::new()
            .has_headers(true)
            .delimiter(b'\t')
            .from_path(output)
            .map_err(|e| anyhow!("Failed to write '{}': {e}", &output.display()))?;

        let mut taken = HashMap::<StringPair, RmBlastOutput>::new();
        for res in reader.records() {
            let record: RmBlastOutput = res?.deserialize(None)?;

            if record.query == record.target {
                continue;
            }

            // The coverage of the hits to either the query or subject must
            // be at least half of the sequence length.
            //let query_span = 1 + record.query_end.abs_diff(record.query_start);
            //let subject_span = 1 + record.subject_end.abs_diff(record.subject_start);
            //let query_coverage = query_span as f64 / record.query_len as f64;
            //let subject_coverage = subject_span as f64 / record.subject_len as f64;
            //let are_equiv_len = min(query_span, subject_span) as f64
            //    >= (config.general.min_len_similarity
            //        * max(query_span, subject_span) as f64);
            //if query_coverage >= config.general.min_align_cover
            //    && subject_coverage >= config.general.min_align_cover
            //    && are_equiv_len
            //{
            debug!(
                "target = {}, query = {}, pident = {}",
                record.target, record.query, record.pident
            );

            if record.pident >= config.general.percent_id_for_components * 100.0 {
                let pair = StringPair(record.query.clone(), record.target.clone());
                debug!("ADDING {:?}", pair);

                if let Some(val) = taken.get_mut(&pair) {
                    if val.score > record.score {
                        *val = record.clone();
                    }
                } else {
                    taken.insert(pair, record.clone());
                }
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
fn copy_fasta<W: Write>(
    wanted_families: &[String],
    source: &PathBuf,
    destination: &mut FastaWriter<W>,
) -> Result<usize> {
    let mut reader = FastaReader::new(open(source)?);
    let mut num_taken = 0;
    for result in reader.records() {
        let record = result?;
        let family = String::from_utf8(record.name().to_vec())?;
        if wanted_families.contains(&family) {
            destination.write_record(&record)?;
            num_taken += 1;
        }
    }

    Ok(num_taken)
}

// --------------------------------------------------
fn merge_component(
    out_dir: &PathBuf,
    alphabet: &SequenceAlphabet,
    component_file: &PathBuf,
    consensi_file: &Path,
    mut family_to_instance: HashMap<String, PathBuf>,
    taken_instances_dir: &PathBuf,
    alignments: &PathBuf,
    config: &Config,
    num_threads: usize,
) -> Result<PathBuf> {
    debug!("Merging component '{}'", component_file.display());

    // The "component-N" file will contain the names of the families
    let families = read_lines(&component_file.to_path_buf())?;

    // Get the alignments underlying this component
    let mut alignment_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(alignments)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &alignments.display()))?;

    let mut flipped_pairs: HashSet<StringPair> = HashSet::new();
    for res in alignment_reader.records() {
        let record: AlignedConsensusPair = res?.deserialize(None)?;
        if record.is_flipped
            && (families.contains(&record.target) || families.contains(&record.query))
        {
            let pair = StringPair::new(record.target.clone(), record.query.clone());
            if !flipped_pairs.contains(&pair) {
                flipped_pairs.insert(pair);
            }
        }
    }

    // Create a working dir with the same name as the component file
    let component_name = component_file
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let batch_dir = out_dir.join(&component_name);
    fs::create_dir_all(&batch_dir)?;

    debug!(
        "==> Component '{component_name}' has {} families <==",
        families.len()
    );

    let outfile = batch_dir.join("final.fa");
    let mut final_writer = FastaWriter::new(BufWriter::new(open_for_write(&outfile)?));
    let start = Instant::now();

    // Create a consensi file for this round containing only the given families
    let mut batch_consensi = batch_dir.join("consensi.fa");
    {
        // Scoped to cause fasta_writer to close
        let mut fasta_writer =
            FastaWriter::new(BufWriter::new(open_for_write(&batch_consensi)?));
        let num_taken =
            copy_fasta(&families, &consensi_file.to_path_buf(), &mut fasta_writer)?;

        if num_taken == 0 {
            bail!(
                "Failed to copy consensi sequences from '{}' to '{}'",
                consensi_file.display(),
                batch_consensi.display()
            );
        }
    }

    // Create query file for families
    let all_seqs_path = &batch_dir.join("all_seqs.fa");
    debug!(
        "Concatenating all instance sequences into '{}'",
        all_seqs_path.display()
    );

    // Put the instances for the given families into a single file
    cat_sequences(taken_instances_dir, &families, all_seqs_path)?;

    // Create a starting directory for the merge iterations.
    let round0_dir = batch_dir.join("round000");
    fs::create_dir_all(&round0_dir)?;

    // As the consensi are merged, the family names are concatenated
    // into Newick-formatted strings to track the merges.
    // These can get quite long, causing `makeblastdb` to fail.
    // Make the sequence IDs simple integers by numbering
    // and move the family names to the description.
    let consensi_path = round0_dir.join("consensi.fa");
    debug!("Writing numbered consensi to '{}'", consensi_path.display());
    let mut consensi_seqs = number_consensi(&batch_consensi, &consensi_path)?;
    batch_consensi = consensi_path;

    // Start merging
    let trailing_semi = Regex::new(";$").unwrap();
    let mut prev_scores: Option<PathBuf> = None;
    for round in 1.. {
        debug!(">>> Round {round} <<<");
        let round_dir = batch_dir.join(format!("round{round:03}"));
        fs::create_dir_all(&round_dir)?;

        // Align all the sequences to the current consensi.
        // On the first round, all the original consensi will be included.
        // On future rounds, only the newly merged consensi will be present.
        let alignment_file = match alphabet {
            SequenceAlphabet::Dna => run_rmblastn(
                &round_dir,
                &batch_consensi,
                all_seqs_path,
                config,
                num_threads,
            )?,
            _ => run_blastp(
                &round_dir,
                &batch_consensi,
                all_seqs_path,
                config,
                num_threads,
            )?,
        };

        // Extract the scores from the alignment file.
        // On the first round, there will be no "prev_scores" file.
        // On the following rounds, the previous round's scores will be
        // included for those families that were not merged.
        let scores_file =
            extract_scores(&alignment_file, &prev_scores, &batch_consensi, &round_dir)?;

        // Save this round's alignment scores for the next
        prev_scores = Some(scores_file.clone());

        // Find the clear/ambiguous winners from the scores
        let winners = call_winners(
            &scores_file,
            config.general.lambda,
            config.general.confidence_margin,
        )?;
        debug!("winners = {:#?}", winners);

        // Determine independence of all pairs
        let pair_independence = independence(winners);

        // Create a lookup for the scores
        let mut score_lookup: HashMap<StringPair, f64> = HashMap::new();
        for pair in &pair_independence {
            score_lookup.insert(
                StringPair(pair.f1.to_string(), pair.f2.to_string()),
                pair.val,
            );
        }

        // Select only those pair lacking independence
        let non_independent: Vec<_> = pair_independence
            .iter()
            .filter(|v| v.val < config.general.independence_threshold)
            .collect();

        debug!("{} family pair not independent.", non_independent.len());
        debug!("{non_independent:#?}");

        // Stop the loop when all families are independent
        if non_independent.is_empty() {
            break;
        }

        // Merge the least independent families.
        // The new consensi file will only contain the merged families.
        let new_consensi_path = &round_dir.join("new-consensi.fa");
        let mut new_consensi = open_for_write(new_consensi_path)?;
        let mut merge_num = 0;
        let mut already_merged: HashSet<String> = HashSet::new();
        for pair in non_independent {
            // The family names may be in Newick format
            let fams1 = parse_newick(&pair.f1);
            let fams2 = parse_newick(&pair.f2);

            // Collect all the family names involved in this merge
            let mut all_families: Vec<_> = fams1.clone();
            for f2 in &fams2 {
                all_families.push(f2.clone());
            }

            // We might see Fam1->Fam2 and later Fam2->Fam1.
            // Or we might later see Fam1->Fam3.
            // Only merge a family once.
            if all_families.iter().any(|f| already_merged.contains(f)) {
                debug!("Already merged one of {}", all_families.join(", "));
                continue;
            }

            // Increment the number of merged pairs
            merge_num += 1;

            // See if the pair has a score from the other direction.
            // The merges should happen in order from least independent
            // to greater, so this A/B merge "val" will be lower than the
            // symmetrical B/A score.
            let other_score = score_lookup
                .get(&StringPair::new(pair.f2.to_string(), pair.f1.to_string()))
                .map_or("".to_string(), |v| format!("/{v:0.04}"));

            debug!(
                "{}: Merge {} => {} Ind: {:0.04}{other_score}",
                merge_num, &pair.f1, &pair.f2, pair.val
            );

            // Place all merge artefacts into a directory.
            let msa_dir = &round_dir.join(format!("msa-{merge_num:02}"));
            let new_consensus_seq = merge_families(
                MergeFamilies {
                    family1: pair.f1.clone(),
                    family2: pair.f2.clone(),
                    outdir: msa_dir.to_path_buf(),
                    taken_instances_dir,
                    num_threads,
                    alphabet: alphabet.clone(),
                    flipped: flipped_pairs
                        .contains(&StringPair::new(pair.f1.clone(), pair.f2.clone())),
                },
                &mut family_to_instance,
            )?;

            let new_family_newick = format!(
                "({},{}):{:0.02};",
                trailing_semi.replace(&pair.f1, ""),
                trailing_semi.replace(&pair.f2, ""),
                pair.val
            );

            let f1_len = consensi_seqs.get(&pair.f1).map_or(0, |v| v.len());
            let f2_len = consensi_seqs.get(&pair.f2).map_or(0, |v| v.len());
            let new_len = new_consensus_seq.len();

            debug!(
                "{} len was {f1_len}, {} len was {f2_len}, \
                    new consensi len is {new_len}",
                pair.f1, pair.f2
            );

            // Update the consensi mapping
            consensi_seqs.remove(&pair.f1);
            consensi_seqs.remove(&pair.f2);
            consensi_seqs.insert(new_family_newick, new_consensus_seq.clone());

            // Keep track of all the families that have been merged
            for family in all_families {
                already_merged.insert(family);
            }
        }

        // Write the merged families to the new consensi file
        debug!(
            "Writing {} consensi to '{}'",
            &consensi_seqs.len(),
            new_consensi_path.display()
        );

        for (family_number, (family_name, seq)) in consensi_seqs.iter().enumerate() {
            writeln!(new_consensi, ">{family_number} {family_name}\n{seq}",)?;
        }

        batch_consensi = new_consensi_path.to_path_buf();
    }

    let mut reader = FastaReader::new(open(&batch_consensi)?);
    let mut new_seqs = 0;
    for result in reader.records() {
        let mut record = result?;
        if let Some(desc) = record.description() {
            record = FastaRecord::new(
                FastaDefinition::new(desc, None),
                record.sequence().clone(),
            )
        }
        final_writer.write_record(&record)?;
        new_seqs += 1;
    }

    debug!(
        "Added {new_seqs} famil{} in {}.",
        if new_seqs == 1 { "y" } else { "ies" },
        format_seconds(start.elapsed().as_secs()),
    );

    Ok(outfile)
}

// --------------------------------------------------
fn parse_alignment(blast_out: &PathBuf) -> Result<Vec<RmBlastOutput>> {
    // BLAST output fails to include headers
    let mut reader = ReaderBuilder::new()
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

// --------------------------------------------------
fn run_blastp(
    out_dir: &PathBuf,
    db: &Path,
    query: &Path,
    config: &Config,
    num_threads: usize,
) -> Result<PathBuf> {
    fs::create_dir_all(out_dir)?;
    let outfile = out_dir.join("blast.tsv");

    if outfile.exists() {
        debug!("Reusing BLAST output file '{}'", outfile.display());
    } else {
        let db_path = &out_dir.join("db");
        let makeblastdb_args = &[
            "-out".to_string(),
            db_path.to_string_lossy().to_string(),
            "-parse_seqids".to_string(),
            "-dbtype".to_string(),
            "prot".to_string(),
            "-in".to_string(),
            db.to_string_lossy().to_string(),
        ];

        debug!("Running 'makeblastdb {}'", &makeblastdb_args.join(" "));
        let makeblastdb =
            which("makeblastdb").map_err(|e| anyhow!("makeblastdb: {e}"))?;
        let res = std::process::Command::new(makeblastdb)
            .args(makeblastdb_args)
            .output()?;

        if !res.status.success() {
            bail!(String::from_utf8(res.stderr)?);
        }

        let blastp = which("blastp").map_err(|e| anyhow!("blastp: {e}"))?;
        let mut cmd = std::process::Command::new(&blastp);
        let mut blastp_args = vec![
            "-db".to_string(),
            db_path.to_string_lossy().to_string(),
            "-query".to_string(),
            query.to_string_lossy().to_string(),
            "-out".to_string(),
            outfile.to_string_lossy().to_string(),
            "-outfmt".to_string(),
            "6 score qseqid sseqid qlen qstart qend slen sstart send cpg_kdiv pident"
                .to_string(),
            "-num_threads".to_string(),
            num_threads.to_string(),
            "-num_alignments".to_string(),
            "9999999".to_string(),
            "-gapopen".to_string(),
            config.blastp.gap_open.to_string(),
            "-gapextend".to_string(),
            config.blastp.gap_extend.to_string(),
            "-word_size".to_string(),
            config.blastp.word_size.to_string(),
            "-xdrop_ungap".to_string(),
            (config.blastp.min_raw_gapped_score * 2).to_string(),
            "-xdrop_gap".to_string(),
            (config.blastp.min_raw_gapped_score / 2).to_string(),
            "-xdrop_gap_final".to_string(),
            config.blastp.min_raw_gapped_score.to_string(),
            "-mask_level".to_string(),
            config.blastp.mask_level.to_string(),
            "-min_raw_gapped_score".to_string(),
            config.blastp.min_raw_gapped_score.to_string(),
            "-dust".to_string(),
            if config.blastp.dust {
                "yes".to_string()
            } else {
                "no".to_string()
            },
        ];

        if config.blastp.complexity_adjust {
            blastp_args.push("-complexity_adjust".to_string());
        }

        if let Some(matrix) = config.blastp.matrix.as_ref() {
            if !matrix.is_file() {
                bail!("blastp matrix '{}' does not exist", matrix.display());
            }

            let matrix_filename = matrix.file_name().unwrap_or_else(|| {
                panic!("Failed to get filename from '{}'", matrix.display())
            });

            let matrix_dir = matrix.parent().unwrap_or_else(|| {
                panic!("Failed to get dirname from '{}'", matrix.display())
            });

            cmd.env("BLASTMAT", matrix_dir);
            blastp_args.extend_from_slice(&[
                "-matrix".to_string(),
                matrix_filename.to_string_lossy().to_string(),
            ]);
        }

        debug!("Running '{} {}'", blastp.display(), blastp_args.join(" "));

        let start = Instant::now();
        let res = cmd.args(blastp_args).output()?;
        if !res.status.success() {
            bail!(String::from_utf8(res.stderr)?);
        }

        debug!(
            "blastp finished in {}",
            format_seconds(start.elapsed().as_secs())
        );
    }

    Ok(outfile)
}

// --------------------------------------------------
fn run_rmblastn(
    out_dir: &PathBuf,
    db: &Path,
    query: &Path,
    config: &Config,
    num_threads: usize,
) -> Result<PathBuf> {
    fs::create_dir_all(out_dir)?;
    let outfile = out_dir.join("blast.tsv");

    if outfile.exists() {
        debug!("Reusing BLAST output file '{}'", outfile.display());
    } else {
        let db_path = &out_dir.join("db");
        let makeblastdb_args = &[
            "-out".to_string(),
            db_path.to_string_lossy().to_string(),
            "-parse_seqids".to_string(),
            "-dbtype".to_string(),
            "nucl".to_string(),
            "-in".to_string(),
            db.to_string_lossy().to_string(),
        ];

        debug!("Running 'makeblastdb {}'", &makeblastdb_args.join(" "));
        let makeblastdb =
            which("makeblastdb").map_err(|e| anyhow!("makeblastdb: {e}"))?;
        let res = std::process::Command::new(makeblastdb)
            .args(makeblastdb_args)
            .output()?;

        if !res.status.success() {
            bail!(String::from_utf8(res.stderr)?);
        }

        let rmblastn = which("rmblastn").map_err(|e| anyhow!("rmblastn: {e}"))?;
        let mut cmd = std::process::Command::new(&rmblastn);
        let mut rmblastn_args = vec![
            "-db".to_string(),
            db_path.to_string_lossy().to_string(),
            "-query".to_string(),
            query.to_string_lossy().to_string(),
            "-out".to_string(),
            outfile.to_string_lossy().to_string(),
            "-outfmt".to_string(),
            "6 score qseqid sseqid qlen qstart qend slen sstart send cpg_kdiv pident"
                .to_string(),
            "-num_threads".to_string(),
            num_threads.to_string(),
            "-num_alignments".to_string(),
            "9999999".to_string(),
            "-mask_level".to_string(),
            config.rmblastn.mask_level.to_string(),
            "-gapopen".to_string(),
            config.rmblastn.gap_open.to_string(),
            "-gapextend".to_string(),
            config.rmblastn.gap_extend.to_string(),
            "-word_size".to_string(),
            config.rmblastn.word_size.to_string(),
            "-min_raw_gapped_score".to_string(),
            config.rmblastn.min_raw_gapped_score.to_string(),
            "-xdrop_ungap".to_string(),
            (config.rmblastn.min_raw_gapped_score * 2).to_string(),
            "-xdrop_gap".to_string(),
            (config.rmblastn.min_raw_gapped_score / 2).to_string(),
            "-xdrop_gap_final".to_string(),
            config.rmblastn.min_raw_gapped_score.to_string(),
            "-dust".to_string(),
            if config.rmblastn.dust {
                "yes".to_string()
            } else {
                "no".to_string()
            },
        ];

        if config.rmblastn.complexity_adjust {
            rmblastn_args.push("-complexity_adjust".to_string());
        }

        if let Some(matrix) = config.rmblastn.matrix.as_ref() {
            if !matrix.is_file() {
                bail!("Rmblast matrix '{}' does not exist", matrix.display());
            }

            let matrix_filename = matrix.file_name().unwrap_or_else(|| {
                panic!("Failed to get filename from '{}'", matrix.display())
            });

            let matrix_dir = matrix.parent().unwrap_or_else(|| {
                panic!("Failed to get dirname from '{}'", matrix.display())
            });

            cmd.env("BLASTMAT", matrix_dir);
            rmblastn_args.extend_from_slice(&[
                "-matrix".to_string(),
                matrix_filename.to_string_lossy().to_string(),
            ]);
        }

        debug!(
            "Running '{} {}'",
            rmblastn.display(),
            rmblastn_args.join(" ")
        );

        let start = Instant::now();
        let res = cmd.args(rmblastn_args).output()?;
        if !res.status.success() {
            bail!(String::from_utf8(res.stderr)?);
        }

        debug!(
            "Rmblastn finished in {}",
            format_seconds(start.elapsed().as_secs())
        );
    }

    Ok(outfile)
}

// --------------------------------------------------
fn bitscore_to_confidence(vals: &[&usize], lambda: f64) -> Result<Vec<f64>> {
    let mut converted: Vec<_> =
        vals.iter().map(|&&v| (v as f64 * lambda).exp2()).collect();

    if converted.iter().any(|v| v.is_infinite()) {
        // Scale all numbers down
        let delta = **vals.iter().max().unwrap() as i64 - 500;
        converted = vals
            .iter()
            .map(|&&v| ((v as i64 - delta) as f64 * lambda).exp2())
            .collect::<Vec<_>>();
    }

    let total: f64 = converted.iter().sum();
    if total > 0. {
        Ok(converted.into_iter().map(|v| v / total).collect())
    } else {
        bail!("Sum of converted values equals zero")
    }
}

// --------------------------------------------------
fn call_winners(
    scores_file: &PathBuf,
    lambda: f64,
    confidence_margin: isize,
) -> Result<Winners> {
    debug!(
        r#"Calling winners from scores file "{}""#,
        scores_file.display()
    );

    let mut reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(scores_file)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &scores_file.display()))?;
    let records = reader.deserialize();

    let mut scores: HashMap<String, HashMap<String, usize>> = HashMap::new();

    for res in records {
        let rec: AlignmentScore = res?;
        scores
            .entry(rec.target)
            .or_default()
            .entry(rec.query)
            .and_modify(|v| *v = max(rec.score, *v))
            .or_insert(rec.score);
    }

    let mut clear_winners: HashMap<String, u32> = HashMap::new();
    let mut winning_sets: HashMap<StringPair, u32> = HashMap::new();
    for (_target, bit_scores) in scores.iter() {
        let families: Vec<_> = bit_scores.keys().collect();
        let raw_scores: Vec<_> = bit_scores.values().collect();
        let conf: Vec<_> = bitscore_to_confidence(&raw_scores, lambda)?;
        let pos: Vec<_> = (0..conf.len()).collect();
        let mut all_comps: Vec<bool> = vec![];
        for &i in &pos {
            let val = conf[i];
            let others: Vec<_> =
                pos.iter().filter(|&&j| j != i).map(|&j| conf[j]).collect();
            let pairs: Vec<_> =
                std::iter::repeat_n(val, others.len()).zip(others).collect();
            let comps: Vec<_> = pairs
                .iter()
                .map(|(x, y)| (x * (1. / confidence_margin as f64)) > *y)
                .collect();
            all_comps.push(comps.iter().all(|&v| v));
        }

        let winners: Vec<String> = families
            .iter()
            .zip(all_comps)
            .filter(|&(_, win)| win)
            .map(|(fam, _)| fam.to_string())
            .collect();

        if winners.len() == 1 {
            let winner = winners.first().unwrap();
            clear_winners
                .entry(winner.to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
        } else {
            let mut fam_comps: Vec<_> = conf.iter().zip(families).collect();
            fam_comps.sort_by(|a, b| {
                b.0.partial_cmp(a.0).unwrap().then_with(|| a.1.cmp(b.1))
            });
            let (&top_conf, _) = fam_comps.first().unwrap();
            let threshold = top_conf * (1. / confidence_margin as f64);
            let winning_set: Vec<_> = fam_comps
                .iter()
                .filter_map(|(&confidence, fam)| {
                    (confidence >= threshold).then_some(fam)
                })
                .collect();

            // The permutations will include A/B and B/A
            // It's important to store the symmetrical keys
            // Even though this is a duplication of the data
            // So don't use StringPair::new, Ken.
            for pair in winning_set.into_iter().permutations(2) {
                if let [&f1, &f2] = pair[..] {
                    let key = StringPair(f1.to_string(), f2.to_string());
                    winning_sets.entry(key).and_modify(|v| *v += 1).or_insert(1);
                }
            }
        }
    }

    Ok(Winners {
        clear_winners,
        winning_sets,
    })
}

// --------------------------------------------------
fn cat_sequences(
    instances_dir: &PathBuf,
    families: &[String],
    outpath: &PathBuf,
) -> Result<()> {
    let mut output = open_for_write(outpath)?;
    let entries = fs::read_dir(instances_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, instances_dir.display()))?;

    for entry in entries {
        let file = entry?;
        let family_name = file
            .path()
            .file_stem()
            .expect("file_stem")
            .to_string_lossy()
            .to_string();

        if !families.contains(&family_name) {
            continue;
        }

        let mut reader = parse_reader(open(&file.path())?)?;
        while let Some(rec) = reader.iter_record()? {
            writeln!(
                output,
                ">{}__{}{}\n{}",
                family_name,
                rec.head(),
                rec.des(),
                rec.seq()
            )?;
        }
    }

    Ok(())
}

// --------------------------------------------------
// This function will first filter the instances for those that
// are most representative of the family, and then it will filter
// the consensi for those that have sufficient supporting instances.
fn check_family_instances(
    consensi: &PathBuf,
    instances_dir: &PathBuf,
    out_dir: &PathBuf,
    taken_instances_dir: &Path,
    config: &Config,
    alphabet: &SequenceAlphabet,
    num_threads: usize,
) -> Result<(PathBuf, HashMap<String, PathBuf>)> {
    debug!("Checking consensi_file '{}'", consensi.display());

    // Find all the input instance files
    let instances: Vec<_> = fs::read_dir(instances_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, instances_dir.display()))?
        .map_while(Result::ok)
        .collect();

    debug!(
        "Found {} instance files in '{}'",
        instances.len(),
        instances_dir.display(),
    );

    let working_dir = out_dir.join("select");
    fs::create_dir_all(&working_dir)?;

    let now = Instant::now();
    instances.par_iter().try_for_each(|entry| -> Result<()> {
        let instance_path = entry.path();
        if let Some(instance_stem) = instance_path.file_stem() {
            // Skip hidden files
            let family_name = instance_stem.to_string_lossy().to_string();
            if !family_name.starts_with(".") {
                let taken_path = taken_instances_dir.join(format!("{family_name}.fa"));
                if let Err(e) = select_instances(
                    consensi,
                    &family_name,
                    &instance_path,
                    &taken_path,
                    &working_dir,
                    alphabet,
                    config,
                    num_threads,
                ) {
                    eprintln!("Warning: {e}");
                }
            }
        } else {
            eprintln!(
                "Error: Cannot get filename for instance {}",
                instance_path.display()
            );
        }
        Ok(())
    })?;

    let family_to_instance = read_instances_dir(taken_instances_dir)?;
    debug!(
        "family_to_instance has {} members",
        family_to_instance.len()
    );
    debug!("Finished instance selection in {:?}", now.elapsed());

    // Create a new consensi file containing only the families that have instances
    let taken_consensi_path = out_dir.join("consensi.fa");
    let mut out_consensi = open_for_write(&taken_consensi_path)?;
    let mut reader = parse_reader(open(consensi)?)?;
    let mut consensi_names: HashMap<String, u32> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        let family = rec.head().to_string();
        if family_to_instance.contains_key(&family) {
            consensi_names
                .entry(rec.head().to_string())
                .and_modify(|v| *v += 1)
                .or_insert(1);
            writeln!(out_consensi, ">{family}\n{}", rec.seq())?;
        } else {
            debug!("Family '{family}' has no instances, removing!");
        }
    }

    let mut dups: Vec<_> = consensi_names
        .iter()
        .flat_map(|(name, &count)| (count > 1).then_some(name))
        .collect();

    if !dups.is_empty() {
        // Have to sort for tests
        dups.sort();
        bail!(
            "The following consensi IDs are duplicated: {}",
            dups.iter().join(", ")
        );
    }

    Ok((taken_consensi_path, family_to_instance))
}

// --------------------------------------------------
fn read_instances_dir(instances_dir: &Path) -> Result<HashMap<String, PathBuf>> {
    // Create hashmap from family name to taken instance file
    let mut family_to_instance: HashMap<String, PathBuf> = HashMap::new();
    let entries = fs::read_dir(instances_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, instances_dir.display()))?;

    for entry in entries {
        let entry = entry?;
        if let Some(stem) = entry.path().file_stem() {
            let family_name = stem.to_string_lossy().to_string();
            if !family_name.starts_with("inst-") {
                family_to_instance.insert(family_name, entry.path().to_path_buf());
            }
        }
    }
    Ok(family_to_instance)
}

// --------------------------------------------------
fn select_instances(
    consensi_path: &PathBuf,
    family_name: &str,
    from_path: &PathBuf,
    to_path: &PathBuf,
    working_dir: &Path,
    alphabet: &SequenceAlphabet,
    config: &Config,
    num_threads: usize,
) -> Result<usize> {
    let min_num_instances = match alphabet {
        SequenceAlphabet::Dna => config.general.min_num_instances_dna,
        _ => config.general.min_num_instances_prot,
    };

    if to_path.is_file() {
        let mut reader = FastaReader::new(open(to_path)?);
        let num = reader.records().count();
        if num < min_num_instances {
            debug!(
                "Removing previous instance file {}, too few ({num})",
                to_path.display()
            );
            fs::remove_file(to_path)?;
            return Ok(0);
        } else {
            debug!("Reusing existing instance file: {}", to_path.display());
            return Ok(num);
        }
    }

    let blast_dir = working_dir.join(family_name);
    fs::create_dir_all(&blast_dir)?;

    // Extract the family's consensus sequence
    let db_path = blast_dir.join("db.fa");
    if !db_path.exists() {
        let mut writer = FastaWriter::new(BufWriter::new(open_for_write(&db_path)?));
        let num_taken =
            copy_fasta(&[family_name.to_string()], consensi_path, &mut writer)?;

        if num_taken != 1 {
            fs::remove_file(db_path)?;
            bail!(
                "Failed to find family '{family_name}' in consensi '{}'",
                consensi_path.display()
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
    let blast_out = match alphabet {
        SequenceAlphabet::Dna => {
            run_rmblastn(&blast_dir, &db_path, from_path, config, num_threads)?
        }
        _ => run_blastp(&blast_dir, &db_path, from_path, config, num_threads)?,
    };
    let alignments = parse_alignment(&blast_out)?;

    // Remove short alignments, find the highest score for each hit
    let mut targets: HashMap<String, usize> = HashMap::new();
    let min_len = match alphabet {
        SequenceAlphabet::Dna => config.general.min_instance_sequence_length_dna,
        _ => config.general.min_instance_sequence_length_prot,
    };

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
        if wanted.len() >= config.general.max_num_instances && coverage_reached
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
            .all(|val| *val >= config.general.min_consensus_coverage);
    }

    let mut num_taken = 0;
    if !wanted.is_empty() {
        let mut reader = FastaReader::new(BufReader::new(open(from_path)?));
        let mut fasta_writer =
            FastaWriter::new(BufWriter::new(open_for_write(to_path)?));
        for record in reader.records().map_while(Result::ok) {
            let name = String::from_utf8(record.name().to_vec())?;
            if wanted.contains(&name) {
                fasta_writer.write_record(&record)?;
                num_taken += 1;
            }
        }
    }

    // Be sure to remove any empty files
    if num_taken < min_num_instances && to_path.exists() {
        debug!("Removing family '{family_name}', too few instances ({num_taken})");
        fs::remove_file(to_path)?;
        num_taken = 0;
    } else {
        debug!("Took {num_taken} instances for {family_name}");
    }

    Ok(num_taken)
}

// --------------------------------------------------
fn downsample(
    fasta: &PathBuf,
    num_wanted: usize,
    rev_comp: bool,
    start_numbering_at: usize,
    mut output: impl Write,
) -> Result<usize> {
    let mut reader = parse_reader(open(fasta)?)?;
    let mut num_taken = 0;

    debug!(
        "Taking {num_wanted} from '{}' ({})",
        fasta.display(),
        if rev_comp { "RevComp" } else { "Normal" }
    );

    while let Some(rec) = reader.iter_record()? {
        if num_taken == num_wanted {
            break;
        }
        writeln!(
            output,
            ">{} {}{}\n{}",
            start_numbering_at + num_taken,
            rec.head(),
            rec.des(),
            if rev_comp {
                String::from_utf8(revcomp(rec.seq().as_bytes()))?
            } else {
                rec.seq().to_string()
            }
        )?;
        num_taken += 1;
    }

    Ok(num_taken)
}

// --------------------------------------------------
fn extract_scores(
    alignment_file: &PathBuf,
    prev_scores_file: &Option<PathBuf>,
    consensi: &PathBuf,
    outdir: &Path,
) -> Result<PathBuf> {
    debug!(
        "Extracting scores from alignment '{}' {:?}",
        alignment_file.display(),
        prev_scores_file
    );

    let mut reader = parse_reader(open(consensi)?)?;
    let mut consensi_names: HashMap<String, String> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensi_names
            .insert(rec.head().trim().to_string(), rec.des().trim().to_string());
    }

    let scores_file = outdir.join("alignment-scores.tsv");
    let mut scores_wtr = WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&scores_file)
        .map_err(|e| anyhow!("Failed to write '{}': {e}", &scores_file.display()))?;

    // This is a hash to remember the family names we handled
    // in this round. When reading the scores from the previous
    // round, we need to skip these queries.
    let mut skip_query: HashSet<String> = HashSet::new();

    // BLAST output fails to include headers
    let mut align_reader = ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(alignment_file)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &alignment_file.display()))?;

    let records = align_reader.records();
    for res in records {
        let res = res?;
        let rec: RmBlastOutput = res.deserialize(None)?;
        match consensi_names.get(&rec.query) {
            Some(query_name) => {
                // Note the families involved in the previous round's
                // merges in order to skip them when adding the
                // previous round's scores.
                for family in parse_newick(query_name) {
                    skip_query.insert(family);
                }

                scores_wtr.serialize(AlignmentScore {
                    score: rec.score,
                    target: rec.target.to_string(),
                    query: query_name.clone(),
                })?
            }
            _ => eprintln!("Cannot find query '{}'", rec.query),
        }
    }

    // Add all the previous scores for queries not yet merged
    if let Some(prev_scores) = prev_scores_file {
        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(true)
            .from_path(prev_scores)
            .map_err(|e| anyhow!("Failed to read '{}': {e}", &prev_scores.display()))?;
        let records = reader.deserialize();
        for res in records {
            let rec: AlignmentScore = res?;

            // Move along if any of this sequence's families have been seen
            let query_families = parse_newick(&rec.query);
            if !query_families.into_iter().any(|v| skip_query.contains(&v)) {
                scores_wtr.serialize(rec)?;
            }
        }
    }

    Ok(scores_file)
}

// --------------------------------------------------
fn find_independence(num_wins: u32, num_shared: u32) -> f64 {
    if num_shared > 0 {
        num_wins as f64 / (num_wins as f64 + num_shared as f64)
    } else {
        1.
    }
}

// --------------------------------------------------
fn format_seconds(seconds: u64) -> String {
    let mut delta = Duration::seconds(seconds as i64);
    let mut ret = vec![];
    let days = delta.num_days();
    if days > 0 {
        ret.push(format!("{days} day{}", if days == 1 { "" } else { "s" }));
    }
    delta -= Duration::seconds(days * 24 * 60 * 60);

    let hours = delta.num_hours();
    if hours > 0 {
        ret.push(format!("{hours} hour{}", if hours == 1 { "" } else { "s" }));
    }
    delta -= Duration::seconds(hours * 60 * 60);

    let minutes = delta.num_minutes();
    if minutes > 0 {
        ret.push(format!(
            "{minutes} minute{}",
            if minutes == 1 { "" } else { "s" }
        ));
    }
    delta -= Duration::seconds(minutes * 60);

    let seconds = delta.num_seconds();
    if seconds > 0 || ret.is_empty() {
        ret.push(format!(
            "{seconds} second{}",
            if seconds == 1 { "" } else { "s" }
        ));
    }

    ret.join(", ")
}

// --------------------------------------------------
fn independence(winners: Winners) -> Vec<Independence> {
    let mut families = HashSet::<&str>::new();
    for StringPair(f1, f2) in winners.winning_sets.keys() {
        families.insert(f1);
        families.insert(f2);
    }
    let mut vals = vec![];

    for pair in families.into_iter().permutations(2) {
        if let [f1, f2] = pair[..] {
            let num_wins = *winners.clear_winners.get(f1).unwrap_or(&0);
            let key = StringPair::new(f1.to_string(), f2.to_string());
            let &num_shared = winners.winning_sets.get(&key).unwrap_or(&0u32);
            let ind = find_independence(num_wins, num_shared);

            vals.push(Independence {
                f1: f1.to_string(),
                f2: f2.to_string(),
                val: ind,
            });
        }
    }

    // Sort from least independent to most
    vals.sort_by(|a, b| {
        a.val
            .partial_cmp(&b.val)
            .unwrap()
            .then_with(|| a.f1.cmp(&b.f1))
            .then_with(|| a.f2.cmp(&b.f2))
    });

    vals
}

// --------------------------------------------------
fn merge_families(
    args: MergeFamilies,
    family_to_instance: &mut HashMap<String, PathBuf>,
) -> Result<String> {
    fs::create_dir_all(&args.outdir)?;

    let fams1 = parse_newick(&args.family1);
    let fams2 = parse_newick(&args.family2);
    let num_fams1 = fams1.len() as f64;
    let num_fams2 = fams2.len() as f64;
    let num_fams_total = num_fams1 + num_fams2;
    let num_seqs_total = 100;
    let num_from1 =
        (num_seqs_total as f64 * (num_fams1 / num_fams_total)).round() as usize;
    let num_from2 =
        (num_seqs_total as f64 * (num_fams2 / num_fams_total)).round() as usize;
    let f1 = fams1.join("::");
    let f2 = fams2.join("::");
    let new_family_name = format!("{f1}::{f2}");
    let new_family_path = tempfile::Builder::new()
        .prefix("inst-")
        .suffix(".fa")
        //.disable_cleanup(true)
        .tempfile_in(args.taken_instances_dir)?
        .path()
        .to_path_buf();

    debug!(
        "Merging {num_from1} from {f1}, {num_from2} from {f2} => {}",
        new_family_path.display()
    );

    // Block to isolate "output" and force close when passing out of scope
    {
        let mut output = open_for_write(&new_family_path)?;
        let mut total_taken: usize = 0;
        let mut one_flipped = false;
        for (fam, num) in &[(f1, num_from1), (f2, num_from2)] {
            let fasta = family_to_instance
                .get(fam)
                .unwrap_or_else(|| panic!("Missing instances for family '{fam}'"));
            let flip = if args.flipped && !one_flipped {
                one_flipped = true;
                true
            } else {
                false
            };
            total_taken += downsample(fasta, *num, flip, total_taken, &mut output)?;
        }

        if total_taken == 0 {
            bail!("Failed to extract any sequences for MSA");
        }
    }

    // Remember this for the next iteration
    let _ = &family_to_instance.insert(new_family_name, new_family_path.clone());

    // Copy the file to the current merge dir
    let msa_input = args.outdir.join("msa-input.fa");
    fs::copy(new_family_path, &msa_input)?;

    let consensus_path = match args.alphabet {
        SequenceAlphabet::Dna => msa_dna(&msa_input, args.num_threads)?,
        _ => msa_protein(&msa_input, args.num_threads)?,
    };

    let mut reader = parse_reader(open(&consensus_path)?)?;
    let consensus_seq = reader
        .iter_record()?
        .map(|rec| rec.seq().to_string())
        .expect("Failed to read consensus file");

    Ok(consensus_seq)
}

// --------------------------------------------------
fn msa_protein(input_file: &Path, num_threads: usize) -> Result<PathBuf> {
    let mafft = which("mafft").map_err(|e| anyhow!("mafft: {e}"))?;
    let hmmbuild = which("hmmbuild").map_err(|e| anyhow!("hmmbuild: {e}"))?;
    let hmmemit = which("hmmemit").map_err(|e| anyhow!("hmmemit: {e}"))?;

    // mafft creates MSA
    let start = Instant::now();
    let mut mafft_cmd = std::process::Command::new(&mafft);
    let mafft_args = &[
        "--auto".to_string(),
        "--thread".to_string(),
        num_threads.to_string(),
        input_file.to_string_lossy().to_string(),
    ];
    debug!(r#"Running "{} {}""#, mafft.display(), mafft_args.join(" "));
    let res = mafft_cmd.args(mafft_args).output()?;

    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    debug!(
        "Mafft finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    let outdir = input_file.parent().unwrap_or_else(|| {
        panic!("Failed to get parent dir for {}", input_file.display())
    });

    let msa_path = outdir.join("msa.fa");
    {
        let mut file = open_for_write(&msa_path)?;
        write!(file, "{}", String::from_utf8(res.stdout)?)?;
    }

    // hmmerbuild create HMM
    let hmm_out = outdir.join("hmm.out");
    let mut hmmbuild_cmd = std::process::Command::new(&hmmbuild);
    let hmmbuild_args = &[
        hmm_out.to_string_lossy().to_string(),
        msa_path.to_string_lossy().to_string(),
    ];
    let res = hmmbuild_cmd.args(hmmbuild_args).output()?;
    debug!(
        r#"Running "{} {}""#,
        hmmbuild.display(),
        hmmbuild_args.join(" ")
    );

    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    // hmmeremit generates a consensus sequence
    let mut hmmemit_cmd = std::process::Command::new(&hmmemit);
    let hmmemit_args = &["-c".to_string(), hmm_out.to_string_lossy().to_string()];
    let res = hmmemit_cmd.args(hmmemit_args).output()?;
    debug!(
        r#"Running "{} {}""#,
        hmmemit.display(),
        hmmemit_args.join(" ")
    );

    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    let consensus_path = outdir.join("consensus.fa");
    {
        let mut file = open_for_write(&consensus_path)?;
        write!(file, "{}", String::from_utf8(res.stdout)?)?;
    }

    Ok(consensus_path)
}

// --------------------------------------------------
fn msa_dna(input_file: &Path, num_threads: usize) -> Result<PathBuf> {
    let refiner = which("Refiner").map_err(|e| anyhow!("Refiner: {e}"))?;

    let mut refiner_args = vec![
        //"-debug".to_string(),
        "-threads".to_string(),
        num_threads.to_string(),
    ];

    //let rmblast = which("rmblastn").map_err(|e| anyhow!("rmblastn: {e}"))?;
    //if let Some(rmblast_dir) = rmblast
    //    .as_path()
    //    .parent()
    //    .map(|path| path.to_string_lossy().to_string())
    //{
    //    refiner_args.extend_from_slice(&["--rmblast_dir".to_string(), rmblast_dir]);
    //}

    refiner_args.push(input_file.to_string_lossy().to_string());
    debug!(
        r#"Running "{} {}""#,
        &refiner.display(),
        &refiner_args.join(" ")
    );

    let start = Instant::now();
    let mut cmd = std::process::Command::new(&refiner);
    let res = cmd.args(refiner_args).output()?;
    if !res.status.success() {
        debug!("{}", String::from_utf8(res.stdout)?);
        bail!(String::from_utf8(res.stderr)?);
    }

    debug!(
        "Refiner finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    let outdir = input_file.parent().unwrap_or_else(|| {
        panic!("Failed to get parent dir for {}", input_file.display())
    });
    let consensus_path = outdir.join("msa-input.fa.refiner_cons");
    if !consensus_path.exists() {
        bail!(
            "Failed to find expected consensus file {}",
            consensus_path.display()
        );
    }

    Ok(consensus_path)
}

// --------------------------------------------------
fn number_consensi(
    consensi: &PathBuf,
    outpath: &PathBuf,
) -> Result<HashMap<String, String>> {
    let mut outfile = open_for_write(outpath)?;
    let mut reader = parse_reader(open(consensi)?)?;
    let mut seqs: HashMap<String, String> = HashMap::new();
    let mut i = 0;

    while let Some(rec) = reader.iter_record()? {
        writeln!(outfile, ">{i} {}\n{}", rec.head(), rec.seq())?;
        seqs.insert(rec.head().to_string(), rec.seq().to_string());
        i += 1;
    }

    Ok(seqs)
}

// --------------------------------------------------
fn open(filename: &PathBuf) -> Result<Box<dyn BufRead>> {
    Ok(Box::new(BufReader::new(File::open(filename).map_err(
        |e| anyhow!("Cannot read {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn open_for_write(filename: &PathBuf) -> Result<Box<dyn Write>> {
    Ok(Box::new(BufWriter::new(File::create(filename).map_err(
        |e| anyhow!("Cannot write {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
fn read_lines(path: &PathBuf) -> Result<Vec<String>> {
    Ok(open(path)?
        .lines()
        .map_while(Result::ok)
        .filter(|line| !line.is_empty())
        .collect())
}

// --------------------------------------------------
fn parse_newick(val: &str) -> Vec<String> {
    let mut ret = vec![];
    match newick::from_string(val) {
        // Incoming value is only Newick when merged
        Ok(trees) => {
            let mut leaves: Vec<_> = trees
                .into_iter()
                .flat_map(|tree| {
                    let leaves = tree.leaves().collect::<Vec<_>>();
                    leaves.into_iter().flat_map(move |leaf| {
                        tree.name(leaf).map(|name| name.to_string())
                    })
                })
                .collect();
            ret.append(&mut leaves);
        }
        // Otherwise return the original string
        Err(_) => ret.push(val.to_string()),
    }

    ret
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use crate::{parse_alignment, RmBlastOutput};

    use super::{
        bitscore_to_confidence, call_winners, cat_sequences, downsample,
        extract_scores, find_independence, format_seconds, independence,
        number_consensi, open, parse_newick, Independence, StringPair, Winners,
    };
    use anyhow::Result;
    use kseq::parse_reader;
    use pretty_assertions::assert_eq;
    use std::{
        collections::HashMap,
        fs::{self, File},
        path::PathBuf,
    };
    use tempfile::{tempdir, NamedTempFile};

    #[test]
    fn test_bitscore_to_confidence() -> Result<()> {
        let res = bitscore_to_confidence(&[&2274, &2234, &2224, &2245, &2295], 0.1227);
        assert!(res.is_ok());
        assert_eq!(
            res.unwrap(),
            [
                0.14088157852127628,
                0.004692442832468759,
                0.002004634432314714,
                0.011959118673026447,
                0.8404622255409138,
            ]
        );

        Ok(())
    }

    #[test]
    fn test_call_winners() -> Result<()> {
        let scores_file = PathBuf::from("tests/outputs/alignment-scores.tsv");
        let lambda = 0.1227;
        let confidence_margin = 3;
        let res = call_winners(&scores_file, lambda, confidence_margin);
        assert!(res.is_ok());

        let winners = res.unwrap();
        let expected_wins: Vec<(&str, u32)> = vec![
            ("AluY", 95),
            ("AluYm1", 72),
            ("AluYb9", 93),
            ("AluYb8", 103),
            ("AluYa5", 98),
        ];

        for (key, val) in expected_wins.iter() {
            assert_eq!(winners.clear_winners.get(*key), Some(val));
        }

        // Winning sets should be symmetrical (A/B, B/A)
        let expected_winning_sets: HashMap<StringPair, u32> = HashMap::from([
            (StringPair("AluYm1".to_string(), "AluY".to_string()), 28),
            (StringPair("AluY".to_string(), "AluYm1".to_string()), 28),
            //
            (StringPair("AluYb8".to_string(), "AluYa5".to_string()), 3),
            (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 3),
            //
            (StringPair("AluYb8".to_string(), "AluYb9".to_string()), 4),
            (StringPair("AluYb9".to_string(), "AluYb8".to_string()), 4),
            //
            (StringPair("AluY".to_string(), "AluYb8".to_string()), 3),
            (StringPair("AluYb8".to_string(), "AluY".to_string()), 3),
            //
            (StringPair("AluY".to_string(), "AluYa5".to_string()), 7),
            (StringPair("AluYa5".to_string(), "AluY".to_string()), 7),
        ]);

        for (pair, score) in expected_winning_sets {
            let res = winners.winning_sets.get(&pair);
            assert!(res.is_some());
            assert_eq!(res.unwrap(), &score);
        }

        Ok(())
    }

    #[test]
    fn test_cat_sequences() -> Result<()> {
        let instances_dir = PathBuf::from("tests/inputs/instances_100");
        let outdir = tempdir()?;

        let outpath = outdir.path().join("all_seqs.fa");
        let res = cat_sequences(
            &instances_dir,
            &[
                "AluY".to_string(),
                "AluYa5".to_string(),
                "AluYb8".to_string(),
                "AluYb9".to_string(),
                "AluYm1".to_string(),
            ],
            &outpath,
        );
        assert!(res.is_ok());
        assert!(outpath.exists());

        let mut reader = parse_reader(open(&outpath)?)?;
        let mut count = 0;
        while (reader.iter_record()?).is_some() {
            count += 1;
        }
        assert_eq!(count, 500);

        Ok(())
    }

    //#[test]
    //fn test_check_family_instances() -> Result<()> {
    //    // The consensi contains duplicated IDs
    //    let consensi = PathBuf::from("tests/inputs/dup_consensi.fa");
    //    let instances_100 = &[PathBuf::from(
    //        "tests/inputs/instances_100/AluY.fa".to_string(),
    //    )];
    //    let res = check_family_instances(&consensi, instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "The following consensi IDs are duplicated: AluYb9, AluYm1"
    //    );
    //
    //    // The consensi has 5 unique IDs but the instances only 2
    //    let consensi = PathBuf::from("tests/inputs/consensi.fa");
    //    let instances_100 = &[
    //        PathBuf::from("tests/inputs/AluY.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
    //    ];
    //    let res = check_family_instances(&consensi, instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "Missing the following instances: AluYa5, AluYb8, AluYb9",
    //    );
    //
    //    // The consensi has 2 unique IDs but the instances have 4
    //    let consensi = PathBuf::from("tests/inputs/two_consensi.fa");
    //    let instances_100 = vec![
    //        PathBuf::from("tests/inputs/AluY.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYa5.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYb9.fa".to_string()),
    //        PathBuf::from("tests/inputs/AluYm1.fa".to_string()),
    //    ];
    //    let res = check_family_instances(&consensi, &instances_100);
    //    assert!(res.is_err());
    //    assert_eq!(
    //        res.unwrap_err().to_string(),
    //        "Missing the following consensi: AluYb9, AluYm1",
    //    );
    //
    //    Ok(())
    //}

    #[test]
    fn test_downsample() -> Result<()> {
        let fasta = PathBuf::from("tests/inputs/AluY_5.fa");
        let tmp_dir = tempdir()?;
        let outpath = tmp_dir.path().join("sub.fa");
        let rev_comp = false;

        // Scoped to force close of output filehandle
        {
            let out = File::create(&outpath)?;

            // Select 3 of the 5 sequences
            let res = downsample(&fasta, 3, rev_comp, 0, &out);
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), 3);
        }

        {
            let sub = File::open(&outpath)?;
            let mut reader = parse_reader(sub)?;
            let mut num = 0;
            while (reader.iter_record()?).is_some() {
                num += 1;
            }

            // Ensure 3 sequences were written
            assert_eq!(num, 3);
        }

        // Try to select more than the 5 sequences
        // Should only get the actual 5
        {
            let out = File::create(&outpath)?;
            let res = downsample(&fasta, 10, rev_comp, 0, &out);
            assert!(res.is_ok());
            assert_eq!(res.unwrap(), 5);
        }

        {
            let sub = File::open(&outpath)?;
            let mut reader = parse_reader(sub)?;
            let mut num = 0;
            while (reader.iter_record()?).is_some() {
                num += 1;
            }

            // Ensure 5 sequences were written
            assert_eq!(num, 5);
        }

        Ok(())
    }

    #[test]
    fn test_extract_scores() -> Result<()> {
        let outdir = tempdir()?;
        let consensi = PathBuf::from("tests/inputs/numbered_consensi.fa".to_string());
        let alignment_file = PathBuf::from("tests/inputs/blast.tsv");
        let prev_scores: Option<PathBuf> = None;
        let res =
            extract_scores(&alignment_file, &prev_scores, &consensi, outdir.path());
        dbg!(&res);
        assert!(res.is_ok());

        let scores_file = res.unwrap();
        assert!(scores_file.exists());

        let actual = fs::read_to_string(scores_file)?;
        let expected = fs::read_to_string("tests/inputs/alignment-scores.tsv")?;
        assert_eq!(actual, expected);
        Ok(())
    }

    #[test]
    fn test_find_independence() -> Result<()> {
        // No shared wins means fully independent
        //                       wins  shared
        let res = find_independence(1, 0);
        assert_eq!(res, 1.);

        // One clear winner, one shared win == 50%
        //                       wins  shared
        let res = find_independence(1, 1);
        assert_eq!(res, 0.5);

        // No clear winner, only shared wins == 0%
        //                       wins  shared
        let res = find_independence(0, 2);
        assert_eq!(res, 0.);

        //                        wins  shared
        let res = find_independence(37, 46);
        assert_eq!(res, 0.4457831325301205);
        Ok(())
    }

    #[test]
    fn test_format_seconds() -> Result<()> {
        let one_hour = 60 * 60;
        let one_day = one_hour * 24;
        assert_eq!(format_seconds(0), "0 seconds");
        assert_eq!(format_seconds(1), "1 second");
        assert_eq!(format_seconds(59), "59 seconds");
        assert_eq!(format_seconds(60), "1 minute");
        assert_eq!(format_seconds(120), "2 minutes");
        assert_eq!(format_seconds(121), "2 minutes, 1 second");
        assert_eq!(format_seconds(one_hour), "1 hour");
        assert_eq!(format_seconds(one_hour + 1), "1 hour, 1 second");
        assert_eq!(
            format_seconds(one_hour + 121),
            "1 hour, 2 minutes, 1 second"
        );
        assert_eq!(format_seconds((one_hour * 4) + 59), "4 hours, 59 seconds");
        assert_eq!(format_seconds(one_day), "1 day");
        assert_eq!(format_seconds(one_day + 2), "1 day, 2 seconds");
        Ok(())
    }

    #[test]
    fn test_independence() -> Result<()> {
        let clear_winners = HashMap::from([
            ("AluY".to_string(), 37),
            ("AluYm1".to_string(), 36),
            ("AluYb8".to_string(), 50),
        ]);

        // Winning sets are expected to have symmetrical keys
        let winning_sets = HashMap::from([
            (StringPair("AluYa5".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluYb8".to_string(), "AluYa5".to_string()), 6),
            //
            (StringPair("AluYa5".to_string(), "AluYm1".to_string()), 4),
            (StringPair("AluYm1".to_string(), "AluYa5".to_string()), 4),
            //
            (StringPair("AluY".to_string(), "AluYb8".to_string()), 6),
            (StringPair("AluYb8".to_string(), "AluY".to_string()), 6),
            //
            (StringPair("AluY".to_string(), "AluYa5".to_string()), 12),
            (StringPair("AluYa5".to_string(), "AluY".to_string()), 12),
            //
            (StringPair("AluY".to_string(), "AluYm1".to_string()), 46),
            (StringPair("AluYm1".to_string(), "AluY".to_string()), 46),
        ]);

        let ind = independence(Winners {
            clear_winners,
            winning_sets,
        });

        let expected = [
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluY".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluYb8".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYa5".to_string(),
                f2: "AluYm1".to_string(),
                val: 0.0,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluY".to_string(),
                val: 0.43902439024390244,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYm1".to_string(),
                val: 0.4457831325301205,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.7551020408163265,
            },
            Independence {
                f1: "AluY".to_string(),
                f2: "AluYb8".to_string(),
                val: 0.8604651162790697,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluY".to_string(),
                val: 0.8928571428571429,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.8928571428571429,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluYa5".to_string(),
                val: 0.9,
            },
            Independence {
                f1: "AluYb8".to_string(),
                f2: "AluYm1".to_string(),
                val: 1.0,
            },
            Independence {
                f1: "AluYm1".to_string(),
                f2: "AluYb8".to_string(),
                val: 1.0,
            },
        ];
        assert_eq!(ind, expected);
        Ok(())
    }

    #[test]
    fn test_number_consensi() -> Result<()> {
        let orig_consensi = PathBuf::from("tests/inputs/consensi.fa");
        let outpath = NamedTempFile::new()?;
        let res = number_consensi(&orig_consensi, &outpath.path().into());
        assert!(res.is_ok());

        let actual = fs::read_to_string(outpath)?;
        let expected = fs::read_to_string("tests/outputs/numbered_consensi.fa")?;
        assert_eq!(actual, expected);

        Ok(())
    }

    #[test]
    fn test_parse_newick() -> Result<()> {
        let res = parse_newick("A");
        assert_eq!(res, vec!["A"]);

        let res = parse_newick("(((B,D):0.1,A):0.3,C):0.2;");
        assert_eq!(res, vec!["B", "D", "A", "C"]);

        let res = parse_newick("((X,Y):0.04,Z):0.25;");
        assert_eq!(res, vec!["X", "Y", "Z"]);
        Ok(())
    }

    #[test]
    fn test_parse_alignment() -> Result<()> {
        let path = PathBuf::from("./tests/inputs/blast-consensi-self.tsv");
        let res = parse_alignment(&path);
        assert!(res.is_ok());

        let alignments = res.unwrap();
        assert_eq!(alignments.len(), 1000);

        let first = alignments.first().unwrap();
        assert_eq!(
            first,
            &RmBlastOutput {
                score: 710,
                target: "Chompy-2a_tua".to_string(),
                query: "Chompy-2a_tua".to_string(),
                query_len: 79,
                query_start: 1,
                query_end: 79,
                subject_len: 79,
                subject_start: 1,
                subject_end: 79,
                cpg_kdiv: 0.0,
                pident: 100.0,
            }
        );

        let last = alignments.last().unwrap();
        assert_eq!(
            last,
            &RmBlastOutput {
                score: 409,
                target: "Harbinger-3-L_tua".to_string(),
                query: "tuafam018707_consensus".to_string(),
                query_len: 6639,
                query_start: 5945,
                query_end: 6261,
                subject_len: 1898,
                subject_start: 173,
                subject_end: 497,
                cpg_kdiv: 41.17,
                pident: 54.571,
            }
        );

        Ok(())
    }
}
