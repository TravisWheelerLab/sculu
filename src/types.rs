use clap::{builder::PossibleValue, Parser, ValueEnum};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fmt, path::PathBuf};

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

#[derive(Parser, Debug)]
pub enum Command {
    /// Generate config TOML
    Config(ConfigArgs),

    /// Align consensus sequences to self and create "components"
    Components(ComponentsArgs),

    /// Cluster output from "components"
    Cluster(ClusterArgs),

    /// Concatenate the singletons and clusters from "components"
    Concat(ConcatArgs),

    /// Run all steps (build components, cluster, concat)
    Run(RunArgs),
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
    /// The filtered FASTA consensus from "components"
    #[arg(long, value_name = "CONSENSUS", required = true)]
    pub consensus_path: PathBuf,

    /// Merged components from "cluster"
    #[arg(long, value_name = "COMPONENTS", num_args = 0..)]
    pub components: Vec<PathBuf>,

    /// Singletons file from "components"
    #[arg(long, value_name = "SINGLETONS")]
    pub singletons: Option<PathBuf>,

    /// Output directory
    #[arg(short, long, value_name = "OUTDIR", default_value = "sculu-out")]
    pub outdir: PathBuf,
}

#[derive(Debug, Parser, Clone)]
#[command()]
pub struct RunArgs {
    /// FASTA file of subfamily consensus
    #[arg(long, value_name = "CONSENSUS", required = true)]
    pub consensus: PathBuf,

    /// Directory of seed alignments (MSA) for each subfamily
    #[arg(long, value_name = "ALIGNMENTS", required = true)]
    pub alignments: PathBuf,

    /// Sequence alphabet
    #[arg(short, long, value_name = "ALPHABET", required = true)]
    pub alphabet: SequenceAlphabet,

    /// Output directory
    #[arg(long, value_name = "OUTDIR", default_value = "sculu-out")]
    pub outdir: PathBuf,

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

    /// FASTA file of subfamily consensus sequences
    #[arg(long, value_name = "CONSENSUS", required = true)]
    pub consensus: PathBuf,

    /// Directory of seed alignments (MSA) for each subfamily
    #[arg(long, value_name = "ALIGNMENTS", required = true)]
    pub alignments: PathBuf,

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

    /// The filtered FASTA consensus sequences from "components"
    #[arg(long, value_name = "CONSENSUS", required = true)]
    pub consensus: PathBuf,

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
    pub general: GeneralConfig,
    pub rmblastn: RmblastnConfig,
    pub blastp: BlastpConfig,
}

/// "general" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct GeneralConfig {
    pub confidence_margin: isize,
    pub independence_threshold: f64,
    pub lambda: f64,
    pub percent_id_for_components: f64,
    pub max_num_instances: usize,
    pub min_align_cover: f64,
    pub min_consensus_coverage: usize,
    pub min_instance_sequence_length_dna: usize,
    pub min_instance_sequence_length_prot: usize,
    pub min_len_similarity: f64,
    pub min_num_instances_dna: usize,
    pub min_num_instances_prot: usize,
}

/// "blastp" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BlastpConfig {
    pub matrix: Option<PathBuf>,
    pub gap_open: usize,
    pub gap_extend: usize,
    pub word_size: usize,
    pub mask_level: usize,
    pub dust: bool,
    pub complexity_adjust: bool,
    pub min_raw_gapped_score: usize,
    pub xdrop_gap: usize,
    pub xdrop_ungap: usize,
    pub xdrop_gap_final: usize,
}

/// "rmblast" section of sculu.toml
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RmblastnConfig {
    pub matrix: Option<PathBuf>,
    pub gap_open: usize,
    pub gap_extend: usize,
    pub word_size: usize,
    pub mask_level: usize,
    pub dust: bool,
    pub complexity_adjust: bool,
    pub min_raw_gapped_score: usize,
    pub xdrop_gap: usize,
    pub xdrop_ungap: usize,
    pub xdrop_gap_final: usize,
}

#[derive(Debug, PartialEq)]
pub enum Partition {
    Top,
    Bottom,
}

#[derive(Debug, Deserialize, Serialize, Clone, PartialEq)]
pub struct RmBlastOutput {
    pub score: usize,
    pub target: String,
    pub query: String,
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub subject_len: usize,
    pub subject_start: usize,
    pub subject_end: usize,
    pub cpg_kdiv: f64,
    pub pident: f64,
}

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct AlignedConsensusPair {
    pub target: String,
    pub query: String,
    pub is_flipped: bool,
}

#[derive(Debug, PartialEq, Clone, Serialize, Deserialize)]
pub enum Direction {
    #[serde(rename = "forward")]
    Forward,
    #[serde(rename = "reverse")]
    Reverse,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct AlignmentScore {
    pub score: usize,
    pub target: String,
    pub query: String,
}

#[derive(Debug, PartialEq)]
pub struct Independence {
    pub f1: String,
    pub f2: String,
    pub val: f64,
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
    pub consensus_path: PathBuf,
    pub instances_dir: PathBuf,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct StringPair(pub String, pub String);

impl StringPair {
    // Stores the strings in ascending order
    pub fn new(s1: String, s2: String) -> StringPair {
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
pub struct Winners {
    pub clear_winners: HashMap<String, u32>,
    pub winning_sets: HashMap<StringPair, u32>,
}

#[derive(Debug)]
pub struct MergeFamilies<'a> {
    pub family1: String,
    pub family2: String,
    pub outdir: PathBuf,
    pub taken_instances_dir: &'a PathBuf,
    pub num_threads: usize,
    pub alphabet: SequenceAlphabet,
    pub flipped: bool,
}

#[derive(Debug)]
pub struct CheckFamilyInstancesArgs<'a> {
    pub consensus_path: &'a PathBuf,
    pub instances_dir: &'a PathBuf,
    pub out_dir: &'a PathBuf,
    pub taken_instances_dir: &'a PathBuf,
    pub config: &'a Config,
    pub alphabet: &'a SequenceAlphabet,
    pub num_threads: usize,
}

#[derive(Debug)]
pub struct CheckedFamilyResult {
    pub consensus_path: PathBuf,
    pub family_to_instance: HashMap<String, PathBuf>,
}

#[derive(Debug)]
pub struct SelectInstancesArgs<'a> {
    pub consensus_path: &'a PathBuf,
    pub family_name: &'a String,
    pub from_path: &'a PathBuf,
    pub to_path: &'a PathBuf,
    pub working_dir: &'a PathBuf,
    pub alphabet: &'a SequenceAlphabet,
    pub config: &'a Config,
    pub num_threads: usize,
}

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

#[derive(Debug)]
pub struct MsaResult {
    pub consensus_seq: String,
    pub consensus_path: PathBuf,
    pub msa_path: Option<PathBuf>,
}

#[derive(Debug)]
pub struct ClusterResult {
    pub consensus_path: PathBuf,
    pub family_to_msa: HashMap<String, PathBuf>,
}
