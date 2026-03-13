use crate::types::{BlastpConfig, Config, GeneralConfig, RmblastnConfig};
use anyhow::{anyhow, bail, Result};
use chrono::Duration;
use log::debug;
use noodles_fasta::{self, io::Reader as FastaReader, io::Writer as FastaWriter};
use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use which::which;

// --------------------------------------------------
pub fn read_instances_dir(instances_dir: &Path) -> Result<HashMap<String, PathBuf>> {
    // Create hashmap from family name to taken instance file
    let mut family_to_instance: HashMap<String, PathBuf> = HashMap::new();
    let entries = fs::read_dir(instances_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, instances_dir.display()))?;

    for entry in entries {
        let entry = entry?;
        if let Some(stem) = entry.path().file_stem() {
            let family_name = stem.to_string_lossy().to_string();
            if !family_name.starts_with(".") && !family_name.starts_with("inst-") {
                family_to_instance.insert(family_name, entry.path().to_path_buf());
            }
        }
    }
    Ok(family_to_instance)
}

// --------------------------------------------------
pub fn copy_fasta<W: Write>(
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
pub fn open(filename: &PathBuf) -> Result<Box<dyn BufRead>> {
    Ok(Box::new(BufReader::new(File::open(filename).map_err(
        |e| anyhow!("Cannot read {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
pub fn open_for_write(filename: &PathBuf) -> Result<Box<dyn Write>> {
    Ok(Box::new(BufWriter::new(File::create(filename).map_err(
        |e| anyhow!("Cannot write {}: {e}", filename.display()),
    )?)))
}

// --------------------------------------------------
pub fn read_lines(path: &PathBuf) -> Result<Vec<String>> {
    Ok(open(path)?
        .lines()
        .map_while(Result::ok)
        .filter(|line| !line.is_empty())
        .collect())
}

// --------------------------------------------------
pub fn run_rmblastn(
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
pub fn run_blastp(
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
pub fn format_seconds(seconds: u64) -> String {
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
pub fn default_config() -> Config {
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
/// Read the (optional) "sculu.toml" file
pub fn get_config(config_file: &Option<PathBuf>) -> Result<Config> {
    match config_file {
        Some(filename) => {
            let mut file = open(filename)?;
            let mut contents = String::new();
            let _bytes = file.read_to_string(&mut contents);
            let config: Config = toml::from_str(&contents)?;
            Ok(config)
        }
        _ => Ok(default_config()),
    }
}
