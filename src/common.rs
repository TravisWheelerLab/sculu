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
    process::{self, Command},
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
    run_blast(BlastParams {
        out_dir,
        db,
        query,
        num_threads,
        db_type: "nucl",
        executable: "rmblastn",
        gap_open: config.rmblastn.gap_open,
        gap_extend: config.rmblastn.gap_extend,
        word_size: config.rmblastn.word_size,
        mask_level: config.rmblastn.mask_level,
        min_raw_gapped_score: config.rmblastn.min_raw_gapped_score,
        dust: config.rmblastn.dust,
        complexity_adjust: config.rmblastn.complexity_adjust,
        matrix: config.rmblastn.matrix.as_ref(),
        run_in_out_dir: true,
    })
}

// --------------------------------------------------
pub fn run_blastp(
    out_dir: &PathBuf,
    db: &Path,
    query: &Path,
    config: &Config,
    num_threads: usize,
) -> Result<PathBuf> {
    run_blast(BlastParams {
        out_dir,
        db,
        query,
        num_threads,
        db_type: "prot",
        executable: "blastp",
        gap_open: config.blastp.gap_open,
        gap_extend: config.blastp.gap_extend,
        word_size: config.blastp.word_size,
        mask_level: config.blastp.mask_level,
        min_raw_gapped_score: config.blastp.min_raw_gapped_score,
        dust: config.blastp.dust,
        complexity_adjust: config.blastp.complexity_adjust,
        matrix: config.blastp.matrix.as_ref(),
        run_in_out_dir: false,
    })
}

// --------------------------------------------------
struct BlastParams<'a> {
    out_dir: &'a PathBuf,
    db: &'a Path,
    query: &'a Path,
    num_threads: usize,
    db_type: &'static str,
    executable: &'static str,
    gap_open: usize,
    gap_extend: usize,
    word_size: usize,
    mask_level: usize,
    min_raw_gapped_score: usize,
    dust: bool,
    complexity_adjust: bool,
    matrix: Option<&'a PathBuf>,
    // rmblastn must run from out_dir so relative matrix paths resolve correctly
    run_in_out_dir: bool,
}

fn run_blast(p: BlastParams) -> Result<PathBuf> {
    fs::create_dir_all(p.out_dir)?;
    let outfile = p.out_dir.join("blast.tsv");

    if outfile.exists() {
        debug!("Reusing BLAST output file '{}'", outfile.display());
        return Ok(outfile);
    }

    let db_path = p.out_dir.join("db");
    let makeblastdb = which("makeblastdb").map_err(|e| anyhow!("makeblastdb: {e}"))?;
    let mut cmd = Command::new(makeblastdb);
    cmd.args([
        "-out",
        &db_path.to_string_lossy(),
        "-parse_seqids",
        "-dbtype",
        p.db_type,
        "-in",
        &p.db.to_string_lossy(),
    ]);
    let _ = run_cmd(cmd)?;

    let blast_exe = which(p.executable).map_err(|e| anyhow!("{}: {e}", p.executable))?;
    let mut cmd = Command::new(blast_exe);
    let mut args = vec![
        "-db".to_string(),
        db_path.to_string_lossy().to_string(),
        "-query".to_string(),
        p.query.to_string_lossy().to_string(),
        "-out".to_string(),
        outfile.to_string_lossy().to_string(),
        "-outfmt".to_string(),
        "6 score qseqid sseqid qlen qstart qend slen sstart send cpg_kdiv pident".to_string(),
        "-num_threads".to_string(),
        p.num_threads.to_string(),
        "-num_alignments".to_string(),
        "9999999".to_string(),
        "-gapopen".to_string(),
        p.gap_open.to_string(),
        "-gapextend".to_string(),
        p.gap_extend.to_string(),
        "-word_size".to_string(),
        p.word_size.to_string(),
        "-mask_level".to_string(),
        p.mask_level.to_string(),
        "-min_raw_gapped_score".to_string(),
        p.min_raw_gapped_score.to_string(),
        "-xdrop_ungap".to_string(),
        (p.min_raw_gapped_score * 2).to_string(),
        "-xdrop_gap".to_string(),
        (p.min_raw_gapped_score / 2).to_string(),
        "-xdrop_gap_final".to_string(),
        p.min_raw_gapped_score.to_string(),
        "-dust".to_string(),
        if p.dust { "yes" } else { "no" }.to_string(),
    ];

    if p.complexity_adjust {
        args.push("-complexity_adjust".to_string());
    }

    if let Some(matrix) = p.matrix {
        if !matrix.is_file() {
            bail!("{} matrix '{}' does not exist", p.executable, matrix.display());
        }
        let matrix_filename = matrix
            .file_name()
            .ok_or_else(|| anyhow!("Failed to get filename from '{}'", matrix.display()))?;
        let matrix_dir = matrix
            .parent()
            .ok_or_else(|| anyhow!("Failed to get dirname from '{}'", matrix.display()))?;
        cmd.env("BLASTMAT", matrix_dir);
        args.extend_from_slice(&[
            "-matrix".to_string(),
            matrix_filename.to_string_lossy().to_string(),
        ]);
    }

    if p.run_in_out_dir {
        cmd.current_dir(p.out_dir);
    }
    cmd.args(args);
    let _ = run_cmd(cmd)?;

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

// --------------------------------------------------
pub fn run_cmd(mut cmd: Command) -> Result<process::Output> {
    let start = Instant::now();
    let res = cmd.output()?;

    if !res.status.success() {
        bail!(String::from_utf8(res.stderr)?);
    }

    debug!(
        "Command {cmd:?} finished in {}",
        format_seconds(start.elapsed().as_secs())
    );

    Ok(res)
}

// --------------------------------------------------
#[cfg(test)]
mod tests {
    use super::format_seconds;
    use anyhow::Result;

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
}
