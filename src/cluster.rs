use crate::{
    common::{
        copy_fasta, format_seconds, get_config, open, open_for_write,
        read_instances_dir, read_lines, run_blastp, run_rmblastn,
    },
    types::{
        AlignedConsensusPair, AlignmentScore, ClusterArgs, Config, Independence,
        MergeFamilies, RmBlastOutput, SequenceAlphabet, StringPair, Winners,
    },
};
use anyhow::{anyhow, bail, Result};
use bio::alphabets::dna::revcomp;
use itertools::Itertools;
use kseq::parse_reader;
use log::debug;
use newick::Newick;
use noodles_fasta::{self, io::Reader as FastaReader, io::Writer as FastaWriter};
use regex::Regex;
use std::{
    cmp::max,
    collections::{HashMap, HashSet},
    fs,
    io::{BufWriter, Write},
    path::{Path, PathBuf},
    time::Instant,
};
use which::which;

// --------------------------------------------------
pub fn cluster_component(args: &ClusterArgs, num_threads: usize) -> Result<PathBuf> {
    let config = get_config(&args.config)?;

    if !&args.outdir.is_dir() {
        fs::create_dir_all(&args.outdir)?;
    }

    // If given an component file to process, assume the
    // consensus and instances have been properly filtered
    // and use given args as-is.
    let family_to_instance = read_instances_dir(&args.instances)?;

    let merged = merge_component(args, &config, family_to_instance, num_threads)?;

    debug!(
        "'{}' merged into '{}'",
        args.component.display(),
        merged.display()
    );

    Ok(merged)
}

// --------------------------------------------------
fn merge_component(
    args: &ClusterArgs,
    config: &Config,
    mut family_to_instance: HashMap<String, PathBuf>,
    num_threads: usize,
) -> Result<PathBuf> {
    debug!("Merging component '{}'", args.component.display());

    let filename_re = Regex::new(r"component-(\d+)").unwrap();
    let component_number =
        match filename_re.captures(&args.component.to_string_lossy().to_string()) {
            Some(caps) => caps.get(1).unwrap().as_str().to_string(),
            _ => bail!(r#"Component "{}" is not valid"#, args.component.display()),
        };

    // The "component-N" file will contain the names of the families
    let families = read_lines(&args.component.to_path_buf())?;

    // Get the alignments underlying this component
    let mut alignment_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(true)
        .from_path(&args.alignments)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &args.alignments.display()))?;

    let mut flipped_pairs: HashSet<StringPair> = HashSet::new();
    for res in alignment_reader.records() {
        let record: AlignedConsensusPair = res?.deserialize(None)?;
        let target = record.target.clone();
        let query = record.query.clone();
        if record.is_flipped
            && [&target, &query].iter().any(|name| families.contains(name))
        {
            let _ = flipped_pairs.insert(StringPair::new(target, query));
        }
    }

    // Create a working dir with the same name as the component file
    let component_name = args
        .component
        .file_name()
        .unwrap()
        .to_string_lossy()
        .to_string();
    let batch_dir = args.outdir.join(&component_name);
    fs::create_dir_all(&batch_dir)?;

    debug!(
        "==> Component '{component_name}' has {} families <==",
        families.len()
    );

    let start = Instant::now();

    // Create a consensus file for this round containing only the given families
    let mut batch_consensus = batch_dir.join("consensus.fa");
    {
        // Scoped to cause fasta_writer to close
        let mut fasta_writer =
            FastaWriter::new(BufWriter::new(open_for_write(&batch_consensus)?));
        let num_taken =
            copy_fasta(&families, &args.consensus.to_path_buf(), &mut fasta_writer)?;

        if num_taken == 0 {
            bail!(
                "Failed to copy consensus sequences from '{}' to '{}'",
                args.consensus.display(),
                batch_consensus.display()
            );
        } else {
            debug!(
                r#"Wrote {num_taken} batch consensus sequences to "{}""#,
                batch_consensus.display()
            );
        }
    }

    // Put the instances for the given families into a single file
    let all_seqs_path = &batch_dir.join("all_seqs.fa");
    let num_instances = cat_sequences(&args.instances, &families, all_seqs_path)?;

    debug!(
        "Concatenated {num_instances} instance sequences into '{}'",
        all_seqs_path.display()
    );

    // Create a starting directory for the merge iterations.
    let round0_dir = batch_dir.join("round000");
    fs::create_dir_all(&round0_dir)?;

    // As the consensus sequences are merged, the family names are concatenated
    // into Newick-formatted strings to track the merges.
    // These can get quite long, causing `makeblastdb` to fail.
    // Make the sequence IDs simple integers by numbering
    // and move the family names to the description.
    let consensus_path = round0_dir.join("consensus.fa");
    debug!(
        "Writing numbered consensus sequences to '{}'",
        consensus_path.display()
    );
    let mut consensus_seqs = number_fasta(&batch_consensus, &consensus_path)?;
    batch_consensus = consensus_path;

    // Start merging
    let trailing_semi = Regex::new(";$").unwrap();
    let mut prev_scores: Option<PathBuf> = None;
    let mut next_family_number = 1;
    let mut sculu_name_to_desc: HashMap<String, String> = HashMap::new();
    //let mut prev_msa_result: Option<PathBuf> = None;
    for round in 1.. {
        debug!(">>> Round {round} <<<");
        let round_dir = batch_dir.join(format!("round{round:03}"));
        fs::create_dir_all(&round_dir)?;

        // Align all the sequences to the current consensus sequences.
        // On the first round, all the original consensus will be included.
        // On future rounds, only the newly merged consensus will be present.
        let alignment_file = match args.alphabet {
            SequenceAlphabet::Dna => run_rmblastn(
                &round_dir,
                &batch_consensus,
                all_seqs_path,
                config,
                num_threads,
            )?,
            _ => run_blastp(
                &round_dir,
                &batch_consensus,
                all_seqs_path,
                config,
                num_threads,
            )?,
        };

        // Extract the scores from the alignment file.
        // On the first round, there will be no "prev_scores" file.
        // On the following rounds, the previous round's scores will be
        // included for those families that were not merged.
        let scores_file = extract_scores(
            &alignment_file,
            &prev_scores,
            &batch_consensus,
            &round_dir,
        )?;

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
        // The new consensus file will only contain the merged families.
        let new_consensus_path = &round_dir.join("new-consensus.fa");
        let mut new_consensus = open_for_write(new_consensus_path)?;
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
                    taken_instances_dir: &args.instances,
                    num_threads,
                    alphabet: args.alphabet.clone(),
                    flipped: flipped_pairs
                        .contains(&StringPair::new(pair.f1.clone(), pair.f2.clone())),
                },
                &mut family_to_instance,
            )?;

            let f1_desc = sculu_name_to_desc.get(&pair.f1).unwrap_or(&pair.f1);
            let f2_desc = sculu_name_to_desc.get(&pair.f2).unwrap_or(&pair.f2);
            let new_family_newick = format!(
                "({},{}):{:0.02};",
                //trailing_semi.replace(&pair.f1, ""),
                //trailing_semi.replace(&pair.f2, ""),
                trailing_semi.replace(&f1_desc, ""),
                trailing_semi.replace(&f2_desc, ""),
                pair.val
            );

            let f1_len = consensus_seqs.get(&pair.f1).map_or(0, |v| v.len());
            let f2_len = consensus_seqs.get(&pair.f2).map_or(0, |v| v.len());
            let new_len = new_consensus_seq.len();

            debug!(
                "{} len was {f1_len}, {} len was {f2_len}, \
                    new consensus len is {new_len}",
                pair.f1, pair.f2
            );

            // Update the consensus mapping
            consensus_seqs.remove(&pair.f1);
            consensus_seqs.remove(&pair.f2);
            consensus_seqs.insert(new_family_newick, new_consensus_seq.clone());

            // Keep track of all the families that have been merged
            for family in all_families {
                already_merged.insert(family);
            }
        }

        // Write the merged families to the new consensus file
        debug!(
            "Writing {} consensus sequences to '{}'",
            &consensus_seqs.len(),
            new_consensus_path.display()
        );

        //for (family_number, (family_name, seq)) in consensus_seqs.iter().enumerate() {
        for (family_name, seq) in consensus_seqs.iter() {
            let sculu_name = format!("sculufam{component_number}-{next_family_number}");
            let _ = sculu_name_to_desc
                .insert(sculu_name.to_string(), family_name.to_string());
            writeln!(new_consensus, ">{sculu_name} {family_name}\n{seq}",)?;
            next_family_number += 1;
        }

        batch_consensus = new_consensus_path.to_path_buf();
    }

    let outfile = batch_dir.join("final.fa");
    let mut final_writer = FastaWriter::new(BufWriter::new(open_for_write(&outfile)?));
    let mut reader = FastaReader::new(open(&batch_consensus)?);
    let mut new_seqs = 0;
    for result in reader.records() {
        //let mut record = result?;
        //if let Some(desc) = record.description() {
        //    record = FastaRecord::new(
        //        FastaDefinition::new(desc, None),
        //        record.sequence().clone(),
        //    )
        //}
        let record = result?;
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
fn cat_sequences(
    instances_dir: &PathBuf,
    families: &[String],
    outpath: &PathBuf,
) -> Result<usize> {
    let mut output = open_for_write(outpath)?;
    let entries = fs::read_dir(instances_dir)
        .map_err(|e| anyhow!(r#"Failed to read "{}": {e}"#, instances_dir.display()))?;
    let mut count = 0;

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
            count += 1;
        }
    }

    Ok(count)
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
        "Taking {num_wanted} from '{}'{}",
        fasta.display(),
        if rev_comp { " (RevComp)" } else { "" }
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
fn call_winners(
    scores_file: &PathBuf,
    lambda: f64,
    confidence_margin: isize,
) -> Result<Winners> {
    debug!(
        r#"Calling winners from scores file "{}""#,
        scores_file.display()
    );

    let mut reader = csv::ReaderBuilder::new()
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
            //let comps: Vec<_> = pairs
            //    .iter()
            //    .map(|(x, y)| (x * (1. / confidence_margin as f64)) > *y)
            //    .collect();
            //all_comps.push(comps.iter().all(|&v| v));
            let all_ok = pairs
                .iter()
                .all(|(x, y)| (x * (1. / confidence_margin as f64)) > *y);
            all_comps.push(all_ok);
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
fn extract_scores(
    alignment_file: &PathBuf,
    prev_scores_file: &Option<PathBuf>,
    consensus_path: &PathBuf,
    outdir: &Path,
) -> Result<PathBuf> {
    debug!(
        "Extracting scores from alignment '{}' {:?}",
        alignment_file.display(),
        prev_scores_file
    );

    let mut reader = parse_reader(open(consensus_path)?)?;
    let mut consensus_names: HashMap<String, String> = HashMap::new();
    while let Some(rec) = reader.iter_record()? {
        consensus_names
            .insert(rec.head().trim().to_string(), rec.des().trim().to_string());
    }
    debug!("consensus_names = {consensus_names:#?}");

    let scores_file = outdir.join("alignment-scores.tsv");
    let mut scores_wtr = csv::WriterBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(&scores_file)
        .map_err(|e| anyhow!("Failed to write '{}': {e}", &scores_file.display()))?;

    // This is a hash to remember the family names we handled
    // in this round. When reading the scores from the previous
    // round, we need to skip these queries.
    let mut skip_query: HashSet<String> = HashSet::new();

    // BLAST output fails to include headers
    let mut align_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_path(alignment_file)
        .map_err(|e| anyhow!("Failed to read '{}': {e}", &alignment_file.display()))?;

    let records = align_reader.records();
    for res in records {
        let res = res?;
        let rec: RmBlastOutput = res.deserialize(None)?;
        match consensus_names.get(&rec.query) {
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
        let mut reader = csv::ReaderBuilder::new()
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
fn number_fasta(
    in_path: &PathBuf,
    out_path: &PathBuf,
) -> Result<HashMap<String, String>> {
    let mut outfile = open_for_write(out_path)?;
    let mut reader = parse_reader(open(in_path)?)?;
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
fn find_independence(num_wins: u32, num_shared: u32) -> f64 {
    if num_shared > 0 {
        num_wins as f64 / (num_wins as f64 + num_shared as f64)
    } else {
        1.
    }
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
