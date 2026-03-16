use crate::{
    common::{copy_fasta, open, open_for_write, read_lines},
    types::ConcatArgs,
};
use anyhow::{bail, Result};
use log::debug;
use noodles_fasta::{
    self,
    io::Reader as FastaReader,
    io::Writer as FastaWriter,
    record::{Definition as FastaDefinition, Record as FastaRecord},
};
use std::{
    fs,
    io::{BufReader, BufWriter},
};

// --------------------------------------------------
pub fn concat_files(args: &ConcatArgs) -> Result<()> {
    let intermediate_dir = args.outdir.join("intermediate");
    let final_seed_alignments_dir = args.outdir.join("seed_alignments");
    fs::create_dir_all(&final_seed_alignments_dir)?;

    let outfile = args.outdir.join("sculu_families.fa");
    let mut fasta_writer = FastaWriter::new(BufWriter::new(open_for_write(&outfile)?));

    if let Some(singletons_file) = &args.singletons {
        // The singletons file will have one family name per line
        let family_names = read_lines(singletons_file)?;
        debug!(
            "Copying {} sequences from singletons file",
            family_names.len()
        );

        // Copy singleton's original consensus sequences to the final output file
        copy_fasta(&family_names, &args.consensus_path, &mut fasta_writer)?;

        // Also copy the original seed alignments from the "build_components" step
        let seed_alignments_dir = intermediate_dir.join("seed_alignments");
        if !seed_alignments_dir.is_dir() {
            bail!(r#"Missing "{}""#, seed_alignments_dir.display());
        }

        for family_name in family_names {
            let filename = format!("{family_name}.stk");
            let seed_alignments = seed_alignments_dir.join(&filename);
            if seed_alignments.is_file() {
                fs::copy(seed_alignments, final_seed_alignments_dir.join(&filename))?;
            } else {
                eprintln!(r#"Missing "{}""#, seed_alignments.display());
            }
        }
    }

    for component_file in &args.components {
        let component_dir = component_file.parent().expect("Cannot get parent dir");
        debug!(r#"Copying from "{}""#, component_file.display());
        let mut reader = FastaReader::new(BufReader::new(open(component_file)?));
        for (record_num, result) in reader.records().enumerate() {
            let mut record = result?;
            let cur_family_name = String::from_utf8(record.name().to_vec())?;
            let new_family_name = format!("sculufam-{}", record_num + 1);
            record = FastaRecord::new(
                FastaDefinition::new(
                    new_family_name.clone(),
                    record.description().map(Into::into),
                ),
                record.sequence().clone(),
            );
            fasta_writer.write_record(&record)?;

            let msa = component_dir.join(format!("{cur_family_name}.stk"));
            let seed_alignments = component_dir.join(&msa);
            if seed_alignments.is_file() {
                let new_filename = format!("{new_family_name}.stk");
                let mut dest_fh =
                    open_for_write(&final_seed_alignments_dir.join(&new_filename))?;
                for line in read_lines(&msa)? {
                    writeln!(
                        dest_fh,
                        "{}",
                        if line.starts_with("#=GF ID ") {
                            format!("#GF ID {new_family_name}")
                        } else {
                            line
                        }
                    )?;
                }
            } else {
                eprintln!(r#"Missing "{}""#, seed_alignments.display());
            }
        }
    }

    debug!(r#"Final output written to "{}""#, outfile.display());

    Ok(())
}
