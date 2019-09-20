use census_proteomics::*;
use clap::{App, Arg, ArgGroup};
use serde_json;
use std::collections::HashMap;
use std::fs;
use std::io::prelude::*;
use std::path::{Path, PathBuf};

fn combine_protein<'a, P: AsRef<Path>>(
    path: P,
    filters: &Filter<'a>,
    average: bool,
) -> std::io::Result<()> {
    let mut outpath = PathBuf::from(path.as_ref());
    if !outpath.set_extension("csv") {
        panic!("Cannot set file extension for {}", outpath.display());
    }

    let file = fs::read_to_string(path)?;
    let data = census_proteomics::read_census(&file)
        .map_err(|_| std::io::Error::from(std::io::ErrorKind::InvalidData))?;
    let data = data.filter(filters);

    let mut file = fs::File::create(outpath)?;
    writeln!(
        file,
        "accession,description,spectral_count,sequence_count,{}",
        (1..=data.channels)
            .map(|i| format!("channel_{}", i))
            .collect::<Vec<String>>()
            .join(",")
    )?;

    for prot in &data.proteins {
        let adj = prot
            .total()
            .into_iter()
            .map(|v| {
                format!(
                    "{}",
                    if average {
                        v / prot.peptides.len() as u32
                    } else {
                        v
                    }
                )
            })
            .collect::<Vec<String>>()
            .join(",");

        writeln!(
            file,
            "{},{},{},{},{}",
            prot.accession,
            prot.description.replace(",", ";"),
            prot.spectral_count,
            prot.sequence_count,
            adj
        )?;
    }

    Ok(())
}

fn combine_peptide<'a, P: AsRef<Path>>(
    path: P,
    filters: &Filter<'a>,
    average: bool,
) -> std::io::Result<()> {
    let mut outpath = PathBuf::from(path.as_ref());
    if !outpath.set_extension("csv") {
        panic!("Cannot set file extension for {}", outpath.display());
    }

    let file = fs::read_to_string(path)?;
    let data = census_proteomics::read_census(&file).expect("Error parsing census file!");
    let data = data.filter(filters);

    let mut file = fs::File::create(outpath)?;
    writeln!(
        file,
        "accession,description,spectral_count,sequence_count,sequence,{}",
        (1..=data.channels)
            .map(|i| format!("channel_{}", i))
            .collect::<Vec<String>>()
            .join(",")
    )?;

    for prot in &data.proteins {
        let mut map: HashMap<&str, Vec<u32>> = HashMap::new();
        let mut cnt: HashMap<&str, u32> = HashMap::new();
        for peptide in &prot.peptides {
            let entry = map
                .entry(peptide.sequence)
                .or_insert((0..data.channels).map(|_| 0).collect::<Vec<u32>>());
            for (idx, val) in peptide.values.iter().enumerate() {
                entry[idx] += *val;
            }
            *cnt.entry(peptide.sequence).or_insert(0) += 1;
        }

        for (sequence, summed_values) in map {
            let adj = summed_values
                .into_iter()
                .map(|v| format!("{}", if average { v / cnt[sequence] } else { v }))
                .collect::<Vec<String>>()
                .join(",");

            writeln!(
                file,
                "{},{},{},{},{},{}",
                prot.accession,
                prot.description.replace(",", ";"),
                prot.spectral_count,
                prot.sequence_count,
                sequence,
                adj
            )?;
        }
    }

    Ok(())
}

fn generate_example() {
    let filter = Filter::default()
        .add_peptide_filter(PeptideFilter::ChannelIntensity(1, 1000))
        .add_peptide_filter(PeptideFilter::ChannelCV(vec![1, 2], 0.6))
        .add_peptide_filter(PeptideFilter::Unique)
        .add_peptide_filter(PeptideFilter::Tryptic)
        .add_peptide_filter(PeptideFilter::TotalIntensity(5000))
        .add_protein_filter(ProteinFilter::ExcludeReverse)
        .add_protein_filter(ProteinFilter::SequenceCounts(2));
    let s = serde_json::to_string_pretty(&filter).expect("Serialization error");
    let mut f = fs::File::create("filter.json")
        .expect("Could not open filter.json for writing")
        .write_all(s.as_bytes())
        .expect("Error writing to filter.json");
}

fn main() {
    let matches = App::new("census2csv")
        .version("0.1")
        .author("Michael R. Lazar <lazear@scripps.edu>")
        .about("Parse, filter, and convert census out files to csv")
        .group(
            ArgGroup::with_name("combine")
                .required(true)
                .args(&["peptide", "protein"]),
        )
        .arg(
            Arg::with_name("peptide")
                .help("Output peptide-level data")
                .long("peptide")
                .short("e")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("protein")
                .help("Output protein-level data")
                .long("protein")
                .short("r")
                .takes_value(false),
        )
        .arg(
            Arg::with_name("filter")
                .help("JSON file containing filters to apply")
                .short("f")
                .long("filter")
                .value_name("FILE")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("average")
                .help("Average results by number of reported spectral matches, default is sum")
                .short("a")
                .long("avg"),
        )
        .arg(
            Arg::with_name("INPUT")
                .help("list of input files to convert")
                .multiple(true),
        )
        .get_matches();

    // declare up here to get around borrowck
    let mut filterbuf: String;

    let filter = match matches.value_of("filter") {
        Some(path) => {
            filterbuf = fs::read_to_string(path).expect(&format!("Cannot read {}", path));
            match serde_json::from_str(&filterbuf) {
                Ok(f) => f,
                Err(e) => {
                    println!("Error while parsing filter.json {:?}", e);
                    std::process::abort();
                }
            }
        }
        None => Filter::default(),
    };

    for f in matches
        .value_of("INPUT")
        .expect("No input files!")
        .split_whitespace()
    {
        let res = if matches.is_present("peptide") {
            combine_peptide(f, &filter, matches.is_present("average"))
        } else {
            combine_protein(f, &filter, matches.is_present("average"))
        };
        if let Err(e) = res {
            println!("Error during processing of file {}: {}", f, e);
        }
    }
}
