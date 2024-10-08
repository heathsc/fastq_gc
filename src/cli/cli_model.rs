use std::path::PathBuf;

use clap::{command, value_parser, Arg, ArgAction, Command};

use crate::utils::{BisulfiteType, LogLevel};

pub(super) fn cli_model() -> Command {
    command!()
        .arg(
            Arg::new("timestamp")
                .short('X')
                .long("timestamp")
                .value_parser(value_parser!(stderrlog::Timestamp))
                .value_name("GRANULARITY")
                .default_value("none")
                .help("Prepend log entries with a timestamp"),
        )
        .arg(
            Arg::new("loglevel")
                .short('l')
                .long("loglevel")
                .value_name("LOGLEVEL")
                .value_parser(value_parser!(LogLevel))
                .ignore_case(true)
                .default_value("info")
                .help("Set log level"),
        )
        .arg(
            Arg::new("quiet")
                .action(ArgAction::SetTrue)
                .long("quiet")
                .conflicts_with("loglevel")
                .help("Silence all output"),
        )
        .arg(
            Arg::new("threads")
                .short('t')
                .long("threads")
                .value_parser(value_parser!(u64).range(1..))
                .value_name("INT")
                .default_value("1")
                .help("Set number of process threads"),
        )
        .arg(
            Arg::new("trim")
                .long("trim")
                .value_parser(value_parser!(usize))
                .value_name("BASES")
                .default_value("0")
                .help("Trim bases from start of read"),
        )
        .arg(
            Arg::new("min_qual")
                .long("min-qual")
                .value_parser(value_parser!(u8))
                .value_name("QUAL")
                .default_value("0")
                .help("Minimum base quality to consider"),
        )
        .arg(
            Arg::new("bisulfite")
                .action(ArgAction::SetTrue)
                .long("bisulfite")
                .short('b')
                .conflicts_with("kmers")
                .help("Set bisulfite mode to Forward"),
        )
        .arg(
            Arg::new("bisulfite_type")
                .short('B')
                .long("bisulfite-type")
                .value_name("TYPE")
                .value_parser(value_parser!(BisulfiteType))
                .ignore_case(true)
                .conflicts_with("bisulfite")
                .conflicts_with("kmers")
                .help("Set specific bisulfite type"),
        )
        .arg(
            Arg::new("flowcell")
                .long("flowcell")
                .value_parser(value_parser!(String))
                .value_name("ID")
                .help("Flowcell ID"),
        )
        .arg(
            Arg::new("library")
                .long("library")
                .value_parser(value_parser!(String))
                .value_name("ID")
                .help("Library ID"),
        )
        .arg(
            Arg::new("sample")
                .long("sample")
                .value_parser(value_parser!(String))
                .value_name("ID")
                .help("Sample ID"),
        )
        .arg(
            Arg::new("barcode")
                .long("barcode")
                .value_parser(value_parser!(String))
                .value_name("BARCODE")
                .help("Sample barcode"),
        )
        .arg(
            Arg::new("lane")
                .long("lane")
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .help("Lane number"),
        )
        .arg(
            Arg::new("index")
                .long("index")
                .value_parser(value_parser!(String))
                .value_name("ID")
                .help("Index ID"),
        )
        .arg(
            Arg::new("read_end")
                .long("read-end")
                .value_parser(value_parser!(u8))
                .value_name("INT")
                .help("Read end"),
        )
        .arg(
            Arg::new("control")
                .long("control-fasta")
                .short('c')
                .value_parser(value_parser!(PathBuf))
                .value_name("FASTA")
                .help("Input FASTA file with control sequences to be excluded"),
        )
        .arg(
            Arg::new("kmers")
                .long("kmers")
                .short('k')
                .value_parser(value_parser!(PathBuf))
                .value_name("KM FILE")
                .help("Input KM file with kmers for coverage estimation"),
        )
        .arg(
            Arg::new("output_prefix")
                .long("output_prefix")
                .short('o')
                .value_parser(value_parser!(String))
                .value_name("OUTPUT PREFIX")
                .default_value("fastq_gc")
                .help("Prefix for output JSON file"),
        )
        .arg(
            Arg::new("input")
                .value_parser(value_parser!(PathBuf))
                .value_name("INPUT")
                .action(ArgAction::Append)
                .help("Input FASTQ file"),
        )
}
