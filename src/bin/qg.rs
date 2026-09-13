//! `qg`: enumeration, signatures and gate checks for quadratula.
//!
//! Subcommands (all paths explicit, all output plain text):
//!   laws                          survey of the law file (operations, variables, shapes)
//!   enumerate <n> <out.tsv.gz>    iso classes of order n with automorphism counts and signatures
//!   count-latin <n>               all Latin squares of order n, by plain backtracking
//!   parastrophe <out.tsv>         induced S3 action on the 990 laws
//!   g4 <maxn>                     sig(Q^s) = sig(Q) o tau_s for every labeled Latin square of order <= maxn
//!   controls <floor.tsv.gz>...    classical counts (groups, commutative, ...) over enumerated floors

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use quadratula::enumerate::{all_latin, count_latin, factorial, iso_classes};
use quadratula::law::{Law, Op, parse_law_body};
use quadratula::parastrophe::{S3, law_map, name};
use quadratula::qg::{Compiled, Engine, Quasigroup};
use quadratula::{Signature, load_laws};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let cmd = args.get(1).map(String::as_str).unwrap_or("");
    let arg = |i: usize| -> &str {
        args.get(i).map(String::as_str).unwrap_or_else(|| {
            eprintln!("missing argument {i} for {cmd}");
            std::process::exit(2)
        })
    };
    match cmd {
        "laws" => laws(),
        "enumerate" => enumerate(arg(2).parse().unwrap(), arg(3)),
        "count-latin" => {
            let n: usize = arg(2).parse().unwrap();
            let t = Instant::now();
            let c = count_latin(n);
            println!("order\t{n}\tlatin_squares\t{c}\tseconds\t{:.2}", t.elapsed().as_secs_f64());
        }
        "parastrophe" => parastrophe(arg(2)),
        "g4" => g4(arg(2).parse().unwrap()),
        "controls" => controls(&args[2..]),
        "check" => check(arg(2)),
        "models" => models(arg(2).parse().unwrap(), arg(3), arg(4)),
        other => {
            eprintln!("unknown command {other:?}; try: laws, enumerate, count-latin, parastrophe, g4, controls, check");
            std::process::exit(2);
        }
    }
}

fn laws() {
    let laws = load_laws();
    let mut op_count: BTreeMap<&str, usize> = BTreeMap::new();
    let mut vars_hist: BTreeMap<usize, usize> = BTreeMap::new();
    let mut ops_hist: BTreeMap<(usize, usize), usize> = BTreeMap::new();
    let mut letters: BTreeMap<char, usize> = BTreeMap::new();
    let mut max_depth = 0;
    for l in &laws {
        for side in [&l.lhs, &l.rhs] {
            side.visit_ops(&mut |o| *op_count.entry(o.token()).or_default() += 1);
            side.visit_vars(&mut |v| *letters.entry((b'a' + v) as char).or_default() += 1);
            max_depth = max_depth.max(side.depth());
        }
        *vars_hist.entry(l.vars().len()).or_default() += 1;
        *ops_hist.entry((l.lhs.ops(), l.rhs.ops())).or_default() += 1;
    }
    let keys: std::collections::HashSet<String> = laws.iter().map(Law::key).collect();
    let lhs_xyz = laws.iter().filter(|l| l.lhs.to_string().replace(['(', ')', ' ', '*', '/', '\\'], "") == "xyz").count();
    println!("laws\t{}", laws.len());
    for (k, v) in &op_count {
        println!("op_occurrences\t{k}\t{v}");
    }
    for (k, v) in &letters {
        println!("variable_occurrences\t{k}\t{v}");
    }
    for (k, v) in &vars_hist {
        println!("laws_with_distinct_variables\t{k}\t{v}");
    }
    for ((a, b), v) in &ops_hist {
        println!("laws_with_ops_lhs_rhs\t{a}\t{b}\t{v}");
    }
    println!("max_side_depth\t{max_depth}");
    println!("distinct_keys_up_to_renaming_and_side_swap\t{}", keys.len());
    println!("tautologies\t{}", laws.iter().filter(|l| l.is_tautology()).count());
    println!("lhs_reads_x_y_z_in_order\t{lhs_xyz}");
    // Round trip: printing reproduces the file text.
    let text = std::fs::read_to_string(quadratula::LAWS_PATH).unwrap();
    let rt = text.lines().zip(&laws).filter(|(line, l)| line.split_once(": ").unwrap().1 == l.to_string()).count();
    println!("round_trip_identical_lines\t{rt}");
}

fn sat_string(sig: &Signature) -> String {
    let ids = sig.satisfied_ids();
    if ids.is_empty() {
        "-".into()
    } else {
        ids.iter().map(u32::to_string).collect::<Vec<_>>().join(",")
    }
}

fn enumerate(n: usize, out: &str) {
    let laws = load_laws();
    let engine = Engine::new(&laws);
    let t0 = Instant::now();
    let classes = iso_classes(n);
    let t_enum = t0.elapsed().as_secs_f64();
    let t1 = Instant::now();
    let sigs: Vec<Signature> = classes
        .par_iter()
        .map(|c| engine.signature(&Quasigroup::from_table(n, &c.table).unwrap()))
        .collect();
    let t_sig = t1.elapsed().as_secs_f64();
    let orbit_sum: u64 = classes.iter().map(|c| factorial(n) / c.aut as u64).sum();
    let f = std::fs::File::create(out).expect("create output");
    let mut w = BufWriter::new(GzEncoder::new(f, Compression::default()));
    writeln!(w, "order\ttable\taut\tsatisfied").unwrap();
    for (c, s) in classes.iter().zip(&sigs) {
        let q = Quasigroup::from_table(n, &c.table).unwrap();
        writeln!(w, "{n}\t{}\t{}\t{}", q.table_string(), c.aut, sat_string(s)).unwrap();
    }
    w.into_inner().unwrap().finish().unwrap();
    println!(
        "order\t{n}\tiso_classes\t{}\tsum_n!/aut\t{orbit_sum}\tenumerate_seconds\t{t_enum:.2}\tsignature_seconds\t{t_sig:.2}\tthreads\t{}",
        classes.len(),
        rayon::current_num_threads()
    );
}

fn parastrophe(out: &str) {
    let laws = load_laws();
    let mut w = BufWriter::new(std::fs::File::create(out).unwrap());
    writeln!(w, "sigma\tlaws_mapped_into_list\tis_bijection\tmap").unwrap();
    for s in S3 {
        let m = law_map(&laws, s);
        let mapped = m.iter().filter(|x| x.is_some()).count();
        let mut seen = vec![false; laws.len()];
        m.iter().flatten().for_each(|&j| seen[j] = true);
        let bij = mapped == laws.len() && seen.iter().all(|&b| b);
        let ids: Vec<String> = m.iter().map(|x| x.map_or(0, |j| j + 1).to_string()).collect();
        writeln!(w, "{}\t{mapped}\t{bij}\t{}", name(s), ids.join(",")).unwrap();
        let ops: Vec<String> = Op::ALL
            .iter()
            .map(|&o| {
                let (o2, swap) = quadratula::parastrophe::translate_op(s, o);
                format!("{}->{}{}", o.token(), o2.token(), if swap { "(swapped)" } else { "" })
            })
            .collect();
        println!("sigma\t{}\tmapped\t{mapped}\tbijection\t{bij}\tops\t{}", name(s), ops.join(" "));
    }
}

fn g4(maxn: usize) {
    let laws = load_laws();
    let engine = Engine::new(&laws);
    let maps: Vec<Vec<Option<usize>>> = S3.iter().map(|&s| law_map(&laws, s)).collect();
    let mut failures = 0u64;
    for n in 1..=maxn {
        let squares = all_latin(n);
        let (checked, bad): (u64, u64) = squares
            .par_iter()
            .map(|t| {
                let q = Quasigroup::from_table(n, t).unwrap();
                let base = engine.signature(&q);
                let mut bad = 0u64;
                let mut checked = 0u64;
                for (si, &s) in S3.iter().enumerate() {
                    let conj = engine.signature(&q.conjugate(s));
                    for (l, tl) in maps[si].iter().enumerate() {
                        if let Some(tl) = tl {
                            checked += 1;
                            if conj.get(l) != base.get(*tl) {
                                bad += 1;
                            }
                        }
                    }
                }
                (checked, bad)
            })
            .reduce(|| (0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
        failures += bad;
        println!("order\t{n}\tlabeled_squares\t{}\tbit_checks\t{checked}\tmismatches\t{bad}", squares.len());
    }
    if failures > 0 {
        std::process::exit(1);
    }
}

/// All quasigroups of order `n` (up to isomorphism) satisfying the listed
/// Schröder laws (`ids` comma-separated, or `-` for none). Same output format
/// as `enumerate`.
fn models(n: usize, ids: &str, out: &str) {
    let laws = load_laws();
    let engine = Engine::new(&laws);
    let required: Vec<Law> = if ids == "-" {
        vec![]
    } else {
        ids.split(',').map(|x| laws[x.parse::<usize>().unwrap() - 1].clone()).collect()
    };
    let t0 = Instant::now();
    let m = quadratula::search::models(n, &required, &engine);
    let secs = t0.elapsed().as_secs_f64();
    let f = std::fs::File::create(out).expect("create output");
    let mut w = BufWriter::new(GzEncoder::new(f, Compression::default()));
    writeln!(w, "order\ttable\taut\tsatisfied").unwrap();
    for model in &m.models {
        let q = Quasigroup::from_table(n, &model.table).unwrap();
        writeln!(w, "{n}\t{}\t{}\t{}", q.table_string(), model.aut, sat_string(&model.signature)).unwrap();
    }
    w.into_inner().unwrap().finish().unwrap();
    println!(
        "order\t{n}\tlaws\t{ids}\tiso_models\t{}\tnodes\t{}\ttasks\t{}\tseconds\t{secs:.2}",
        m.models.len(),
        m.nodes,
        m.tasks
    );
}

/// Signatures of externally supplied tables. Input lines: `key<TAB>n<TAB>entries`,
/// entries comma-separated row-major. Output lines: `key<TAB>satisfied`.
fn check(path: &str) {
    let laws = load_laws();
    let engine = Engine::new(&laws);
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    for line in text.lines().filter(|l| !l.trim().is_empty()) {
        let cols: Vec<&str> = line.split('\t').collect();
        let n: usize = cols[1].parse().unwrap();
        let table: Vec<u8> = cols[2].split(',').map(|x| x.trim().parse().unwrap()).collect();
        match Quasigroup::from_table(n, &table) {
            Ok(q) => println!("{}\t{}", cols[0], sat_string(&engine.signature(&q))),
            Err(e) => println!("{}\tNOT_LATIN:{e}", cols[0]),
        }
    }
}

pub fn read_floor(path: &str) -> Vec<(usize, Vec<u8>, u32, Vec<u32>)> {
    let f = std::fs::File::open(path).unwrap_or_else(|e| panic!("open {path}: {e}"));
    let mut out = Vec::new();
    for line in BufReader::new(GzDecoder::new(f)).lines().skip(1) {
        let line = line.unwrap();
        let cols: Vec<&str> = line.split('\t').collect();
        let n: usize = cols[0].parse().unwrap();
        let table = cols[1].chars().map(|c| c.to_digit(36).unwrap() as u8).collect();
        let aut = cols[2].parse().unwrap();
        let sat = if cols[3] == "-" { vec![] } else { cols[3].split(',').map(|x| x.parse().unwrap()).collect() };
        out.push((n, table, aut, sat));
    }
    out
}

/// Classical counts, computed two ways: with the control law written directly
/// and, where Le Floch names a Schröder equivalent, from the stored signature.
fn controls(paths: &[String]) {
    let laws = load_laws();
    let index: std::collections::HashMap<String, u32> = laws.iter().map(|l| (l.key(), l.id)).collect();
    let direct = |s: &str| {
        let (l, r) = parse_law_body(s).unwrap();
        Compiled::new(&Law { id: 0, lhs: l, rhs: r })
    };
    let assoc = direct("x * (y * z) = (x * y) * z");
    let comm = direct("x * y = y * x");
    let semi = direct("x = y * (x * y)");
    let idem = direct("x * x = x");
    let key_of = |s: &str| {
        let (l, r) = parse_law_body(s).unwrap();
        Law { id: 0, lhs: l, rhs: r }.key()
    };
    let assoc_id = index.get(&key_of("x * (y * z) = (x * y) * z")).copied();
    println!("associativity_schroeder_id\t{}", assoc_id.map_or("absent".into(), |i| i.to_string()));
    println!("order\tall\tgroups\tabelian_groups\tcommutative\tsemisymmetric\ttotally_symmetric\tidempotent\tsch_assoc\tsch1_commutative\tsch5_semisymmetric");
    for p in paths {
        let floor = read_floor(p);
        let n = floor[0].0;
        let mut c = [0u64; 10];
        for (_, t, _, sat) in &floor {
            let q = Quasigroup::from_table(n, t).unwrap();
            let (a, m, s, i) = (assoc.holds(&q), comm.holds(&q), semi.holds(&q), idem.holds(&q));
            c[0] += 1;
            c[1] += a as u64;
            c[2] += (a && m) as u64;
            c[3] += m as u64;
            c[4] += s as u64;
            c[5] += (s && m) as u64;
            c[6] += i as u64;
            c[7] += assoc_id.is_some_and(|id| sat.contains(&id)) as u64;
            c[8] += sat.contains(&1) as u64;
            c[9] += sat.contains(&5) as u64;
        }
        let cols: Vec<String> = c.iter().map(u64::to_string).collect();
        println!("{n}\t{}", cols.join("\t"));
    }
}
