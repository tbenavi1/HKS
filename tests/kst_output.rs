// `hks lookup --format kst` against `--format bed`: the .kst must decode to the BED,
// and be byte-identical to `karyozip -e <effort> compress` of that BED.
//
// The byte comparison needs the karyozip CLI: $KARYOZIP, else `karyozip` on PATH. Without
// it those comparisons are skipped with a note; the decode comparison always runs, through
// the karyozip library hks links.

use std::{path::{Path, PathBuf}, process::{Command, Output}, sync::OnceLock};

const HKS: &str = env!("CARGO_BIN_EXE_hks");
const K: usize = 15;
const N_READS: usize = 10000;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 { // xorshift64*
        self.0 ^= self.0 >> 12; self.0 ^= self.0 << 25; self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: usize) -> usize { (self.next() % n as u64) as usize }
    fn base(&mut self) -> u8 { b"ACGT"[self.below(4)] }
    fn seq(&mut self, n: usize) -> Vec<u8> { (0..n).map(|_| self.base()).collect() }
}

fn revcomp(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|&c| match c { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', x => x }).collect()
}

fn write_fasta(path: &Path, recs: &[(String, Vec<u8>)]) {
    let mut s = Vec::new();
    for (name, seq) in recs {
        s.extend_from_slice(format!(">{name}\n").as_bytes());
        s.extend_from_slice(seq);
        s.push(b'\n');
    }
    std::fs::write(path, s).unwrap();
}

fn run(cmd: &mut Command) -> Output {
    let out = cmd.output().unwrap_or_else(|e| panic!("could not run {cmd:?}: {e}"));
    assert!(out.status.success(), "{cmd:?} failed:\n{}", String::from_utf8_lossy(&out.stderr));
    out
}

/// An index of four related genomes and two query files, built once for all tests.
fn fixture() -> &'static Path {
    static DIR: OnceLock<PathBuf> = OnceLock::new();
    DIR.get_or_init(|| {
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR")).join("kst_fixture");
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let mut rng = Rng(0x9E3779B97F4A7C15);

        // B shares half of A and C carries a piece of D, so runs get internal labels too
        let a = rng.seq(20000);
        let mut b = a[..10000].to_vec(); b.extend(rng.seq(10000));
        let d = rng.seq(20000);
        let mut c = rng.seq(15000); c.extend_from_slice(&d[2000..7000]);
        let genomes = vec![("A".to_string(), a), ("B".to_string(), b), ("C".to_string(), c), ("D".to_string(), d)];
        write_fasta(&dir.join("genomes.fa"), &genomes);
        std::fs::write(dir.join("hierarchy.txt"), "A clade1\nB clade1\nclade1 clade2\nC clade2\nclade2 root\nD root\n").unwrap();

        // Reads from either strand with 4% substitutions; some unrelated, some shorter than k
        let mut reads = |n: usize, prefix: &str| -> Vec<(String, Vec<u8>)> {
            (0..n).map(|i| {
                let len = if rng.below(50) == 0 { 1 + rng.below(K) } else { 50 + rng.below(1450) };
                let mut s = if rng.below(20) == 0 { rng.seq(len) } else {
                    let g = &genomes[rng.below(4)].1;
                    let p = rng.below(g.len() - len);
                    g[p..p + len].to_vec()
                };
                if rng.below(2) == 0 { s = revcomp(&s); }
                for x in s.iter_mut() { if rng.below(25) == 0 { *x = rng.base(); } }
                (format!("{prefix}_{i} some description"), s)
            }).collect()
        };
        let q1 = reads(N_READS, "read");
        let q2 = reads(500, "other");
        write_fasta(&dir.join("q1.fa"), &q1);
        write_fasta(&dir.join("q2.fa"), &q2);

        run(Command::new(HKS).current_dir(&dir).args(["build-base", "-s", &K.to_string(), "--input", "genomes.fa", "-o", "base.hksb", "-t", "2"]));
        run(Command::new(HKS).current_dir(&dir).args(["add-feature-set", "-i", "base.hksb", "-o", "f.hksf",
            "--feature-per-seq-file", "genomes.fa", "--feature-hierarchy", "hierarchy.txt", "--feature-set-name", "t", "-t", "2"]));
        dir
    })
}

/// `hks lookup` on the fixture index, with KaryoScope's flags plus `extra`.
fn lookup(extra: &[&str]) -> Output {
    let dir = fixture();
    let mut cmd = Command::new(HKS);
    cmd.current_dir(dir).args(["lookup", "-i", "base.hksb", "--feature-set-file", "f.hksf", "--no-header"]).args(extra);
    cmd.output().unwrap()
}

fn lookup_ok(extra: &[&str]) {
    let out = lookup(extra);
    assert!(out.status.success(), "hks lookup {extra:?} failed:\n{}", String::from_utf8_lossy(&out.stderr));
}

fn karyozip_cli() -> Option<String> {
    let bin = std::env::var("KARYOZIP").unwrap_or_else(|_| "karyozip".to_string());
    match Command::new(&bin).arg("--version").output() {
        Ok(o) if o.status.success() => Some(bin),
        _ => {
            eprintln!("note: no karyozip CLI ($KARYOZIP or on PATH): skipping the byte comparison with karyozip compress");
            None
        }
    }
}

fn decode(kst: &Path) -> Vec<u8> {
    let mut r = karyozip::reader::Reader::open_with(kst.to_str().unwrap(), &Default::default());
    let mut out = Vec::new();
    r.decode(&mut out);
    out
}

fn read(name: &str) -> Vec<u8> {
    std::fs::read(fixture().join(name)).unwrap_or_else(|e| panic!("{name}: {e}"))
}

/// The .kst decodes to the BED, and equals karyozip's compress of that BED.
fn check_against_bed(bed: &str, kst: &str, effort: &str) {
    let dir = fixture();
    let bed_bytes = read(bed);
    assert!(!bed_bytes.is_empty());
    assert!(decode(&dir.join(kst)) == bed_bytes, "{kst} does not decode to {bed}");
    if let Some(kz) = karyozip_cli() {
        let reference = format!("{kst}.karyozip");
        run(Command::new(kz).current_dir(dir).args(["-e", effort, "compress", "-i", bed, "-o", &reference]));
        assert!(read(kst) == read(&reference), "{kst} differs from karyozip -e {effort} compress -i {bed}");
    }
}

fn n_lines(name: &str) -> usize { read(name).iter().filter(|&&c| c == b'\n').count() }

#[test]
fn kst_matches_bed_ranks_names_and_label_ids() {
    let cases: [(&str, &[&str]); 3] = [
        ("ranks", &[]),
        ("names", &["--report-query-names"]),
        ("ids", &["--report-label-ids"]),
    ];
    for (tag, flags) in cases {
        let bed = format!("q1.{tag}.bed");
        let mut args = vec!["--report-misses", "--miss-label", "novel", "-t", "4", "-q", "q1.fa"];
        args.extend_from_slice(flags);
        let mut with_out = args.clone();
        with_out.extend(["-o", &bed]);
        lookup_ok(&with_out);
        // More runs than a -e low block holds, so the track has several blocks
        assert!(n_lines(&bed) > 262144, "{bed}: only {} runs", n_lines(&bed));
        for effort in ["low", "medium"] {
            let kst = format!("q1.{tag}.{effort}.kst");
            let mut with_out = args.clone();
            with_out.extend(["--format", "kst", "--kst-effort", effort, "-o", &kst]);
            lookup_ok(&with_out);
            check_against_bed(&bed, &kst, effort);
        }
    }
}

#[test]
fn kst_default_effort_is_low_and_default_miss_label() {
    let common = ["--report-misses", "-t", "4", "-q", "q2.fa"];
    lookup_ok(&[&common[..], &["-o", "q2.default.bed"]].concat());
    lookup_ok(&[&common[..], &["--format", "kst", "-o", "q2.default.kst"]].concat());
    assert!(read("q2.default.bed").windows(5).any(|w| w == b"\tnone"), "expected misses labelled none");
    check_against_bed("q2.default.bed", "q2.default.kst", "low");
}

#[test]
fn kst_batch_matches_single_runs() {
    let common = ["--report-misses", "--miss-label", "novel", "--report-query-names", "-t", "4"];
    lookup_ok(&[&common[..], &["-q", "q1.fa", "-o", "b1.bed", "-q", "q2.fa", "-o", "b2.bed"]].concat());
    lookup_ok(&[&common[..], &["--format", "kst", "-q", "q1.fa", "-o", "b1.kst", "-q", "q2.fa", "-o", "b2.kst"]].concat());
    lookup_ok(&[&common[..], &["--format", "kst", "-q", "q2.fa", "-o", "s2.kst"]].concat());
    check_against_bed("b1.bed", "b1.kst", "low");
    check_against_bed("b2.bed", "b2.kst", "low");
    assert!(read("b2.kst") == read("s2.kst"), "batch and single-query output differ");
}

#[test]
fn kst_does_not_depend_on_thread_count() {
    let common = ["--report-misses", "--miss-label", "novel", "--format", "kst", "-q", "q1.fa"];
    lookup_ok(&[&common[..], &["-t", "1", "-o", "t1.kst"]].concat());
    lookup_ok(&[&common[..], &["-t", "8", "-o", "t8.kst"]].concat());
    assert!(read("t1.kst") == read("t8.kst"));
}

#[test]
fn kst_index_matches_karyozip_index() {
    lookup_ok(&["--report-misses", "--report-query-names", "--format", "kst", "--kst-index", "-q", "q2.fa", "-o", "ix.kst"]);
    assert!(fixture().join("ix.kst.kzi").exists());
    if let Some(kz) = karyozip_cli() {
        std::fs::copy(fixture().join("ix.kst"), fixture().join("ix2.kst")).unwrap();
        run(Command::new(kz).current_dir(fixture()).args(["index", "-i", "ix2.kst"]));
        assert!(read("ix.kst.kzi") == read("ix2.kst.kzi"));
    }
}

#[test]
fn kst_refuses_without_report_misses() {
    let out = lookup(&["--format", "kst", "-q", "q2.fa", "-o", "refused.kst"]);
    assert!(!out.status.success());
    let err = String::from_utf8_lossy(&out.stderr);
    assert!(err.contains("--format kst needs --report-misses"), "unexpected error: {err}");
    assert!(!fixture().join("refused.kst").exists());
    // and before the index is loaded
    assert!(!err.contains("Loading the index"), "refused only after loading the index: {err}");
}

#[test]
fn kst_refuses_stdout_and_kst_flags_without_kst() {
    let out = lookup(&["--format", "kst", "--report-misses", "-q", "q2.fa"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("--format kst needs --output"));

    let out = lookup(&["--kst-effort", "medium", "--report-misses", "-q", "q2.fa", "-o", "x.bed"]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("need --format kst"));
}

#[test]
fn kst_replaces_existing_output_and_leaves_no_temporary_file() {
    // Its own directory: the other tests write beside their outputs at the same time
    let dir = fixture().join("over");
    std::fs::create_dir_all(&dir).unwrap();
    std::fs::write(dir.join("over.kst"), b"old").unwrap();
    lookup_ok(&["--report-misses", "--format", "kst", "-q", "q2.fa", "-o", "over/over.kst"]);
    assert!(read("over/over.kst").starts_with(b"KST00001"));
    let names: Vec<_> = std::fs::read_dir(&dir).unwrap().filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().into_owned()).collect();
    assert_eq!(names, ["over.kst"], "unexpected files beside the output");
}
