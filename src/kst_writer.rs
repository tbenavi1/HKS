//! Writes lookup runs straight into a karyozip `.kst` track, with no BED in between.
//!
//! The rows handed to karyozip are exactly the rows `OutputWriter` would print with
//! `--report-misses --no-header`: same sequence names or ranks, same label names or
//! ids, same miss label. A `.kst` written here is therefore byte-identical to
//! `karyozip -e <effort> compress` of the BED that `--format bed` writes.

use std::{fs::File, io::{BufWriter, Write}, ops::Range, path::{Path, PathBuf}};

use crate::parallel_queries::RunWriter;

/// karyozip's `-e` levels that are offered here. `high` decides on the whole file and
/// is meant for archives, so it is left to the karyozip CLI.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum KstEffort {
    /// Blocks of 256 Ki runs: fast to write, and a query by read decodes one small block (karyozip -e low).
    #[default]
    Low,
    /// Whole-read blocks of about 4 Mi runs: smaller, but a query by read decodes a whole block (karyozip -e medium).
    Medium,
}

impl KstEffort {
    /// The writer options the karyozip CLI uses for the same `-e` (karyozip rust/src/main.rs).
    fn options(self) -> karyozip::writer::Options {
        use karyozip::writer::Options;
        match self {
            KstEffort::Low => Options { block: 262144, slow: false, whole: false, ..Default::default() },
            KstEffort::Medium => Options { block: 4194304, varblk: true, ..Default::default() },
        }
    }
}

/// The one place where karyozip's errors become ours. karyozip returns `String` today
/// and will return its own error type; anything that prints works here unchanged.
fn kz<T, E: std::fmt::Display>(r: Result<T, E>, path: &Path) -> Result<T, String> {
    r.map_err(|e| format!("Could not write {}: karyozip: {e}", path.display()))
}

pub struct KstRunWriter {
    path: PathBuf,
    writer: Option<karyozip::writer::Writer>,
    seq_names: Option<Vec<String>>,
    /// The label column for every label id: its name, or the id as text.
    labels: Vec<Vec<u8>>,
    miss_label: Vec<u8>,
    /// Reused for the decimal rank of each run's sequence.
    rank_buf: Vec<u8>,
    /// The first error from a row. `RunWriter` cannot return one, so it waits for `finish`.
    err: Option<String>,
}

impl KstRunWriter {
    /// `seq_names` and `color_names` as for `OutputWriter`: None prints ranks and label ids.
    /// The miss label is required: a `.kst` holds gapless runs, so misses must be written.
    pub fn new(path: &Path, seq_names: Option<Vec<String>>, color_names: Option<&[String]>, n_colors: usize, miss_label: &str, effort: KstEffort) -> Result<Self, String> {
        // The staging file karyozip keeps while blocks accumulate (unlinked at once) lives
        // beside the output, which is where the space for the result is anyway.
        let dir = parent_dir(path);
        let writer = karyozip::writer::Writer::new(effort.options(), &dir)
            .map_err(|e| format!("Could not create a temporary file in {} for {}: {e}", dir.display(), path.display()))?;
        let labels = match color_names {
            Some(names) => names.iter().map(|n| n.as_bytes().to_vec()).collect(),
            None => (0..n_colors).map(|c| c.to_string().into_bytes()).collect(),
        };
        Ok(Self {
            path: path.to_path_buf(), writer: Some(writer), seq_names, labels,
            miss_label: miss_label.as_bytes().to_vec(), rank_buf: Vec::new(), err: None,
        })
    }

    /// Writes the track: into a temporary file beside the output, renamed over it only
    /// once it is complete and synced. On failure the temporary file is removed and any
    /// earlier file at the output path is left as it was. With `index`, also writes the
    /// `.kzi` names index beside it, as `karyozip index` does.
    pub fn finish(mut self, index: bool) -> Result<(), String> {
        if let Some(e) = self.err.take() { return Err(e); }
        let index_path = match index {
            true => Some(self.path.to_str().ok_or_else(|| format!("--kst-index needs a UTF-8 output path, got {}", self.path.display()))?.to_string()),
            false => None,
        };
        let writer = self.writer.take().expect("finish is called once");
        let tmp = TempOutput::create(&self.path)?;
        let mut out = BufWriter::with_capacity(1 << 21, tmp.file.try_clone().map_err(|e| format!("Could not write {}: {e}", tmp.path.display()))?);
        let stats = kz(writer.close(&mut out), &self.path)?;
        out.into_inner().map_err(|e| e.into_error()).and_then(|f| f.sync_all())
            .map_err(|e| format!("Could not write {}: {e}", tmp.path.display()))?;
        tmp.commit(&self.path)?;
        log::info!("Wrote {} runs of {} sequences in {} blocks to {}", stats.rows, stats.reads, stats.blocks, self.path.display());

        if let Some(p) = index_path {
            let opts = karyozip::reader::ReadOpts { index: true, ..Default::default() };
            karyozip::reader::Reader::open_with(&p, &opts).build_index(&p);
        }
        Ok(())
    }

    fn row(&mut self, seq_id: isize, label: Option<usize>, range: Range<usize>) -> Result<(), String> {
        let name: &[u8] = match &self.seq_names {
            Some(names) => names[seq_id as usize].as_bytes(),
            None => {
                self.rank_buf.clear();
                write!(self.rank_buf, "{seq_id}").unwrap(); // Unwrap ok: writing to a Vec
                &self.rank_buf
            }
        };
        let label: &[u8] = match label {
            Some(c) => &self.labels[c],
            None => &self.miss_label,
        };
        let w = self.writer.as_mut().expect("no rows after finish");
        kz(w.row(name, range.len() as u64, label, range.start as u64), &self.path)
    }
}

impl RunWriter for KstRunWriter {
    fn write_header(&mut self) {} // A .kst has no header line

    fn write_run(&mut self, seq_id: isize, run_color: Option<usize>, range: Range<usize>) {
        if range.is_empty() || self.err.is_some() { return; }
        if let Err(e) = self.row(seq_id, run_color, range) {
            self.err = Some(e);
        }
    }

    fn flush(&mut self) {} // Everything is written by finish()
}

fn parent_dir(path: &Path) -> PathBuf {
    match path.parent() {
        Some(p) if !p.as_os_str().is_empty() => p.to_path_buf(),
        _ => PathBuf::from("."),
    }
}

/// A file beside the output, removed on drop unless it was renamed into place.
struct TempOutput {
    path: PathBuf,
    file: File,
    committed: bool,
}

impl TempOutput {
    fn create(out: &Path) -> Result<Self, String> {
        let base = out.file_name().ok_or_else(|| format!("Output path {} has no file name", out.display()))?;
        let dir = parent_dir(out);
        for attempt in 0..100u32 {
            let mut name = std::ffi::OsString::from(".");
            name.push(base);
            name.push(format!(".tmp.{}.{attempt}", std::process::id()));
            let path = dir.join(name);
            match std::fs::OpenOptions::new().write(true).create_new(true).open(&path) {
                Ok(file) => return Ok(Self { path, file, committed: false }),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => return Err(format!("Could not create output file {}: {e}", path.display())),
            }
        }
        Err(format!("Could not create a temporary file beside {}", out.display()))
    }

    fn commit(mut self, out: &Path) -> Result<(), String> {
        std::fs::rename(&self.path, out).map_err(|e| format!("Could not rename {} to {}: {e}", self.path.display(), out.display()))?;
        self.committed = true;
        Ok(())
    }
}

impl Drop for TempOutput {
    fn drop(&mut self) {
        if !self.committed { let _ = std::fs::remove_file(&self.path); }
    }
}
