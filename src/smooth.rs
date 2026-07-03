//! Hierarchy-aware smoothing for HKS lookup output.
//!
//! Implements Algorithm S1 from the HKS paper: a two-phase scan that identifies
//! windows exhibiting a specific → general → specific pattern in the category
//! hierarchy and reassigns interior intervals to the LCA of the flanking anchors.

use std::collections::HashSet;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};

use crate::lca_tree::LcaTree;

// ---------------------------------------------------------------------------
// Interval representation
// ---------------------------------------------------------------------------

#[derive(Clone)]
pub struct Interval {
    pub start: u64,
    pub end: u64,
    pub feature: usize,     // node ID in the hierarchy
    pub originally_none: bool, // true if this interval was "none" / unmatched in the input
}

// ---------------------------------------------------------------------------
// Core smoothing algorithm (port of Algorithm S1)
// ---------------------------------------------------------------------------

/// Smooth a single query's intervals in-place using the hierarchy.
/// Returns the number of feature reassignments made.
pub fn smooth_intervals(intervals: &mut Vec<Interval>, tree: &LcaTree, max_gap: u64) -> u64 {
    if intervals.len() < 2 {
        return 0;
    }
    let mut total_reassignments = 0u64;
    let n = intervals.len();

    loop {
        let mut changed = false;
        let mut was_related = vec![false; n];
        let mut i = 0;

        while i < n {
            let window_start = i;
            let mut related: Vec<usize> = vec![i];
            let mut first_unrelated: Option<usize> = None;

            // ------------------------------------------------------------------
            // Ascending phase: scan rightward accepting features that are
            // ancestors of last_rel_feat (i.e. more general / closer to root).
            // Unrelated features on other branches are skipped but their ancestor
            // paths are added to `disallowed` so we stop if we'd cross into them.
            // ------------------------------------------------------------------
            let mut disallowed: HashSet<usize> = HashSet::new();
            let mut last_rel_idx = i;
            let mut last_rel_feat = intervals[i].feature;
            let mut last_rel_end = intervals[i].end;
            let mut j = i;

            while j + 1 < n {
                let nf = intervals[j + 1].feature;
                let ns = intervals[j + 1].start;

                if ns > last_rel_end + max_gap {
                    break;
                }

                if tree.is_ancestor(nf, last_rel_feat) {
                    // nf IS an ancestor of last_rel_feat → nf is more general → extend
                    if disallowed.contains(&nf) {
                        break;
                    }
                    last_rel_feat = nf;
                    last_rel_end = intervals[j + 1].end;
                    last_rel_idx = j + 1;
                    related.push(j + 1);
                    j += 1;
                } else if !tree.is_ancestor(last_rel_feat, nf) {
                    // Neither is ancestor of the other → different branches
                    disallowed.extend(tree.ancestors(nf));
                    if first_unrelated.is_none() && !was_related[j + 1] {
                        first_unrelated = Some(j + 1);
                    }
                    j += 1;
                } else {
                    // last_rel_feat IS ancestor of nf → nf is more specific → stop
                    break;
                }
            }

            let mut window_end = last_rel_idx;

            // ------------------------------------------------------------------
            // Descending phase: continue rightward from the peak, accepting
            // features that are descendants of last_rel_feat (more specific).
            // ------------------------------------------------------------------
            let peak_feat = intervals[window_end].feature;
            let mut k = window_end;

            while k + 1 < n {
                let nf = intervals[k + 1].feature;
                let ns = intervals[k + 1].start;

                if ns > last_rel_end + max_gap {
                    break;
                }

                if tree.is_ancestor(last_rel_feat, nf) {
                    // last_rel_feat IS ancestor of nf → nf is more specific → extend
                    last_rel_feat = nf;
                    last_rel_end = intervals[k + 1].end;
                    last_rel_idx = k + 1;
                    related.push(k + 1);
                    k += 1;
                } else if !tree.is_ancestor(nf, last_rel_feat) {
                    // Neither is ancestor → unrelated
                    if tree.is_ancestor(peak_feat, nf) {
                        // nf is a descendant of peak → would restart a new ascending window → stop
                        break;
                    }
                    if first_unrelated.is_none() && !was_related[k + 1] {
                        first_unrelated = Some(k + 1);
                    }
                    k += 1;
                } else {
                    // nf IS ancestor of last_rel_feat → nf is more general → stop
                    break;
                }
            }

            window_end = last_rel_idx;

            // Drop the last element from related (boundary stays unchanged)
            if related.len() >= 2 {
                related.pop();
            }
            for &idx in &related {
                was_related[idx] = true;
            }

            // ------------------------------------------------------------------
            // Reassignment: replace interior features with LCA(left, right) when
            // they are strictly more general (shallower) than the LCA.
            // ------------------------------------------------------------------
            if window_end > window_start {
                let left = intervals[window_start].feature;
                let right = intervals[window_end].feature;
                let lca = tree.lca(left, right);
                for w in (window_start + 1)..window_end {
                    let orig = intervals[w].feature;
                    // is_ancestor(orig, lca) means orig IS an ancestor of lca,
                    // i.e. orig is more general than lca → replace with lca
                    if tree.is_ancestor(orig, lca) && orig != lca {
                        intervals[w].feature = lca;
                        intervals[w].originally_none = false;
                        changed = true;
                        total_reassignments += 1;
                    }
                }
            }

            // Advance i
            i = match first_unrelated {
                Some(fu) if fu > window_end && window_end > window_start => window_end,
                Some(fu) => fu,
                None if window_end > i => window_end,
                None => {
                    let mut next = i;
                    while next < n && was_related[next] {
                        next += 1;
                    }
                    next
                }
            };
        }

        if !changed {
            break;
        }
    }
    total_reassignments
}

// ---------------------------------------------------------------------------
// Merge adjacent contiguous intervals with the same feature
// ---------------------------------------------------------------------------

/// Merges adjacent intervals that have the same feature and are contiguous.
/// Returns (merged intervals, number of intervals eliminated).
pub fn merge_intervals(intervals: Vec<Interval>) -> (Vec<Interval>, u64) {
    if intervals.is_empty() {
        return (intervals, 0);
    }
    let n_in = intervals.len();
    let mut out: Vec<Interval> = Vec::with_capacity(n_in);
    let mut cur = intervals.into_iter();
    let mut current = cur.next().unwrap();
    for next in cur {
        if next.feature == current.feature
            && next.start == current.end
            && next.originally_none == current.originally_none
        {
            current.end = next.end;
        } else {
            out.push(current);
            current = next;
        }
    }
    out.push(current);
    let eliminated = (n_in - out.len()) as u64;
    (out, eliminated)
}

// ---------------------------------------------------------------------------
// Streaming smooth processor: parse TSV → smooth per query → write TSV
// ---------------------------------------------------------------------------

#[derive(Default)]
pub struct SmoothStats {
    pub reads_processed: u64,
    pub intervals_in: u64,
    pub intervals_smoothed: u64,
    pub intervals_merged: u64,
    pub intervals_out: u64,
}

/// Resolve a feature token from the TSV input into a node ID.
/// Returns (node_id, originally_none).
fn resolve_feature(
    token: &str,
    name_to_id: &std::collections::HashMap<String, usize>,
    root_id: usize,
    uses_names: bool,
) -> (usize, bool) {
    if uses_names {
        if token == "none" {
            (root_id, true)
        } else {
            let id = name_to_id.get(token)
                .unwrap_or_else(|| panic!("Unknown feature name in input: {:?}", token));
            (*id, false)
        }
    } else {
        // Numeric ID mode
        if token == "-" {
            (root_id, true)
        } else {
            let id: usize = token.parse()
                .unwrap_or_else(|_| panic!("Cannot parse feature ID: {:?}", token));
            (id, false)
        }
    }
}

/// Format a feature for output.
fn format_feature(
    feature: usize,
    originally_none: bool,
    names: &[String],
    root_id: usize,
    uses_names: bool,
) -> String {
    if originally_none && feature == root_id {
        // Was none and smoothing didn't resolve it → keep as none
        if uses_names { "none".to_string() } else { "-".to_string() }
    } else if uses_names {
        names[feature].to_string()
    } else {
        feature.to_string()
    }
}

/// Flush a completed query: smooth → merge → write.
fn flush_query(
    query_id: &str,
    buf: &mut Vec<Interval>,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    uses_names: bool,
    max_gap: u64,
    writer: &mut impl Write,
    stats: &mut SmoothStats,
) {
    let n_in = buf.len() as u64;
    let reassigned = smooth_intervals(buf, tree, max_gap);
    let (merged, eliminated) = merge_intervals(std::mem::take(buf));
    let n_out = merged.len() as u64;

    for iv in &merged {
        let feat_str = format_feature(iv.feature, iv.originally_none, names, root_id, uses_names);
        writeln!(writer, "{}\t{}\t{}\t{}", query_id, iv.start, iv.end, feat_str)
            .expect("write error");
    }

    stats.reads_processed += 1;
    stats.intervals_in += n_in;
    stats.intervals_smoothed += reassigned;
    stats.intervals_merged += eliminated;
    stats.intervals_out += n_out;
}

/// Run the smoothing pipeline on TSV input.
pub fn run_smooth(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
) -> SmoothStats {
    // Build name → id lookup
    let name_to_id: std::collections::HashMap<String, usize> = names.iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id))
        .collect();

    let reader = BufReader::new(input);
    let mut writer = BufWriter::new(output);
    let mut stats = SmoothStats::default();
    let mut buf: Vec<Interval> = Vec::new();
    let mut current_query = String::new();
    let mut uses_names = false; // determined from header
    let mut header_seen = false;

    for line in reader.lines() {
        let line = line.expect("IO error reading input");
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }

        // Parse header
        if !header_seen {
            if trimmed.starts_with("query_rank") || trimmed.starts_with("query_name") {
                header_seen = true;
                uses_names = trimmed.contains("label_name");
                // Pass header through
                writeln!(writer, "{}", trimmed).expect("write error");
                continue;
            }
            // No header line — treat as data
            header_seen = true;
        }

        let mut cols = trimmed.splitn(4, '\t');
        let query_id = cols.next().expect("missing query column");
        let start: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad start coordinate");
        let end: u64 = cols.next()
            .and_then(|s| s.parse().ok())
            .expect("bad end coordinate");
        let feat_token = cols.next().expect("missing feature column");

        let (feature, originally_none) = resolve_feature(feat_token, &name_to_id, root_id, uses_names);

        // Flush on query boundary
        if query_id != current_query {
            if !buf.is_empty() {
                log::info!("Smoothing {}", current_query);
                flush_query(&current_query, &mut buf, tree, names, root_id, uses_names, max_gap, &mut writer, &mut stats);
            }
            current_query = query_id.to_string();
        }

        buf.push(Interval { start, end, feature, originally_none });
    }

    // Flush final query
    if !buf.is_empty() {
        log::info!("Smoothing {}", current_query);
        flush_query(&current_query, &mut buf, tree, names, root_id, uses_names, max_gap, &mut writer, &mut stats);
    }

    writer.flush().expect("flush error");
    stats
}

impl SmoothStats {
    /// Fold one group's partial stats into this accumulator.
    fn merge(&mut self, o: &SmoothStats) {
        self.reads_processed += o.reads_processed;
        self.intervals_in += o.intervals_in;
        self.intervals_smoothed += o.intervals_smoothed;
        self.intervals_merged += o.intervals_merged;
        self.intervals_out += o.intervals_out;
    }
}

/// Smooth one query group given as a raw byte slice of the input TSV.
///
/// Every line in `slice` shares one query id. Parses the slice, runs
/// smooth + merge, and renders the result to a freshly allocated byte buffer.
/// Doing the parse *and* the output formatting here (not just the smoothing
/// compute) is what lets the parallel path actually scale: parsing and
/// formatting tens of millions of lines dominate the wall time, so both must
/// run on the worker thread. Parsing borrows byte subslices of the input —
/// no per-line `String` allocation, unlike `reader.lines()`.
fn smooth_group(
    slice: &[u8],
    tree: &LcaTree,
    names: &[String],
    name_to_id: &std::collections::HashMap<String, usize>,
    root_id: usize,
    uses_names: bool,
    max_gap: u64,
) -> (SmoothStats, Vec<u8>) {
    let mut query_id: &str = "";
    let mut have_qid = false;
    let mut intervals: Vec<Interval> = Vec::new();

    for raw in slice.split(|&b| b == b'\n') {
        let line = raw.trim_ascii();
        if line.is_empty() {
            continue;
        }
        let mut cols = line.splitn(4, |&b| b == b'\t');
        let qid_b = cols.next().expect("missing query column");
        let start_b = cols.next().expect("missing start coordinate");
        let end_b = cols.next().expect("missing end coordinate");
        let feat_b = cols.next().expect("missing feature column");

        if !have_qid {
            query_id = std::str::from_utf8(qid_b).expect("non-UTF8 query name");
            have_qid = true;
        }
        let start: u64 = std::str::from_utf8(start_b).ok()
            .and_then(|s| s.parse().ok())
            .expect("bad start coordinate");
        let end: u64 = std::str::from_utf8(end_b).ok()
            .and_then(|s| s.parse().ok())
            .expect("bad end coordinate");
        let feat_token = std::str::from_utf8(feat_b).expect("non-UTF8 feature name");

        let (feature, originally_none) = resolve_feature(feat_token, name_to_id, root_id, uses_names);
        intervals.push(Interval { start, end, feature, originally_none });
    }

    let n_in = intervals.len() as u64;
    let reassigned = smooth_intervals(&mut intervals, tree, max_gap);
    let (merged, eliminated) = merge_intervals(intervals);
    let n_out = merged.len() as u64;

    // Grow from empty rather than reserving slice.len(): smoothing+merging
    // often collapses a huge input group (e.g. a whole chromosome) to a handful
    // of intervals, so reserving the input size would mmap hundreds of MB per
    // worker for a few bytes of output — concurrent large reservations serialize
    // on the kernel's mmap lock and cause a hard regression at high thread counts.
    let mut out: Vec<u8> = Vec::new();
    for iv in &merged {
        let feat_str = format_feature(iv.feature, iv.originally_none, names, root_id, uses_names);
        writeln!(out, "{}\t{}\t{}\t{}", query_id, iv.start, iv.end, feat_str)
            .expect("formatting error");
    }

    let stats = SmoothStats {
        reads_processed: 1,
        intervals_in: n_in,
        intervals_smoothed: reassigned,
        intervals_merged: eliminated,
        intervals_out: n_out,
    };
    (stats, out)
}

/// Parallel smoothing pipeline: smooth each query (sequence) on its own thread.
///
/// Each worker **parses, smooths, and formats its own slice** of the input.
/// The main thread does only a cheap boundary scan (locate per-query byte
/// ranges — no integer parsing or allocation) up front, then emits the results.
/// This moves the two dominant serial costs — parsing and output formatting of
/// tens of millions of lines — into the parallel region (the earlier version
/// left both on the main thread and parallelised only the smoothing compute,
/// which capped speedup near ~1.4x regardless of thread count).
///
/// Two output-ordering modes:
///
/// * `preserve_order == true` (assemblies: few, long sequences where output
///   order must match the input). Worker buffers are collected and written in
///   input order — byte-for-byte identical to [`run_smooth`].
/// * `preserve_order == false` (reads: many short sequences where order is
///   irrelevant). Worker buffers are streamed out in completion order via a
///   channel, so the writer never blocks on a straggler and peak memory is
///   bounded by the in-flight backlog rather than the whole output.
///
/// `n_threads` controls the worker parallelism (a local rayon pool). The whole
/// input is read into memory once and parsed from borrowed byte slices.
pub fn run_smooth_parallel(
    input: impl Read,
    output: impl Write,
    tree: &LcaTree,
    names: &[String],
    root_id: usize,
    max_gap: u64,
    n_threads: usize,
    preserve_order: bool,
) -> SmoothStats {
    use rayon::prelude::*;

    // Build name → id lookup.
    let name_to_id: std::collections::HashMap<String, usize> = names.iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id))
        .collect();

    // ---- Read the whole input into memory (one bulk read; I/O bound) ----
    let mut reader = BufReader::new(input);
    let mut data: Vec<u8> = Vec::new();
    reader.read_to_end(&mut data).expect("IO error reading input");
    let mut writer = BufWriter::new(output);

    // ---- Boundary scan: extract the header and one byte range per maximal
    //      run of lines sharing a query id. Cheap — only locates newlines and
    //      the first tab of each line; no integer parsing or allocation. ----
    let n = data.len();
    let mut groups: Vec<(usize, usize)> = Vec::new();
    let mut header: Option<Vec<u8>> = None;
    let mut uses_names = false;
    let mut cur_qid: Option<&[u8]> = None;
    let mut cur_start = 0usize;
    let mut first_content = true;
    let mut i = 0usize;
    while i < n {
        let j = match data[i..].iter().position(|&b| b == b'\n') {
            Some(k) => i + k,
            None => n,
        };
        let line = data[i..j].trim_ascii();
        let next_i = j + 1;
        if line.is_empty() {
            i = next_i;
            continue;
        }
        if first_content {
            first_content = false;
            if line.starts_with(b"query_rank") || line.starts_with(b"query_name") {
                uses_names = line.windows(b"label_name".len()).any(|w| w == b"label_name");
                header = Some(line.to_vec());
                i = next_i;
                continue;
            }
            // Not a header — fall through and treat this line as data.
        }
        let qid = match line.iter().position(|&b| b == b'\t') {
            Some(k) => &line[..k],
            None => line,
        };
        if cur_qid != Some(qid) {
            if cur_qid.is_some() {
                groups.push((cur_start, i));
            }
            cur_qid = Some(qid);
            cur_start = i;
        }
        i = next_i;
    }
    if cur_qid.is_some() {
        groups.push((cur_start, n));
    }

    // Header goes out first, matching run_smooth.
    if let Some(h) = &header {
        writer.write_all(h).expect("write error");
        writer.write_all(b"\n").expect("write error");
    }

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_threads)
        .build()
        .expect("failed to build rayon thread pool");

    let mut stats = SmoothStats::default();

    if preserve_order {
        // Assembly mode: collect worker buffers, emit in input order.
        let results: Vec<(SmoothStats, Vec<u8>)> = pool.install(|| {
            groups
                .par_iter()
                .map(|&(s, e)| {
                    smooth_group(&data[s..e], tree, names, &name_to_id, root_id, uses_names, max_gap)
                })
                .collect()
        });
        for (st, bytes) in &results {
            writer.write_all(bytes).expect("write error");
            stats.merge(st);
        }
    } else {
        // Reads mode: stream worker buffers out in completion order. The output
        // handle isn't Send, so it stays on this thread; workers hand back
        // (stats, bytes) over a channel while a scoped thread drives the pool.
        let data_ref = &data;
        let groups_ref = &groups;
        let nti = &name_to_id;
        std::thread::scope(|scope| {
            let (tx, rx) = std::sync::mpsc::channel::<(SmoothStats, Vec<u8>)>();
            scope.spawn(move || {
                pool.install(|| {
                    groups_ref.par_iter().for_each_with(tx, |tx, &(s, e)| {
                        let res = smooth_group(
                            &data_ref[s..e], tree, names, nti, root_id, uses_names, max_gap,
                        );
                        tx.send(res).expect("channel send error");
                    });
                });
            });
            for (st, bytes) in rx {
                writer.write_all(&bytes).expect("write error");
                stats.merge(&st);
            }
        });
    }

    writer.flush().expect("flush error");
    stats
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lca_tree::LcaTree;

    /// Build a small tree:   root(3) → A(2) → { B(0), C(1) }
    fn cousin_tree() -> LcaTree {
        // Edges: B→A, C→A, A→root
        LcaTree::new(4, vec![(0, 2), (1, 2), (2, 3)]).unwrap()
    }

    fn iv(feature: usize, start: u64, end: u64) -> Interval {
        Interval { start, end, feature, originally_none: feature == 3 }
    }

    /// Canonical case: B, root, C  →  B, A, C
    /// (a too-general interior interval is promoted to LCA of its specific flankers)
    #[test]
    fn promotes_general_interior_to_lca() {
        let tree = cousin_tree();
        let root = tree.root(); // node 3
        let a = 2usize;
        let b = 0usize;
        let c = 1usize;

        let mut intervals = vec![
            iv(b, 0, 100),
            iv(root, 100, 200),
            iv(c, 200, 300),
        ];

        let reassigned = smooth_intervals(&mut intervals, &tree, 1000);
        assert_eq!(reassigned, 1, "expected exactly one reassignment");
        assert_eq!(intervals[0].feature, b,    "left anchor unchanged");
        assert_eq!(intervals[1].feature, a,    "interior promoted to LCA(B,C) = A");
        assert_eq!(intervals[2].feature, c,    "right anchor unchanged");
    }

    /// Already-at-LCA interior should not be touched.
    #[test]
    fn no_change_when_interior_already_at_lca() {
        let tree = cousin_tree();
        let a = 2usize;
        let b = 0usize;
        let c = 1usize;

        let mut intervals = vec![
            iv(b, 0, 100),
            iv(a, 100, 200),
            iv(c, 200, 300),
        ];

        let reassigned = smooth_intervals(&mut intervals, &tree, 1000);
        assert_eq!(reassigned, 0, "nothing to promote when interior is already LCA");
        assert_eq!(intervals[1].feature, a);
    }

    /// The parallel pipeline must produce byte-identical output to the
    /// sequential one across multiple queries (some needing promotion, one
    /// with a `none` run, one single-interval).
    #[test]
    fn parallel_matches_sequential() {
        let tree = cousin_tree(); // B(0), C(1) → A(2) → root(3)
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        let input = "query_name\tfrom_kmer\tto_kmer\tlabel_name\n\
                     seq1\t0\t100\tB\n\
                     seq1\t100\t200\troot\n\
                     seq1\t200\t300\tC\n\
                     seq2\t0\t50\tnone\n\
                     seq2\t50\t150\tB\n\
                     seq2\t150\t250\tC\n\
                     seq3\t0\t100\tA\n";

        let mut out_seq: Vec<u8> = Vec::new();
        let s1 = run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000);

        for nt in [1usize, 2, 4] {
            let mut out_par: Vec<u8> = Vec::new();
            let s2 =
                run_smooth_parallel(input.as_bytes(), &mut out_par, &tree, &names, root, 1000, nt, true);
            assert_eq!(out_seq, out_par, "parallel output differs at n_threads={nt}");
            assert_eq!(s1.reads_processed, s2.reads_processed);
            assert_eq!(s1.intervals_in, s2.intervals_in);
            assert_eq!(s1.intervals_smoothed, s2.intervals_smoothed);
            assert_eq!(s1.intervals_out, s2.intervals_out);
        }
    }

    /// In `--no-preserve-order` mode the output can appear in any per-sequence
    /// order, but the header must be preserved and the *set* of data lines must
    /// be identical to the sequential run.
    #[test]
    fn unordered_matches_sequential_as_set() {
        let tree = cousin_tree();
        let root = tree.root();
        let names = vec![
            "B".to_string(),
            "C".to_string(),
            "A".to_string(),
            "root".to_string(),
        ];
        let input = "query_name\tfrom_kmer\tto_kmer\tlabel_name\n\
                     seq1\t0\t100\tB\n\
                     seq1\t100\t200\troot\n\
                     seq1\t200\t300\tC\n\
                     seq2\t0\t50\tnone\n\
                     seq2\t50\t150\tB\n\
                     seq2\t150\t250\tC\n\
                     seq3\t0\t100\tA\n";

        let mut out_seq: Vec<u8> = Vec::new();
        run_smooth(input.as_bytes(), &mut out_seq, &tree, &names, root, 1000);
        let mut out_un: Vec<u8> = Vec::new();
        run_smooth_parallel(input.as_bytes(), &mut out_un, &tree, &names, root, 1000, 4, false);

        let seq = String::from_utf8(out_seq).unwrap();
        let un = String::from_utf8(out_un).unwrap();
        // Header identical and first.
        assert_eq!(seq.lines().next(), un.lines().next(), "header differs");
        // Data-line sets equal (order-insensitive).
        let sorted = |s: &str| {
            let mut v: Vec<&str> = s.lines().skip(1).collect();
            v.sort_unstable();
            v.join("\n")
        };
        assert_eq!(sorted(&seq), sorted(&un), "unordered data-line set differs");
    }
}
