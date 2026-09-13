//! Latin squares of order `n`: full counts, and one representative per
//! isomorphism class.
//!
//! Isomorphism is simultaneous relabelling of rows, columns and symbols by one
//! permutation `p`: `T'[p(x)][p(y)] = p(T[x][y])`. The representative is the
//! row-major lexicographically least table in its class. `permutations` and
//! the canonicity idea are adapted from parsimagma (`src/finite.rs` at commit
//! 2a75a1b, Apache-2.0); the early-exit comparison, the automorphism count and
//! the two pruning rules below are new here.
//!
//! Pruning, both *necessary* conditions for a least table, so no class is
//! lost:
//! 1. `T[0][0] <= 1`. If some `x` is idempotent, relabel it to 0 and the first
//!    cell becomes 0. Otherwise send any `x` to 0 and `x * x` to 1.
//! 2. Row 0 is least among its conjugates `p r p^-1` with `p(0) = 0`, since
//!    such a relabelling rewrites row 0 alone as that conjugate.

use rayon::prelude::*;

pub const MAX_N: usize = 8;

pub struct PermTable {
    pub n: usize,
    pub p: Vec<Vec<u8>>,
    pub inv: Vec<Vec<u8>>,
}

/// All permutations of `0..n` (from parsimagma).
pub fn permutations(n: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut cur: Vec<u8> = (0..n as u8).collect();
    permute(&mut cur, 0, &mut out);
    out
}

fn permute(cur: &mut Vec<u8>, k: usize, out: &mut Vec<Vec<u8>>) {
    if k == cur.len() {
        out.push(cur.clone());
        return;
    }
    for i in k..cur.len() {
        cur.swap(k, i);
        permute(cur, k + 1, out);
        cur.swap(k, i);
    }
}

impl PermTable {
    pub fn new(n: usize) -> PermTable {
        let p = permutations(n);
        let inv = p
            .iter()
            .map(|q| {
                let mut v = vec![0u8; n];
                for (i, &x) in q.iter().enumerate() {
                    v[x as usize] = i as u8;
                }
                v
            })
            .collect();
        PermTable { n, p, inv }
    }
}

/// If `t` is the least table in its isomorphism class, the order of its
/// automorphism group; otherwise `None`.
pub fn canonical_aut(t: &[u8], pt: &PermTable) -> Option<u32> {
    let n = pt.n;
    let mut aut = 0;
    'perm: for (p, inv) in pt.p.iter().zip(&pt.inv) {
        for i in 0..n {
            let x = inv[i] as usize * n;
            for j in 0..n {
                let v = p[t[x + inv[j] as usize] as usize];
                let w = t[i * n + j];
                if v < w {
                    return None;
                }
                if v > w {
                    continue 'perm;
                }
            }
        }
        aut += 1;
    }
    Some(aut)
}

/// The relabelled table `T'` for permutation index `k`.
pub fn relabel(t: &[u8], pt: &PermTable, k: usize) -> Vec<u8> {
    let n = pt.n;
    let p = &pt.p[k];
    let mut out = vec![0u8; n * n];
    for x in 0..n {
        for y in 0..n {
            out[p[x] as usize * n + p[y] as usize] = p[t[x * n + y] as usize];
        }
    }
    out
}

/// Rows `r` with `r[0] <= 1` that are least among their conjugates `p r p^-1`
/// over permutations `p` fixing 0. Computed orbit by orbit, so the cost is
/// (number of orbits) x (n-1)! rather than n! x (n-1)!.
pub fn row0_candidates(pt: &PermTable) -> Vec<Vec<u8>> {
    let n = pt.n;
    let stab: Vec<(&Vec<u8>, &Vec<u8>)> = pt.p.iter().zip(&pt.inv).filter(|(p, _)| p[0] == 0).collect();
    let mut seen: std::collections::HashSet<Vec<u8>> = std::collections::HashSet::new();
    let mut out = Vec::new();
    for r in &pt.p {
        if r[0] > 1 || seen.contains(r) {
            continue;
        }
        let mut least = r.clone();
        for (p, inv) in &stab {
            let conj: Vec<u8> = (0..n).map(|i| p[r[inv[i] as usize] as usize]).collect();
            if conj < least {
                least = conj.clone();
            }
            seen.insert(conj);
        }
        if least[0] <= 1 && !out.contains(&least) {
            out.push(least);
        }
    }
    out.sort();
    out
}

/// Prefixes (rows 0 and 1) used to split the search into parallel tasks.
fn tasks(n: usize, row0s: &[Vec<u8>], all_rows: &[Vec<u8>]) -> Vec<Vec<u8>> {
    if n == 1 {
        return row0s.to_vec();
    }
    let mut out = Vec::new();
    for r0 in row0s {
        for r1 in all_rows {
            if r0.iter().zip(r1).all(|(a, b)| a != b) {
                let mut v = r0.clone();
                v.extend_from_slice(r1);
                out.push(v);
            }
        }
    }
    out
}

struct Search {
    n: usize,
    full: u16,
    rows: [u16; MAX_N],
    cols: [u16; MAX_N],
    t: [u8; MAX_N * MAX_N],
}

impl Search {
    fn new(n: usize, prefix: &[u8]) -> Search {
        assert!(n >= 1 && n <= MAX_N);
        let mut s = Search {
            n,
            full: ((1u32 << n) - 1) as u16,
            rows: [0; MAX_N],
            cols: [0; MAX_N],
            t: [0; MAX_N * MAX_N],
        };
        for (pos, &v) in prefix.iter().enumerate() {
            let (i, j) = (pos / n, pos % n);
            assert!(s.rows[i] & (1 << v) == 0 && s.cols[j] & (1 << v) == 0, "bad prefix");
            s.rows[i] |= 1 << v;
            s.cols[j] |= 1 << v;
            s.t[pos] = v;
        }
        s
    }

    fn walk(&mut self, pos: usize, leaf: &mut impl FnMut(&[u8])) {
        let n = self.n;
        if pos == n * n {
            leaf(&self.t[..n * n]);
            return;
        }
        let (i, j) = (pos / n, pos % n);
        let mut free = !(self.rows[i] | self.cols[j]) & self.full;
        while free != 0 {
            let v = free.trailing_zeros() as u8;
            free &= free - 1;
            let bit = 1u16 << v;
            self.rows[i] |= bit;
            self.cols[j] |= bit;
            self.t[pos] = v;
            self.walk(pos + 1, leaf);
            self.rows[i] &= !bit;
            self.cols[j] &= !bit;
        }
    }
}

/// A representative of an isomorphism class.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Canon {
    pub table: Vec<u8>,
    pub aut: u32,
}

/// One least table per isomorphism class of quasigroups of order `n`, sorted.
pub fn iso_classes(n: usize) -> Vec<Canon> {
    let pt = PermTable::new(n);
    let row0s = row0_candidates(&pt);
    let jobs = tasks(n, &row0s, &pt.p);
    let mut out: Vec<Canon> = jobs
        .par_iter()
        .flat_map_iter(|prefix| {
            let mut s = Search::new(n, prefix);
            let mut found = Vec::new();
            s.walk(prefix.len(), &mut |t| {
                if let Some(aut) = canonical_aut(t, &pt) {
                    found.push(Canon { table: t.to_vec(), aut });
                }
            });
            found
        })
        .collect();
    out.sort();
    out
}

/// Number of Latin squares of order `n`, by plain backtracking (no symmetry).
pub fn count_latin(n: usize) -> u64 {
    let pt = PermTable::new(n);
    let jobs = tasks(n, &pt.p, &pt.p);
    jobs.par_iter()
        .map(|prefix| {
            let mut s = Search::new(n, prefix);
            let mut c = 0u64;
            s.walk(prefix.len(), &mut |_| c += 1);
            c
        })
        .sum()
}

/// Every Latin square of order `n`, row-major. Only sensible for small `n`.
pub fn all_latin(n: usize) -> Vec<Vec<u8>> {
    let mut out = Vec::new();
    let mut s = Search::new(n, &[]);
    s.walk(0, &mut |t| out.push(t.to_vec()));
    out
}

pub fn factorial(n: usize) -> u64 {
    (1..=n as u64).product()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn small_counts() {
        let iso: Vec<usize> = (1..=4).map(|n| iso_classes(n).len()).collect();
        assert_eq!(iso, vec![1, 1, 5, 35]);
        let lat: Vec<u64> = (1..=4).map(count_latin).collect();
        assert_eq!(lat, vec![1, 2, 12, 576]);
    }

    /// The orbit-wise candidate filter agrees with the direct definition.
    #[test]
    fn row0_candidates_match_definition() {
        for n in 1..=6 {
            let pt = PermTable::new(n);
            let direct: Vec<Vec<u8>> = pt
                .p
                .iter()
                .filter(|r| r[0] <= 1)
                .filter(|r| {
                    pt.p.iter().zip(&pt.inv).filter(|(p, _)| p[0] == 0).all(|(p, inv)| {
                        let conj: Vec<u8> = (0..n).map(|i| p[r[inv[i] as usize] as usize]).collect();
                        conj.as_slice() >= r.as_slice()
                    })
                })
                .cloned()
                .collect();
            let mut direct = direct;
            direct.sort();
            assert_eq!(direct, row0_candidates(&pt), "order {n}");
        }
    }

    /// Pruning must not lose classes: compare with an unpruned sweep.
    #[test]
    fn pruning_is_sound_at_order_4() {
        let pt = PermTable::new(4);
        let mut brute: Vec<Canon> = all_latin(4)
            .into_iter()
            .filter_map(|t| canonical_aut(&t, &pt).map(|aut| Canon { table: t, aut }))
            .collect();
        brute.sort();
        assert_eq!(brute, iso_classes(4));
    }
}
