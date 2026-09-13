//! Exhaustive search, up to isomorphism, for quasigroups of one order that
//! satisfy a given set of laws.
//!
//! Cells are filled row-major. Placing `x * y = z` defines one entry in each of
//! the three operation tables (`y // z = x`, `z \\ x = y`). After every
//! placement, every instance of every required law whose two sides can both be
//! computed from defined entries is checked, and the branch is cut on the first
//! disagreement. This only removes partial tables that no completion can rescue,
//! so the search is complete. Symmetry is broken with the enumerator's necessary
//! conditions on row 0, and each complete table must pass the full
//! lexicographic canonicity test, so every isomorphism class of models appears
//! exactly once, with its automorphism count.

use crate::enumerate::{PermTable, canonical_aut, row0_candidates};
use crate::law::Law;
use crate::qg::{Compiled, Engine, Ins, MAX_VARS, Quasigroup, STACK};
use crate::sig::Signature;
use rayon::prelude::*;

pub const MAX_ORDER: usize = 9;
const UNK: u8 = u8::MAX;

#[derive(Clone)]
struct State<'a> {
    n: usize,
    full: u16,
    rows: [u16; MAX_ORDER],
    cols: [u16; MAX_ORDER],
    ops: [[u8; MAX_ORDER * MAX_ORDER]; 3],
    laws: &'a [Compiled],
}

impl<'a> State<'a> {
    fn new(n: usize, laws: &'a [Compiled]) -> State<'a> {
        assert!((1..=MAX_ORDER).contains(&n), "order {n} outside 1..={MAX_ORDER}");
        State {
            n,
            full: ((1u32 << n) - 1) as u16,
            rows: [0; MAX_ORDER],
            cols: [0; MAX_ORDER],
            ops: [[UNK; MAX_ORDER * MAX_ORDER]; 3],
            laws,
        }
    }

    #[inline]
    fn place(&mut self, pos: usize, v: u8) {
        let n = self.n;
        let (i, j) = (pos / n, pos % n);
        let bit = 1u16 << v;
        self.rows[i] |= bit;
        self.cols[j] |= bit;
        self.ops[0][pos] = v;
        self.ops[1][j * n + v as usize] = i as u8; // j // v = i  iff  i * j = v
        self.ops[2][v as usize * n + i] = j as u8; // v \\ i = j  iff  i * j = v
    }

    #[inline]
    fn unplace(&mut self, pos: usize) {
        let n = self.n;
        let (i, j) = (pos / n, pos % n);
        let v = self.ops[0][pos];
        let bit = 1u16 << v;
        self.rows[i] &= !bit;
        self.cols[j] &= !bit;
        self.ops[0][pos] = UNK;
        self.ops[1][j * n + v as usize] = UNK;
        self.ops[2][v as usize * n + i] = UNK;
    }

    #[inline]
    fn eval(&self, prog: &[Ins], env: &[u8; MAX_VARS]) -> u8 {
        let mut st = [0u8; STACK];
        let mut sp = 0;
        for ins in prog {
            match *ins {
                Ins::Var(s) => {
                    st[sp] = env[s as usize];
                    sp += 1;
                }
                Ins::Op(o) => {
                    let b = st[sp - 1] as usize;
                    let a = st[sp - 2] as usize;
                    let v = self.ops[o as usize][a * self.n + b];
                    if v == UNK {
                        return UNK;
                    }
                    sp -= 1;
                    st[sp - 1] = v;
                }
            }
        }
        st[0]
    }

    /// No instance of a required law is already violated.
    fn consistent(&self) -> bool {
        let n = self.n as u8;
        for law in self.laws {
            let mut env = [0u8; MAX_VARS];
            'inst: loop {
                let l = self.eval(&law.lhs, &env);
                if l != UNK {
                    let r = self.eval(&law.rhs, &env);
                    if r != UNK && r != l {
                        return false;
                    }
                }
                let mut i = 0;
                loop {
                    if i == law.k {
                        break 'inst;
                    }
                    env[i] += 1;
                    if env[i] < n {
                        break;
                    }
                    env[i] = 0;
                    i += 1;
                }
            }
        }
        true
    }

    fn table(&self) -> Vec<u8> {
        self.ops[0][..self.n * self.n].to_vec()
    }

    /// Depth-first over cells `pos..`; `leaf` receives complete tables.
    fn walk(&mut self, pos: usize, nodes: &mut u64, leaf: &mut impl FnMut(&State)) {
        *nodes += 1;
        if !self.consistent() {
            return;
        }
        let n = self.n;
        if pos == n * n {
            leaf(self);
            return;
        }
        let (i, j) = (pos / n, pos % n);
        let mut free = !(self.rows[i] | self.cols[j]) & self.full;
        while free != 0 {
            let v = free.trailing_zeros() as u8;
            free &= free - 1;
            self.place(pos, v);
            self.walk(pos + 1, nodes, leaf);
            self.unplace(pos);
        }
    }
}

pub struct Model {
    pub table: Vec<u8>,
    pub aut: u32,
    pub signature: Signature,
}

pub struct Models {
    pub n: usize,
    pub models: Vec<Model>,
    pub nodes: u64,
    pub tasks: usize,
}

/// Every quasigroup of order `n` satisfying all of `required`, one per
/// isomorphism class, sorted by table, with full signatures over `engine`.
pub fn models(n: usize, required: &[Law], engine: &Engine) -> Models {
    let pt = PermTable::new(n);
    let laws: Vec<Compiled> = required.iter().map(Compiled::new).collect();
    let split = (2 * n).min(n * n);

    // Prefixes: a row-0 candidate followed by consistent cells up to `split`.
    let row0s = row0_candidates(&pt);
    let prefix_nodes = std::sync::atomic::AtomicU64::new(0);
    let prefixes: Vec<Vec<u8>> = row0s
        .par_iter()
        .flat_map_iter(|r0| {
            let mut s = State::new(n, &laws);
            for (j, &v) in r0.iter().enumerate() {
                s.place(j, v);
            }
            let mut out = Vec::new();
            let mut nodes = 0u64;
            if s.consistent() {
                collect(&mut s, n, split, &mut nodes, &mut out);
            }
            prefix_nodes.fetch_add(nodes, std::sync::atomic::Ordering::Relaxed);
            out
        })
        .collect();

    let (found, nodes): (Vec<(Vec<u8>, u32)>, u64) = prefixes
        .par_iter()
        .map(|prefix| {
            let mut s = State::new(n, &laws);
            for (pos, &v) in prefix.iter().enumerate() {
                s.place(pos, v);
            }
            let mut found = Vec::new();
            let mut nodes = 0u64;
            s.walk(prefix.len(), &mut nodes, &mut |st: &State| {
                let t = st.table();
                if let Some(aut) = canonical_aut(&t, &pt) {
                    found.push((t, aut));
                }
            });
            (found, nodes)
        })
        .reduce(
            || (Vec::new(), 0),
            |mut a, b| {
                a.0.extend(b.0);
                (a.0, a.1 + b.1)
            },
        );

    let mut models: Vec<Model> = found
        .into_par_iter()
        .map(|(table, aut)| {
            let signature = engine.signature(&Quasigroup::from_table(n, &table).unwrap());
            Model { table, aut, signature }
        })
        .collect();
    models.sort_by(|a, b| a.table.cmp(&b.table));
    Models {
        n,
        models,
        nodes: nodes + prefix_nodes.into_inner(),
        tasks: prefixes.len(),
    }
}

fn collect(s: &mut State, pos: usize, split: usize, nodes: &mut u64, out: &mut Vec<Vec<u8>>) {
    *nodes += 1;
    if !s.consistent() {
        return;
    }
    if pos == split {
        out.push(s.ops[0][..pos].to_vec());
        return;
    }
    let n = s.n;
    let (i, j) = (pos / n, pos % n);
    let mut free = !(s.rows[i] | s.cols[j]) & s.full;
    while free != 0 {
        let v = free.trailing_zeros() as u8;
        free &= free - 1;
        s.place(pos, v);
        collect(s, pos + 1, split, nodes, out);
        s.unplace(pos);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enumerate::iso_classes;

    /// With no laws, the search is the enumerator; with laws, it must agree
    /// with filtering the enumerator's output.
    #[test]
    fn agrees_with_enumeration_through_order_5() {
        let laws = crate::load_laws();
        let engine = Engine::new(&laws);
        for n in 1..=5 {
            let all = iso_classes(n);
            assert_eq!(models(n, &[], &engine).models.len(), all.len(), "order {n}, no laws");
            for id in [1usize, 5, 23, 38, 54, 102, 408] {
                let law = &laws[id - 1];
                let c = Compiled::new(law);
                let expect: Vec<Vec<u8>> = all
                    .iter()
                    .filter(|m| c.holds(&Quasigroup::from_table(n, &m.table).unwrap()))
                    .map(|m| m.table.clone())
                    .collect();
                let got: Vec<Vec<u8>> = models(n, std::slice::from_ref(law), &engine).models.into_iter().map(|m| m.table).collect();
                assert_eq!(got, expect, "order {n}, Sch-{id}");
            }
        }
    }
}
