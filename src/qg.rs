//! Quasigroups as Latin squares, their three operations, and law evaluation.

use crate::law::{Law, Op, Term};
use crate::sig::Signature;

/// A finite quasigroup on `{0..n-1}` with all three of Le Floch's operations
/// tabulated row-major: `ops[o][x * n + y] = x o y`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Quasigroup {
    pub n: usize,
    pub ops: [Vec<u8>; 3],
}

impl Quasigroup {
    /// Build from a multiplication table, rejecting anything that is not a
    /// Latin square.
    pub fn from_table(n: usize, table: &[u8]) -> Result<Self, String> {
        if table.len() != n * n {
            return Err(format!("table has {} cells, expected {}", table.len(), n * n));
        }
        let unset = u8::MAX;
        let mut slash = vec![unset; n * n];
        let mut back = vec![unset; n * n];
        for c in 0..n {
            for a in 0..n {
                let b = table[c * n + a] as usize;
                if b >= n {
                    return Err(format!("entry {b} out of range"));
                }
                // a // b = c  iff  c * a = b
                if slash[a * n + b] != unset {
                    return Err(format!("column {a} repeats value {b}"));
                }
                slash[a * n + b] = c as u8;
            }
        }
        for b in 0..n {
            for c in 0..n {
                let a = table[b * n + c] as usize;
                // a \\ b = c  iff  b * c = a
                if back[a * n + b] != unset {
                    return Err(format!("row {b} repeats value {a}"));
                }
                back[a * n + b] = c as u8;
            }
        }
        Ok(Quasigroup { n, ops: [table.to_vec(), slash, back] })
    }

    #[inline]
    pub fn op(&self, o: Op, x: u8, y: u8) -> u8 {
        self.ops[o.index()][x as usize * self.n + y as usize]
    }

    pub fn table(&self) -> &[u8] {
        &self.ops[0]
    }

    /// The parastrophe for a permutation `s` of triple positions: the
    /// multiplication triples `(t0, t1, t2)` with `t0 * t1 = t2` are sent to
    /// `u` with `u[s[k]] = t[k]`.
    pub fn conjugate(&self, s: [usize; 3]) -> Quasigroup {
        let n = self.n;
        let mut t = vec![0u8; n * n];
        for x in 0..n {
            for y in 0..n {
                let tri = [x as u8, y as u8, self.ops[0][x * n + y]];
                let mut u = [0u8; 3];
                for k in 0..3 {
                    u[s[k]] = tri[k];
                }
                t[u[0] as usize * n + u[1] as usize] = u[2];
            }
        }
        Quasigroup::from_table(n, &t).expect("a parastrophe of a Latin square is Latin")
    }

    pub fn table_string(&self) -> String {
        self.ops[0].iter().map(|&v| char::from_digit(v as u32, 36).unwrap()).collect()
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum Ins {
    Var(u8),
    Op(u8),
}

/// A law compiled to two postfix programs over variable slots `0..k`.
#[derive(Clone, Debug)]
pub struct Compiled {
    pub(crate) k: usize,
    pub(crate) lhs: Vec<Ins>,
    pub(crate) rhs: Vec<Ins>,
}

pub(crate) const STACK: usize = 16;
pub(crate) const MAX_VARS: usize = 8;

impl Compiled {
    pub fn new(law: &Law) -> Compiled {
        let vars = law.vars();
        assert!(vars.len() <= MAX_VARS, "too many variables in {law}");
        fn emit(t: &Term, vars: &[u8], out: &mut Vec<Ins>) {
            match t {
                Term::Var(v) => out.push(Ins::Var(vars.iter().position(|w| w == v).unwrap() as u8)),
                Term::App(o, l, r) => {
                    emit(l, vars, out);
                    emit(r, vars, out);
                    out.push(Ins::Op(o.index() as u8));
                }
            }
        }
        let mut lhs = Vec::new();
        let mut rhs = Vec::new();
        emit(&law.lhs, &vars, &mut lhs);
        emit(&law.rhs, &vars, &mut rhs);
        assert!(law.lhs.depth() < STACK && law.rhs.depth() < STACK);
        Compiled { k: vars.len(), lhs, rhs }
    }

    #[inline]
    fn eval(prog: &[Ins], env: &[u8; MAX_VARS], q: &Quasigroup) -> u8 {
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
                    sp -= 1;
                    st[sp - 1] = q.ops[o as usize][a * q.n + b];
                }
            }
        }
        st[0]
    }

    /// Does `q` satisfy the law for every assignment of its variables?
    pub fn holds(&self, q: &Quasigroup) -> bool {
        let n = q.n as u8;
        let mut env = [0u8; MAX_VARS];
        loop {
            if Self::eval(&self.lhs, &env, q) != Self::eval(&self.rhs, &env, q) {
                return false;
            }
            let mut i = 0;
            loop {
                if i == self.k {
                    return true;
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
}

/// Signature engine over a fixed law list.
pub struct Engine {
    laws: Vec<Compiled>,
}

impl Engine {
    pub fn new(laws: &[Law]) -> Engine {
        Engine { laws: laws.iter().map(Compiled::new).collect() }
    }

    pub fn n_laws(&self) -> usize {
        self.laws.len()
    }

    /// Bit `i` is set iff `q` satisfies the law on line `i + 1`.
    pub fn signature(&self, q: &Quasigroup) -> Signature {
        let mut s = Signature::zeros(self.laws.len());
        for (i, c) in self.laws.iter().enumerate() {
            if c.holds(q) {
                s.set(i);
            }
        }
        s
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn z3() -> Quasigroup {
        // x * y = x + y mod 3
        let t: Vec<u8> = (0..9).map(|i| ((i / 3 + i % 3) % 3) as u8).collect();
        Quasigroup::from_table(3, &t).unwrap()
    }

    #[test]
    fn rejects_non_latin() {
        assert!(Quasigroup::from_table(2, &[0, 0, 1, 1]).is_err());
        assert!(Quasigroup::from_table(2, &[0, 1, 0, 1]).is_err());
    }

    /// Le Floch's quasigroup axioms (5), exhaustively on a non-commutative
    /// quasigroup of order 3 (x * y = 2x + y mod 3).
    #[test]
    fn le_floch_axioms_hold() {
        let t: Vec<u8> = (0..9).map(|i| ((2 * (i / 3) + i % 3) % 3) as u8).collect();
        let q = Quasigroup::from_table(3, &t).unwrap();
        let (s, sl, bk) = (Op::Star, Op::Slash, Op::Back);
        for x in 0..3 {
            for y in 0..3 {
                assert_eq!(q.op(s, q.op(sl, y, x), y), x);
                assert_eq!(q.op(sl, y, q.op(s, x, y)), x);
                assert_eq!(q.op(s, x, q.op(bk, y, x)), y);
                assert_eq!(q.op(bk, q.op(s, x, y), x), y);
                assert_eq!(q.op(bk, x, q.op(sl, y, x)), y);
                assert_eq!(q.op(sl, q.op(bk, x, y), x), y);
            }
        }
    }

    #[test]
    fn transpose_is_opposite() {
        let t: Vec<u8> = (0..9).map(|i| ((2 * (i / 3) + i % 3) % 3) as u8).collect();
        let q = Quasigroup::from_table(3, &t).unwrap();
        let op = q.conjugate([1, 0, 2]);
        for x in 0..3u8 {
            for y in 0..3u8 {
                assert_eq!(op.op(Op::Star, x, y), q.op(Op::Star, y, x));
            }
        }
    }

    #[test]
    fn evaluates_simple_laws() {
        let q = z3();
        let law = |s: &str| {
            let (l, r) = crate::law::parse_law_body(s).unwrap();
            Compiled::new(&Law { id: 0, lhs: l, rhs: r })
        };
        assert!(law("x * (y * z) = (x * y) * z").holds(&q));
        assert!(law("x * (y * z) = x * (z * y)").holds(&q));
        assert!(!law("x * (y * z) = x // (y * z)").holds(&q));
    }
}
