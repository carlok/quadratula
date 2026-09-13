//! Parsing of Schröder's 990 quasigroup laws in Le Floch's syntax
//! (`data/lefloch/quasigroup-equations.txt`, arXiv:2603.29909, CC BY 4.0).
//!
//! Each line is `Sch-N: lhs = rhs`. Terms are built from variables (single
//! lowercase letters) and three infix operations written `*`, `//` and `\\`
//! (two backslash characters). `//` and `\\` are Le Floch's *conjugate*
//! divisions: `x // y = y / x` and `x \\ y = y \ x`, where `/` and `\` are the
//! usual right and left divisions. Nothing about which operations occur is
//! assumed here: the parser accepts exactly these three tokens and rejects
//! anything else.
//!
//! Grammar, deliberately strict: `expr := atom (op atom)?`, `atom := var |
//! '(' expr ')'`. A side with two operations and no parentheses (`a * b * c`)
//! is rejected rather than read under some associativity convention. This
//! follows parsimagma's `law.rs` rule of rejecting rather than guessing.

use std::fmt;

/// The three quasigroup operations, in Le Floch's cyclic order.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, PartialOrd, Ord)]
pub enum Op {
    /// `x * y`
    Star,
    /// `x // y`, the `z` with `z * x = y`
    Slash,
    /// `x \\ y`, the `z` with `y * z = x`
    Back,
}

impl Op {
    pub const ALL: [Op; 3] = [Op::Star, Op::Slash, Op::Back];

    pub fn token(self) -> &'static str {
        match self {
            Op::Star => "*",
            Op::Slash => "//",
            Op::Back => "\\\\",
        }
    }

    pub fn index(self) -> usize {
        self as usize
    }

    pub fn from_index(i: usize) -> Op {
        Op::ALL[i]
    }
}

/// A term. Variables are stored by the letter's offset from `a`, so the
/// original names survive and printing round-trips.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub enum Term {
    Var(u8),
    App(Op, Box<Term>, Box<Term>),
}

impl Term {
    pub fn ops(&self) -> usize {
        match self {
            Term::Var(_) => 0,
            Term::App(_, l, r) => 1 + l.ops() + r.ops(),
        }
    }

    pub fn depth(&self) -> usize {
        match self {
            Term::Var(_) => 0,
            Term::App(_, l, r) => 1 + l.depth().max(r.depth()),
        }
    }

    pub fn visit_ops(&self, f: &mut impl FnMut(Op)) {
        if let Term::App(o, l, r) = self {
            f(*o);
            l.visit_ops(f);
            r.visit_ops(f);
        }
    }

    pub fn visit_vars(&self, f: &mut impl FnMut(u8)) {
        match self {
            Term::Var(v) => f(*v),
            Term::App(_, l, r) => {
                l.visit_vars(f);
                r.visit_vars(f);
            }
        }
    }

    /// Apply `f` to every variable.
    pub fn map_vars(&self, f: &impl Fn(u8) -> u8) -> Term {
        match self {
            Term::Var(v) => Term::Var(f(*v)),
            Term::App(o, l, r) => Term::App(*o, Box::new(l.map_vars(f)), Box::new(r.map_vars(f))),
        }
    }

    /// Replace every application `o(l, r)` by `g(o)`, where `g` returns the
    /// new operation and whether to swap the arguments.
    pub fn map_ops(&self, g: &impl Fn(Op) -> (Op, bool)) -> Term {
        match self {
            Term::Var(v) => Term::Var(*v),
            Term::App(o, l, r) => {
                let (o2, swap) = g(*o);
                let (l2, r2) = (l.map_ops(g), r.map_ops(g));
                if swap {
                    Term::App(o2, Box::new(r2), Box::new(l2))
                } else {
                    Term::App(o2, Box::new(l2), Box::new(r2))
                }
            }
        }
    }
}

impl fmt::Display for Term {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        fn go(t: &Term, top: bool, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            match t {
                Term::Var(v) => write!(f, "{}", (b'a' + v) as char),
                Term::App(o, l, r) => {
                    if !top {
                        write!(f, "(")?;
                    }
                    go(l, false, f)?;
                    write!(f, " {} ", o.token())?;
                    go(r, false, f)?;
                    if !top {
                        write!(f, ")")?;
                    }
                    Ok(())
                }
            }
        }
        go(self, true, f)
    }
}

#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct Law {
    /// 1-based Schröder number as printed in the file (`Sch-N`).
    pub id: u32,
    pub lhs: Term,
    pub rhs: Term,
}

impl fmt::Display for Law {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} = {}", self.lhs, self.rhs)
    }
}

impl Law {
    /// Distinct variables, in order of first appearance (lhs then rhs).
    pub fn vars(&self) -> Vec<u8> {
        let mut out = Vec::new();
        let mut push = |v: u8| {
            if !out.contains(&v) {
                out.push(v)
            }
        };
        self.lhs.visit_vars(&mut push);
        self.rhs.visit_vars(&mut push);
        out
    }

    /// A key identifying the law up to renaming variables and swapping the
    /// two sides: rename variables to `a, b, c, ...` by first appearance and
    /// take the smaller printed form of the two orientations.
    pub fn key(&self) -> String {
        let orient = |l: &Term, r: &Term| {
            let tmp = Law { id: 0, lhs: l.clone(), rhs: r.clone() };
            let vs = tmp.vars();
            let ren = |v: u8| vs.iter().position(|&w| w == v).unwrap() as u8;
            format!("{} = {}", l.map_vars(&ren), r.map_vars(&ren))
        };
        let a = orient(&self.lhs, &self.rhs);
        let b = orient(&self.rhs, &self.lhs);
        a.min(b)
    }

    pub fn is_tautology(&self) -> bool {
        self.lhs == self.rhs
    }
}

#[derive(Debug)]
pub struct ParseError {
    pub line: usize,
    pub msg: String,
}

impl fmt::Display for ParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "line {}: {}", self.line, self.msg)
    }
}

impl std::error::Error for ParseError {}

/// Parse the whole law file. Line `k` (1-based) must be `Sch-k: ...`.
pub fn parse_laws(text: &str) -> Result<Vec<Law>, ParseError> {
    let mut out = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        let err = |msg: String| ParseError { line: i + 1, msg };
        let (label, body) = line
            .split_once(':')
            .ok_or_else(|| err("missing `Sch-N:` label".into()))?;
        let id: u32 = label
            .strip_prefix("Sch-")
            .and_then(|s| s.parse().ok())
            .ok_or_else(|| err(format!("bad label {label:?}")))?;
        if id as usize != out.len() + 1 {
            return Err(err(format!("expected Sch-{}, found Sch-{id}", out.len() + 1)));
        }
        let (lhs, rhs) = parse_law_body(body).map_err(err)?;
        out.push(Law { id, lhs, rhs });
    }
    Ok(out)
}

/// Parse `lhs = rhs`.
pub fn parse_law_body(body: &str) -> Result<(Term, Term), String> {
    let toks = tokenize(body)?;
    let mut p = Parser { toks: &toks, pos: 0 };
    let lhs = p.expr()?;
    p.expect(&Tok::Eq)?;
    let rhs = p.expr()?;
    if p.pos != toks.len() {
        return Err(format!("trailing tokens after position {}", p.pos));
    }
    Ok((lhs, rhs))
}

#[derive(Debug, PartialEq)]
enum Tok {
    Var(u8),
    Op(Op),
    LParen,
    RParen,
    Eq,
}

fn tokenize(s: &str) -> Result<Vec<Tok>, String> {
    let b = s.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        match b[i] {
            b' ' | b'\t' => i += 1,
            b'(' => {
                out.push(Tok::LParen);
                i += 1
            }
            b')' => {
                out.push(Tok::RParen);
                i += 1
            }
            b'=' => {
                out.push(Tok::Eq);
                i += 1
            }
            b'*' => {
                out.push(Tok::Op(Op::Star));
                i += 1
            }
            b'/' if b.get(i + 1) == Some(&b'/') => {
                out.push(Tok::Op(Op::Slash));
                i += 2
            }
            b'\\' if b.get(i + 1) == Some(&b'\\') => {
                out.push(Tok::Op(Op::Back));
                i += 2
            }
            c @ b'a'..=b'z' => {
                if b.get(i + 1).is_some_and(|d| d.is_ascii_alphanumeric()) {
                    return Err(format!("multi-character identifier at byte {i}"));
                }
                out.push(Tok::Var(c - b'a'));
                i += 1
            }
            c => return Err(format!("unexpected character {:?} at byte {i}", c as char)),
        }
    }
    Ok(out)
}

struct Parser<'a> {
    toks: &'a [Tok],
    pos: usize,
}

impl Parser<'_> {
    fn expect(&mut self, t: &Tok) -> Result<(), String> {
        if self.toks.get(self.pos) == Some(t) {
            self.pos += 1;
            Ok(())
        } else {
            Err(format!("expected {t:?} at token {}", self.pos))
        }
    }

    fn expr(&mut self) -> Result<Term, String> {
        let l = self.atom()?;
        if let Some(Tok::Op(o)) = self.toks.get(self.pos) {
            let o = *o;
            self.pos += 1;
            let r = self.atom()?;
            if let Some(Tok::Op(_)) = self.toks.get(self.pos) {
                return Err(format!("unparenthesised operation chain at token {}", self.pos));
            }
            Ok(Term::App(o, Box::new(l), Box::new(r)))
        } else {
            Ok(l)
        }
    }

    fn atom(&mut self) -> Result<Term, String> {
        match self.toks.get(self.pos) {
            Some(Tok::Var(v)) => {
                self.pos += 1;
                Ok(Term::Var(*v))
            }
            Some(Tok::LParen) => {
                self.pos += 1;
                let e = self.expr()?;
                self.expect(&Tok::RParen)?;
                Ok(e)
            }
            other => Err(format!("expected variable or '(' at token {}, found {other:?}", self.pos)),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_all_three_ops() {
        let (l, r) = parse_law_body("x // (y \\\\ z) = (z * y) // x").unwrap();
        assert_eq!(l.to_string(), "x // (y \\\\ z)");
        assert_eq!(r.to_string(), "(z * y) // x");
    }

    #[test]
    fn rejects_chains_and_junk() {
        assert!(parse_law_body("x * y * z = x").is_err());
        assert!(parse_law_body("x / y = x").is_err());
        assert!(parse_law_body("x * y = xy").is_err());
    }

    #[test]
    fn key_is_invariant_under_renaming_and_swap() {
        let a = parse_law_body("x * (y * z) = y * (x * z)").unwrap();
        let b = parse_law_body("q * (p * r) = p * (q * r)").unwrap();
        let la = Law { id: 0, lhs: a.0, rhs: a.1 };
        let lb = Law { id: 0, lhs: b.1, rhs: b.0 };
        assert_eq!(la.key(), lb.key());
    }
}
