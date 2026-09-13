//! Law signatures: the bit vector recording which laws a quasigroup satisfies.
//!
//! Vendored from parsimagma (`src/sig.rs` at commit 2a75a1b,
//! https://github.com/carlok/parsimagma, Apache-2.0, Copyright 2026 Carlo
//! Perassi). Changes: documentation only; the width is a parameter, here the
//! 990 Schröder laws.
//!
//! Every separation question reduces to a bit test. `Q` separates law `i`
//! from law `j` exactly when `sig[i] && !sig[j]`.

use std::fmt;

/// A packed bit vector over a fixed number of laws.
#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Signature {
    nbits: usize,
    words: Vec<u64>,
}

impl Signature {
    pub fn zeros(nbits: usize) -> Self {
        Signature {
            nbits,
            words: vec![0u64; nbits.div_ceil(64)],
        }
    }

    #[inline]
    pub fn nbits(&self) -> usize {
        self.nbits
    }

    #[inline]
    pub fn words(&self) -> &[u64] {
        &self.words
    }

    #[inline]
    pub fn set(&mut self, i: usize) {
        debug_assert!(i < self.nbits);
        self.words[i >> 6] |= 1u64 << (i & 63);
    }

    #[inline]
    pub fn clear(&mut self, i: usize) {
        debug_assert!(i < self.nbits);
        self.words[i >> 6] &= !(1u64 << (i & 63));
    }

    #[inline]
    pub fn get(&self, i: usize) -> bool {
        debug_assert!(i < self.nbits);
        self.words[i >> 6] >> (i & 63) & 1 == 1
    }

    pub fn count(&self) -> u32 {
        self.words.iter().map(|w| w.count_ones()).sum()
    }

    /// Indices of the set bits, ascending.
    pub fn iter_set(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(wi, &w)| {
            let mut w = w;
            std::iter::from_fn(move || {
                if w == 0 {
                    None
                } else {
                    let b = w.trailing_zeros() as usize;
                    w &= w - 1;
                    Some(wi * 64 + b)
                }
            })
        })
    }

    /// 1-based ETP equation ids of the satisfied laws.
    pub fn satisfied_ids(&self) -> Vec<u32> {
        self.iter_set().map(|i| i as u32 + 1).collect()
    }

    /// Does this magma separate some pair? True when it satisfies at least
    /// one law and refutes at least one.
    pub fn is_separating(&self) -> bool {
        let c = self.count() as usize;
        c > 0 && c < self.nbits
    }

    /// Number of ordered pairs `(i, j)` this signature discharges, i.e.
    /// `|{i : sig[i]}| * |{j : !sig[j]}|`, excluding the diagonal (which
    /// cannot occur: a law is never both satisfied and refuted).
    pub fn separations(&self) -> u64 {
        let sat = self.count() as u64;
        sat * (self.nbits as u64 - sat)
    }

    /// Little-endian byte serialisation, `ceil(nbits/8)` bytes.
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.nbits.div_ceil(8));
        for (i, &w) in self.words.iter().enumerate() {
            let b = w.to_le_bytes();
            let take = (self.nbits.div_ceil(8) - i * 8).min(8);
            out.extend_from_slice(&b[..take]);
        }
        out
    }

    pub fn from_bytes(nbits: usize, bytes: &[u8]) -> Self {
        assert_eq!(bytes.len(), nbits.div_ceil(8), "wrong signature length");
        let mut s = Signature::zeros(nbits);
        for (i, &b) in bytes.iter().enumerate() {
            s.words[i / 8] |= (b as u64) << ((i % 8) * 8);
        }
        s
    }
}

impl fmt::Debug for Signature {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Signature({}/{} laws)", self.count(), self.nbits)
    }
}
