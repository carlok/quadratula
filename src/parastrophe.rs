//! The S3 action on quasigroups (parastrophes) and the induced action on laws.
//!
//! A quasigroup is a set of triples `(t0, t1, t2)` with `t0 * t1 = t2`, one
//! per pair of positions. Le Floch's three operations each recover one
//! position from the other two, read cyclically:
//!
//! | output position `k` | operation | arguments       |
//! |---------------------|-----------|-----------------|
//! | 2                   | `*`       | `(t0, t1)`      |
//! | 0                   | `//`      | `(t1, t2)`      |
//! | 1                   | `\\`      | `(t2, t0)`      |
//!
//! so `op_k(t_{k+1}, t_{k+2}) = t_k` (indices mod 3). The parastrophe for `s`
//! has triples `u` with `u[s[k]] = t[k]`. Its operation with output position
//! `m` is `op_m^s(u_{m+1}, u_{m+2}) = u_m = t_k` where `s[k] = m`, which is
//! `op_k` applied to `(u_{s[k+1]}, u_{s[k+2]})`, i.e. to the same two
//! arguments, swapped unless `s[k+1] = m+1`. Rewriting every operation of a
//! law this way gives a law `tau_s(L)` over the original operations with
//! `Q^s |= L  iff  Q |= tau_s(L)`.

use crate::law::{Law, Op};
use std::collections::HashMap;

/// The six permutations of triple positions. Index 0 is the identity.
pub const S3: [[usize; 3]; 6] = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [1, 0, 2], [0, 2, 1], [2, 1, 0]];

pub fn name(s: [usize; 3]) -> String {
    s.iter().map(|d| char::from_digit(*d as u32, 10).unwrap()).collect()
}

fn output_position(o: Op) -> usize {
    match o {
        Op::Star => 2,
        Op::Slash => 0,
        Op::Back => 1,
    }
}

fn op_at(k: usize) -> Op {
    [Op::Slash, Op::Back, Op::Star][k]
}

/// How an operation of `Q^s` is written with the operations of `Q`.
pub fn translate_op(s: [usize; 3], o: Op) -> (Op, bool) {
    let m = output_position(o);
    let k = (0..3).find(|&k| s[k] == m).unwrap();
    let swap = s[(k + 1) % 3] != (m + 1) % 3;
    (op_at(k), swap)
}

/// `tau_s(L)`.
pub fn translate_law(s: [usize; 3], law: &Law) -> Law {
    let g = |o: Op| translate_op(s, o);
    Law { id: 0, lhs: law.lhs.map_ops(&g), rhs: law.rhs.map_ops(&g) }
}

/// For each law (0-based index), the 0-based index of `tau_s(L)` in the list,
/// or `None` if the translated law is not in the list up to renaming and
/// swapping sides.
pub fn law_map(laws: &[Law], s: [usize; 3]) -> Vec<Option<usize>> {
    let index: HashMap<String, usize> = laws.iter().enumerate().map(|(i, l)| (l.key(), i)).collect();
    laws.iter().map(|l| index.get(&translate_law(s, l).key()).copied()).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identity_translates_to_itself() {
        for o in Op::ALL {
            assert_eq!(translate_op(S3[0], o), (o, false));
        }
    }

    #[test]
    fn transpose_swaps_star_arguments() {
        assert_eq!(translate_op([1, 0, 2], Op::Star), (Op::Star, true));
    }
}
