//! Fast versions of the verification gates. The full-size runs (order 6 for
//! G1, order 5 for G4) are driven by the Makefile through the `qg` binary.

use quadratula::enumerate::{count_latin, factorial, iso_classes};
use quadratula::load_laws;
use quadratula::parastrophe::{S3, law_map};
use quadratula::qg::{Engine, Quasigroup};

/// OEIS A057991 (quasigroups up to isomorphism) and A002860 (Latin squares).
#[test]
fn g1_counts_through_order_5() {
    let iso = [1usize, 1, 5, 35, 1411];
    let latin = [1u64, 2, 12, 576, 161280];
    for n in 1..=5 {
        let classes = iso_classes(n);
        assert_eq!(classes.len(), iso[n - 1], "iso classes at order {n}");
        let orbit_sum: u64 = classes.iter().map(|c| factorial(n) / c.aut as u64).sum();
        assert_eq!(orbit_sum, latin[n - 1], "orbit-stabiliser sum at order {n}");
        assert_eq!(count_latin(n), latin[n - 1], "plain count at order {n}");
    }
}

#[test]
fn laws_parse_and_are_distinct() {
    let laws = load_laws();
    assert_eq!(laws.len(), 990);
    let keys: std::collections::HashSet<_> = laws.iter().map(|l| l.key()).collect();
    assert_eq!(keys.len(), 990, "two laws coincide up to renaming and side swap");
}

#[test]
fn g4_equivariance_through_order_3() {
    let laws = load_laws();
    let engine = Engine::new(&laws);
    for s in S3 {
        let m = law_map(&laws, s);
        for n in 1..=3 {
            for t in quadratula::enumerate::all_latin(n) {
                let q = Quasigroup::from_table(n, &t).unwrap();
                let base = engine.signature(&q);
                let conj = engine.signature(&q.conjugate(s));
                for (l, tl) in m.iter().enumerate() {
                    if let Some(tl) = tl {
                        assert_eq!(conj.get(l), base.get(*tl), "sigma {s:?} law {} order {n}", l + 1);
                    }
                }
            }
        }
    }
}
