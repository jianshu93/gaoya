// fast_sim_hash.rs -----------------------------------------------------------
/// Charikar SimHash accelerated with the bit-hack from “Speeding up SimHash by 10× using a bit hack”
/// see webpage here: https://www.dynatrace.com/engineering/blog/speeding-up-simhash-by-10x-using-a-bit-hack/

use core::marker::PhantomData;
use std::hash::Hash;
use wyhash::WyRng;
use rand_core::{RngCore, SeedableRng};
use crate::simhash::SimHashBits;
use crate::simhash::sim_hasher::SimHasher;
use crate::simhash::BitArray;

pub struct FastSimHash<H, S, const L: usize, const BULK: usize = 3>
where
    H: SimHasher<T = u64>,        // we seed the PRNG with a 64-bit hash
    S: SimHashBits,               // u64  → up to 64 components
                                  // u128 → up to 128 components
{
    hasher: H,
    marker: PhantomData<S>,
}

impl<H, S, const L: usize, const BULK: usize> FastSimHash<H, S, L, BULK>
where
    H: SimHasher<T = u64>,
    S: SimHashBits,
{

    // 2^BULK counters are updated in one go
    const PER_LANE: usize = 1 << BULK;               // 8
    // Each counter uses 2^(6-BULK) bits inside a u64 “lane”
    const BITS_PER_COUNTER: usize = 1 << (6 - BULK); // 8
    // Bit-mask with the low BITS_PER_COUNTER bits set
    const BULK_MASK: u64 = Self::calc_bulk_mask();
    /// Max value a temporary counter can take before flushing to `counts`
    const TMP_LIMIT: u64 = Self::calc_tmp_limit();

    // Pre-compute repeating bitmask  …001001001…
    const fn calc_bulk_mask() -> u64 {
        let mut mask = 1u64;
        let mut i = 0;
        while i < BULK {
            // ore each 1<<(1<<(5-i)) pattern
            mask |= mask << (1 << (5 - i));
            i += 1;
        }
        mask
    }

    // calculateTemporaryCounterLimit
    const fn calc_tmp_limit() -> u64 {
        (1u64 << 1 << ((1 << (6 - BULK)) - 1)) - 1
    }

    pub fn new(hasher: H) -> Self {
        assert!(L <= S::bit_length(), "Signature length too large for SimHashBits");
        Self { hasher, marker: PhantomData }
    }

    /// Produce the L-bit SimHash signature for an iterator of features.
    /// * `T` is any iterator
    /// * `U: Hash` is the feature type (tokens, shingles, integers, …)
    pub fn create_signature<T, U>(&self, iter: T) -> S
    where
        T: Iterator<Item = U>,
        U: Hash,
    {
        // permanent & temporary counters
        let mut counts: [i32; L] = [0; L];
        let tmp_len = (L + (63 >> (6 - BULK))) >> BULK;   // identical to Java
        let mut tmp = vec![0u64; tmp_len];

        // pre-compute constants used below
        let num_chunks = tmp_len >> (6 - BULK);        // full 64-bit groups
        let num_tail_lanes = tmp_len & (0x3f >> BULK);
        let mut processed = 0u64;

        for feature in iter {
            // hash feature → seed PRNG (WyRand = one 64-bit state, two ops per draw)
            let seed = self.hasher.hash(&feature);
            let mut prng = WyRng::seed_from_u64(seed);

            // update tmp counters in 64-bit chunks
            for h in 0..num_chunks {
                let rnd = prng.next_u64();
                let off = h << (6 - BULK);                    // multiply by lanes-per-chunk
                for j in 0..Self::BITS_PER_COUNTER {         // 8 lanes inside `rnd`
                    tmp[off + j] += (rnd >> j) & Self::BULK_MASK;
                }
            }
            if num_tail_lanes > 0 {
                let rnd = prng.next_u64();
                let off = num_chunks << (6 - BULK);
                for j in 0..num_tail_lanes {
                    tmp[off + j] += (rnd >> j) & Self::BULK_MASK;
                }
            }

            processed += 1;
            if processed == Self::TMP_LIMIT {                 // flush before overflow
                Self::flush_tmp::<L, BULK>(&mut counts, &mut tmp[..]);
                processed = 0;
            }
        }
        // final flush (handles empty iter gracefully)
        Self::flush_tmp::<L, BULK>(&mut counts, &mut tmp[..]);

        // threshold to signature 
        let n = processed as i32; // number of elements in last (possibly empty) block
        let total = (counts[0] as i64) + (n as i64); // only used to silence dead-code warning

        let num_elements = counts
            .iter()
            .map(|&x| if x < 0 { -x } else { x })
            .sum::<i32>() as u32; // not strictly necessary, just sanity

        let limit = (num_elements >> 1) as i32;
        let mut signature = S::zero();

        for i in 0..L {
            // tie-breaker identical to the Java reference
            if counts[i] + ((i & ((!num_elements) as usize & 1)) as i32) > limit {
                signature |= S::one() << i;
            }
        }
        signature
    }

    /// Merge the temporary packed counters into the final `counts` array.
    fn flush_tmp<const L2: usize, const B: usize>(
        counts: &mut [i32; L2],
        tmp: &mut [u64],
    ) {
        // full BULK-sized bundles
        let full = counts.len() >> B;
        for h in 0..full {
            let t = core::mem::take(&mut tmp[h]);
            let off = h << B;
            for g in 0..(1 << B) {
                counts[off + g] +=
                    ((t >> (g << (6 - B))) & Self::TMP_LIMIT) as i32;
            }
        }
        // tail (< B counters)
        for h in full..tmp.len() {
            let t  = core::mem::take(&mut tmp[h]);
            let off = h << B;
            for g in 0..(counts.len() - off) {
                counts[off + g] +=
                    ((t >> (g << (6 - B))) & Self::TMP_LIMIT) as i32;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::FastSimHash;
    use crate::simhash::BitArray;
    use crate::simhash::SimHashBits;
    use num_traits::real::Real;
    use super::SimHasher;
    use crate::simhash::sim_hasher::{Xxh3Hasher64, Xxh3Hasher128};
    use std::time::Instant;
    fn whitespace_split(s: &str) -> impl Iterator<Item=&str> { s.split_whitespace() }

    static S1: &str = "SimHash is a technique used for detecting near-duplicates ...";
    static S2: &str = "SimHash is a technique used for detecting near-duplicates ... utilized ...";

    #[test]
    fn fast_simhash() {
        type Bits = u128;
        const L: usize = 128;

        let fsh = FastSimHash::<Xxh3Hasher64, Bits, L>::new(Xxh3Hasher64::new());
        let h1 = fsh.create_signature(whitespace_split(S1));
        let h2 = fsh.create_signature(whitespace_split(S2));

        // Should differ by roughly 12 % of the bits
        assert!(h1.hamming_distance(&h2) < 16);
    }
    #[test]
    fn fast_simhash_large() {
        use rand::{Rng, SeedableRng};
        use rand::rngs::StdRng;

        // Helper identical to the small test
        fn whitespace_split(s: &str) -> impl Iterator<Item = &str> { s.split_whitespace() }

        type Bits = u128;
        const L: usize = 128;
        // build deterministic test data 
        const N: usize = 10_000;
        let mut rng = StdRng::seed_from_u64(42);

        // S1: 10 000 random u64 numbers separated by spaces
        let data1: Vec<u64> = (0..N).map(|_| rng.gen()).collect();
        let s1 = data1.iter().map(u64::to_string).collect::<Vec<_>>().join(" ");

        // S2: clone + tweak every 20th element  (≈5 % difference)
        let mut data2 = data1.clone();
        for i in (0..N).step_by(20) {
            data2[i] = data2[i].wrapping_add(1);
        }
        let s2 = data2.iter().map(u64::to_string).collect::<Vec<_>>().join(" ");

        let fsh = FastSimHash::<Xxh3Hasher64, Bits, L>::new(Xxh3Hasher64::new());
        let t1 = Instant::now();
        let h1 = fsh.create_signature(whitespace_split(&s1));
        let h2 = fsh.create_signature(whitespace_split(&s2));
        let dur_fast = t1.elapsed();
        println!("fast  SimHash: {:?}", dur_fast);
        // 5 % token change → expect roughly 5 % differing bits (≈6 of 128)
        assert!(h1.hamming_distance(&h2) < 6);
    }
    #[test]
    fn fast_simhash_bitarray() {
        use rand::{Rng, SeedableRng};
        use rand::rngs::StdRng;

        // Helper identical to the small test
        fn whitespace_split(s: &str) -> impl Iterator<Item = &str> { s.split_whitespace() }

        type Bits = BitArray<2>;
        const L: usize = 128;
        // build deterministic test data 
        const N: usize = 10_000;
        let mut rng = StdRng::seed_from_u64(42);

        // S1: 10 000 random u64 numbers separated by spaces
        let data1: Vec<u64> = (0..N).map(|_| rng.gen()).collect();
        //let s1 = data1.iter().map(u64::to_string).collect::<Vec<_>>().join(" ");
        let s1 = &data1;
        // S2: clone + tweak every 20th element  (≈5 % difference)
        let mut data2 = data1.clone();
        for i in (0..N).step_by(20) {
            data2[i] = data2[i].wrapping_add(1);
        }
        //let s2 = data2.iter().map(u64::to_string).collect::<Vec<_>>().join(" ");
        let s2 = &data2;
        let fsh = FastSimHash::<Xxh3Hasher64, Bits, L, 4>::new(Xxh3Hasher64::new());
        let t1 = Instant::now();
        let h1 = fsh.create_signature(s1.iter().cloned());
        let h2 = fsh.create_signature(s2.iter().cloned());
        let dur_fast = t1.elapsed();
        println!("fast  SimHash: {:?}", dur_fast);
        // 5 % token change → expect roughly 5 % differing bits
        assert!(h1.hamming_distance(&h2) < 17);
    }
}