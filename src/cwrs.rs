use entdec;

/// 4.3.4.2. PVQ Decoding
/// 
/// Decoding of PVQ vectors is implemented in decode_pulses(). The unique codeword index is decoded as a uniformly distributed integer value between 0 and V(N, K) - 1, where V(N, K) is the number of possible combinations of K pulses in N samples. The index is then converted to a vector in the same way specified in [PVQ]. The indexing is based on the calculation of V(N, K) (denoted N(L, K) in [PVQ]).
/// 
/// The number of combinations can be computed recursively as V(N, K) = V(N - 1, K) + V(N, K - 1) + V(N - 1, K - 1), with V(N, 0) = 1 and V(0, K) = 0, K ≠ 0. There are many different ways to compute V(N, K), including precomputed tables and direct use of the recursive formulation. The reference implementation applies the recursive formulation one line (or column) at a time to save on memory use, along with an alternate, univariate recurrence to initialize an arbitrary line, and direct polynomial solutions for small N. All of these methods are equivalent, and have different trade-offs in speed, memory usage, and code size. Implementations MAY use any methods they like, as long as they are equivalent to the mathematical definition.
/// 
/// The decoded vector X is recovered as follows. Let i be the index decoded with the procedure in Section 4.1.5 with ft = V(N, K), so that 0 <= i < V(N, K). Let k = K. Then, for j ∈ [0, N - 1], do:
/// 
/// 1. Let p = (V(N - j - 1, k) + V(N - j, k)) / 2.
/// 2. If i < p, then let sgn = 1, else let sgn = -1 and set i = i - p.
/// 3. Let k0 = k and set p = p - V(N - j - 1, k).
/// 4. While p > i, set k = k - 1 and p = p - V(N - j - 1, k).
/// 5. Set X[j] = sgn * (k0 - k) and i = i - p.
/// 
/// The decoded vector X is then normalized such that its L2-norm equals one.
///
/// * [PVQ] Fischer, T., "A Pyramid Vector Quantizer", _IEEE Transactions on Information Theory,_ Vol. 32, pp. 568-583, July 1986.
pub fn decode_pulses(x: &mut [f32],
                     mut k: usize,
                     ec: &mut entdec::EntropyCoder,
                     v: &Vec<Vec<Option<u32>>>) {
    let n = x.len();
    let mut i = ec.decode_uint(v[n][k].unwrap());
    for j in 0..n {
        let mut p = ((v[n - j][k].unwrap() as u64 + v[n - j - 1][k].unwrap() as u64) / 2) as u32;
        let sgn = if i < p {
            1
        } else {
            i -= p;
            -1
        };
        let k0 = k;
        p -= v[n - j - 1][k].unwrap();
        while p > i {
            k -= 1;
            p -= v[n - j - 1][k].unwrap();
        }
        x[j] = sgn as f32 * (k0 - k) as f32;
        i -= p;
    }
}
