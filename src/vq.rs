use cwrs;
use entdec;
use std;
use utils;

/// 4.3.4.3. Spreading
///
/// The normalized vector decoded in Section 4.3.4.2 is then rotated for the purpose of avoiding tonal artifacts. The rotation gain is equal to
///
///     g_r = N / (N + f_r * K)
///
/// where N is the number of dimensions, K is the number of pulses, and f_r depends on the value of the "spread" parameter in the bitstream.
///
///     +--------------+------------------------+
///     | Spread value | f_r                    |
///     +--------------+------------------------+
///     | 0            | infinite (no rotation) |
///     | 1            | 15                     |
///     | 2            | 10                     |
///     | 3            | 5                      |
///     +--------------+------------------------+
///            Table 59: Spreading Values
///
/// The rotation angle is then calculated as
///
///     θ = (π * g_r ^ 2) / 4
///
/// A 2-D rotation R(i,j) between points x_i and x_j is defined as:
///
///     x_i =  cos(θ) * x_i + sin(θ) * x_j
///     x_j = -sin(θ) * x_i + cos(θ) * x_j
///
/// An N-D rotation is then achieved by applying a series of 2-D rotations back and forth, in the following order:
///
///     R(x_1, x_2), R(x_2, x_3), ..., R(x_N-2, X_N-1), R(x_N-1, X_N), R(x_N-2, X_N-1), ..., R(x_1, x_2).
///
/// If the decoded vector represents more than one time block, then this spreading process is applied separately on each time block. Also, if each block represents 8 samples or more, then another N-D rotation, by (π / 2 - θ), is applied _before_ the rotation described above. This extra rotation is applied in an interleaved manner with a stride equal to round(sqrt(N / nb_blocks)), i.e., it is applied independently for each set of sample
///
///     S_k = {stride * n + k}, n = 0 .. N / stride - 1
fn spread_vector(x: &mut [f32], num_blocks: usize, k: usize, spread: usize) {
    let n = x.len();
    const SPREAD_FACTOR: [usize; 4] = [std::usize::MAX, 15, 10, 5];
    let gain = n as f32 / (n + SPREAD_FACTOR[spread] * k) as f32;
    let theta = std::f32::consts::PI * gain.powi(2) / 4.0;

    let n_per_block = n / num_blocks;
    if n_per_block >= 8 {
        let stride = (n_per_block as f32).sqrt().round() as usize;
        for i in 0..num_blocks {
            rotate_block(&mut x[n_per_block * i..n_per_block * (i + 1)], stride, std::f32::consts::PI / 2.0 - theta);
        }
    }
    for i in 0..num_blocks {
        rotate_block(&mut x[n_per_block * i..n_per_block * (i + 1)], 1, theta);
    }
}

fn rotate_block(x: &mut [f32], stride: usize, theta: f32) {
    let c = theta.cos();
    let s = theta.sin();
    if stride > x.len() {
        return;
    }
    for i in 0..x.len() - stride {
        let temp = x[i];
        x[i] = x[i] * c - x[i + stride] * s;
        x[i + stride] = temp * s + x[i + stride] * c;
    }
    if stride * 2 > x.len() {
        return;
    }
    for i in (0..x.len() - stride * 2).rev() {
        let temp = x[i];
        x[i] = x[i] * c - x[i + stride] * s;
        x[i + stride] = temp * s + x[i + stride] * c;
    }
}


fn extract_collapse_mask(x: &mut [f32], stride: usize) -> u32 {
    if stride <= 1 {
        return 1;
    }
    let n = x.len() / stride;
    let mut collapse_mask = 0;
    for i in 0..stride {
        let mut tmp = 0;
        for j in 0..n {
            tmp |= x[i * n + j] as i32;
        }
        if tmp != 0 {
            collapse_mask |= 1 << i;
        }
    }
    return collapse_mask;
}

pub fn alg_unquant(x: &mut [f32],
                   k: usize,
                   spread: usize,
                   stride: usize,
                   ec: &mut entdec::EntropyCoder,
                   gain: f32,
                   v: &Vec<Vec<Option<u32>>>)
                   -> u32 {
    cwrs::decode_pulses(x, k, ec, v);
    let mask = extract_collapse_mask(x, stride);
    utils::renormalise(x, gain);
    if 2 * k < x.len() { // TODO: why?
        spread_vector(x, stride, k, spread);
    }
    return mask;
}
