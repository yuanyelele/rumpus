/// Linear congruential generator.
///
///     x = (a * x + c) mod m
///
/// In _Numerical Recipes_, m = 2^32 , a = 1664525, c = 1013904223
pub fn lcg_rand(seed: u32) -> u32 {
    return (1664525 * seed as u64 + 1013904223) as u32;
}

/// Inner product.
pub fn inner_product(x: &[f32], y: &[f32]) -> f32 {
    assert_eq!(x.len(), y.len());
    let mut xy = 0.0;
    for i in 0..x.len() {
        xy += x[i] * y[i];
    }
    return xy;
}

/// Normalise.
///
///     x *= gain / |x|
pub fn renormalise(x: &mut [f32], gain: f32) {
    let scale = gain / inner_product(x, x).sqrt();
    for i in x {
        *i *= scale;
    }
}
