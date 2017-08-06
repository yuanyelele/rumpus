use consts;

pub fn denormalise_bands(x: &mut [f32], bands: &[f32]) {
    /* Mean energy in each band quantized in Q4 */
    let e_means = vec![103, 100, 92, 85, 81, 77, 72, 70, 78, 75, 73, 71, 78, 74, 69, 72, 70, 74,
                       76, 71, 60];
    for i in 0..21 {
        let g = (bands[i] + e_means[i] as f32 / 16.0).exp2();
        for j in 8 * consts::BANDS[i]..8 * consts::BANDS[i + 1] {
            x[j] *= g;
        }
    }
    for i in 8 * consts::BANDS[21]..960 {
        x[i] = 0.0;
    }
}
