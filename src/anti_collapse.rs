use consts;
use std;
use utils;

pub fn anti_collapse(x: &mut [f32],
                     collapse_masks: &[u8],
                     log_e: &[f32],
                     prev1log_e: &[f32],
                     prev2log_e: &[f32],
                     pulses: &[i32],
                     mut seed: u32) {
    for i in 0..pulses.len() {
        let n0 = consts::BANDS[i + 1] - consts::BANDS[i];
        let depth = (1 + pulses[i]) / n0 as i32 / 8; // depth in 1/8 bits
        let thresh = (-depth as f32 / 8.0).exp2() / 2.0;
        for c in 0..2 {
            let e_diff = prev1log_e[21 * c + i].min(prev2log_e[21 * c + i]) - log_e[21 * c + i];
            let temp = x.len() / 2 * c;
            let x_off = &mut x[temp + 8 * consts::BANDS[i]..temp + 8 * consts::BANDS[i + 1]];
            let r = (e_diff.exp2() * 2.0 * std::f32::consts::SQRT_2).min(thresh) /
                    (x_off.len() as f32).sqrt();
            let mut renormalize = false;
            for k in 0..8 {
                if collapse_masks[i * 2 + c] & (1 << k) == 0 {
                    /* Fill with noise */
                    for j in 0..n0 as usize {
                        seed = utils::lcg_rand(seed);
                        x_off[k + j * 8] = if seed & 0x8000 != 0 { r } else { -r };
                    }
                    renormalize = true;
                }
            }
            /* We just added some energy, so we need to renormalise */
            if renormalize {
                utils::renormalise(x_off, 1.0);
            }

        }
    }
}
