use entdec;

const LAPLACE_NMIN: u32 = 16;

fn ec_laplace_decode(ec: &mut entdec::EntropyCoder, fs: u32, decay: u32) -> i32 {
    let fm = ec.decode(32768);
    if fm < fs {
        ec.update(0, fs, 32768);
        return 0;
    }

    let mut fl = fs;
    let mut val = 1;
    let mut fs = (32768 - 2 * LAPLACE_NMIN - fs) * (16384 - decay) / 32768 + 1;
    while fs > 1 && fm >= fl + fs * 2 {
        fl = fl + fs * 2;
        fs = (fs * 2 - 2) * decay / 32768 + 1;
        val += 1;
    }
    if fs <= 1 {
        let di = (fm - fl) / 2;
        val += di as i32;
        fl += 2 * di;
    }
    if fm < fl + fs {
        val = -val;
    } else {
        fl += fs;
    }
    ec.update(fl, fl + fs, 32768);
    return val;
}

const ALPHA_INTRA: f32 = 0.0;
const ALPHA_INTER: f32 = 0.5;

const BETA_INTRA: f32 = 27853.0 / 32768.0; // 0.85
const BETA_INTER: f32 = 26214.0 / 32768.0; // 0.8

const PROB_INTRA: [u32; 21] = [22, 63, 74, 84, 92, 103, 96, 96, 101, 107, 113, 118, 125, 118, 117,
                               135, 137, 157, 145, 97, 77];
const PROB_INTER: [u32; 21] = [42, 96, 108, 111, 117, 123, 120, 119, 127, 134, 139, 147, 152, 158,
                               154, 166, 173, 184, 184, 150, 139];

const DECAY_INTRA: [u32; 21] = [178, 114, 82, 83, 82, 62, 72, 67, 73, 72, 55, 52, 52, 52, 55, 49,
                                39, 32, 29, 33, 40];
const DECAY_INTER: [u32; 21] = [121, 66, 43, 40, 44, 32, 36, 33, 33, 34, 21, 23, 20, 25, 26, 21,
                                16, 13, 10, 13, 15];

pub fn unquant_coarse_energy(bands: &mut [f32], intra: bool, ec: &mut entdec::EntropyCoder) {
    let (alpha, beta) = if intra {
        (ALPHA_INTRA, BETA_INTRA)
    } else {
        (ALPHA_INTER, BETA_INTER)
    };

    let mut prev = [0.0; 2];
    for i in 0..bands.len() / 2 {
        for c in 0..2 {
            let q = ec_laplace_decode(ec,
                                      if intra {
                                          PROB_INTRA[i] << 7
                                      } else {
                                          PROB_INTER[i] << 7
                                      },
                                      if intra {
                                          DECAY_INTRA[i] << 6
                                      } else {
                                          DECAY_INTER[i] << 6
                                      });
            if bands[21 * c + i] < -9.0 {
                bands[21 * c + i] = -9.0;
            }
            bands[21 * c + i] = alpha * bands[21 * c + i] + prev[c] + q as f32;
            prev[c] = prev[c] + beta * q as f32;
        }
    }
}

pub fn unquant_fine_energy(bands: &mut [f32], fine_quant: &[i32], ec: &mut entdec::EntropyCoder) {
    for i in 0..bands.len() / 2 {
        for c in 0..2 {
            let q = ec.decode_bits(fine_quant[i] as usize);
            bands[21 * c + i] += (q as f32 + 0.5) / (1 << fine_quant[i]) as f32 - 0.5;
        }
    }
}

pub fn unquant_energy_finalise(bands: &mut [f32],
                               fine_quant: &[i32],
                               fine_priority: &[u32],
                               mut bits_left: i32,
                               ec: &mut entdec::EntropyCoder) {
    for prio in 0..2 {
        for i in 0..bands.len() / 2 {
            if bits_left < 2 {
                return;
            }
            if fine_priority[i] != prio {
                continue;
            }
            for c in 0..2 {
                let q = ec.decode_bits(1);
                bands[21 * c + i] += (q as f32 - 0.5) / (1 << (fine_quant[i] + 1)) as f32;
                bits_left -= 1;
            }
        }
    }
}
