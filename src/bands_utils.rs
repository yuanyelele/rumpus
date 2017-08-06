use consts;
use entdec;
use std;
use utils;

pub fn frac_mul16(a: i32, b: i32) -> i32 {
    return (16384 + a * b) >> 15;
}

pub fn get_pulses(i: usize) -> usize {
    return if i < 8 {
               i
           } else {
               (8 + (i & 7)) << ((i >> 3) - 1)
           };
}

// TODO
pub fn bitexact_cos(x: i16) -> i16 {
    if x == 0 {
        return 32767;
    }
    if x == 16384 {
        return 0;
    }
    let mut x2 = (4096 + x as i32 * x as i32) >> 13;
    x2 = 32767 - x2 + frac_mul16(x2, -7651 + frac_mul16(x2, 8277 + frac_mul16(-626, x2)));
    return 1 + x2 as i16;
}

// TODO
// log2(isin / icos) * 2048
pub fn bitexact_log2tan(mut isin: i32, mut icos: i32) -> i32 {
    if isin == 0 {
        return -32768;
    }
    if icos == 0 {
        return 32768;
    }
    let lc = (icos as f32).log2() as i32 + 1;
    let ls = (isin as f32).log2() as i32 + 1;
    icos <<= 15 - lc;
    isin <<= 15 - ls;
    return (ls - lc) * (1 << 11) + frac_mul16(isin, frac_mul16(isin, -2597) + 7932) -
           frac_mul16(icos, frac_mul16(icos, -2597) + 7932);
}

pub fn stereo_merge(x: &mut [f32], y: &mut [f32], m: f32) {
    let xp = utils::inner_product(x, y) * m;
    let side = utils::inner_product(y, y);
    let el = ((m * m) as f32 + side - 2.0 * xp).sqrt();
    let er = ((m * m) as f32 + side + 2.0 * xp).sqrt();
    for i in 0..x.len() as usize {
        let l = m * x[i];
        let r = y[i];
        x[i] = 1.0 / el * (l - r);
        y[i] = 1.0 / er * (l + r);
    }
}

const ORDERY_TABLE: [usize; 31] = [0, 1, 0, 3, 0, 2, 1, 7, 0, 4, 3, 6, 1, 5, 2, 15, 0, 8, 7, 12,
                                   3, 11, 4, 14, 1, 9, 6, 13, 2, 10, 5];

pub fn interleave_hadamard(x: &mut [f32], stride: usize, hadamard: bool) {
    let mut y = vec![0.0; x.len()];
    let n0 = x.len() / stride;
    let ordery = &ORDERY_TABLE[stride - 1..stride * 2 - 1];
    for i in 0..stride {
        for j in 0..n0 {
            y[stride * j + i] = x[n0 * if hadamard { ordery[i] } else { i } + j];
        }
    }
    for i in 0..x.len() {
        x[i] = y[i];
    }
}

pub fn deinterleave_hadamard(y: &mut [f32], stride: usize, hadamard: bool) {
    let mut x = vec![0.0; y.len()];
    let n0 = y.len() / stride;
    let ordery = &ORDERY_TABLE[stride - 1..stride * 2 - 1];
    for i in 0..stride {
        for j in 0..n0 {
            x[n0 * if hadamard { ordery[i] } else { i } + j] = y[stride * j + i];
        }
    }
    for i in 0..y.len() {
        y[i] = x[i];
    }
}

pub fn haar1(x: &mut [f32], stride: usize) {
    let n0 = x.len() / stride;
    for i in 0..stride {
        for j in 0..n0 / 2 {
            let tmp1 = std::f32::consts::SQRT_2 / 2.0 * x[stride * (2 * j) + i];
            let tmp2 = std::f32::consts::SQRT_2 / 2.0 * x[stride * (2 * j + 1) + i];
            x[stride * (2 * j) + i] = tmp1 + tmp2;
            x[stride * (2 * j + 1) + i] = tmp1 - tmp2;
        }
    }
}

const QTHETA_OFFSET: i32 = 4;

pub fn compute_qn(band_width: usize, n: i32, b: i32, lm: i32) -> i32 {
    let pulse_cap = ((band_width as f32).log2() * 8.0).ceil() as i32 + lm * 8;
    let offset = pulse_cap / 2 - QTHETA_OFFSET;
    let qb = b / (2 * n - 1) + offset;
    let qn = ((qb as f32 / 8.0).exp2() / 2.0).round() as i32 * 2;
    return std::cmp::min(qn, 256);
}

pub fn get_theta(ec: &mut entdec::EntropyCoder, qn: i32, b0: i32, is_stereo: bool) -> i32 {
    let itheta;
    if is_stereo {
        let p0 = 3;
        let x0 = qn / 2;
        let ft = p0 * (x0 + 1) + x0;
        let fs = ec.decode(ft as u32) as i32;
        if fs < (x0 + 1) * p0 {
            itheta = fs / p0;
        } else {
            itheta = x0 + 1 + (fs - (x0 + 1) * p0);
        }
        ec.update(if itheta <= x0 {
                      (p0 * itheta) as u32
                  } else {
                      ((itheta - 1 - x0) + (x0 + 1) * p0) as u32
                  },
                  if itheta <= x0 {
                      (p0 * (itheta + 1)) as u32
                  } else {
                      ((itheta - x0) + (x0 + 1) * p0) as u32
                  },
                  ft as u32);
    } else if b0 > 1 {
        itheta = ec.decode_uint(qn as u32 + 1) as i32;
    } else {
        let ft = (qn / 2 + 1) * (qn / 2 + 1);
        let fm = ec.decode(ft as u32) as i32;
        let fl;
        let fs;
        if fm < qn * (qn / 2 + 1) / 4 {
            itheta = (((8 * fm + 1) as f32).sqrt() as i32 - 1) / 2;
            fs = itheta + 1;
            fl = fs * (fs - 1) / 2;
        } else {
            itheta = (2 * (qn + 1) - (((8 * (ft - fm - 1) + 1) as f32).sqrt() as i32)) / 2;
            fs = qn + 1 - itheta;
            fl = ft - fs * (fs + 1) / 2;
        }
        ec.update(fl as u32, (fl + fs) as u32, ft as u32);
    }
    return itheta;
}

const SPREAD_AGGRESSIVE: i32 = 3;

pub fn get_estimate(offset: usize,
                    collapse_masks: &[u8],
                    spread: i32,
                    tf_change: i32,
                    n: usize,
                    transient: bool,
                    x_cm: &mut u8,
                    y_cm: &mut u8,
                    lowband: &mut i32) {
    if offset == 0 || (spread == SPREAD_AGGRESSIVE && tf_change >= 0) {
        *lowband = -1;
        if transient {
            *x_cm = 255;
            *y_cm = 255;
        } else {
            *x_cm = 1;
            *y_cm = 1;
        }
    } else {
        *lowband = 8 * consts::BANDS[offset] as i32 - n as i32;
        *x_cm = 0;
        *y_cm = 0;
        let mut i = offset;
        while 8 * consts::BANDS[i] as i32 > *lowband {
            *x_cm |= collapse_masks[i * 2];
            *y_cm |= collapse_masks[i * 2 + 1];
            i -= 1;
        }
    }
}
