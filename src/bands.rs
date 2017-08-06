use bands_utils;
use consts;
use entdec;
use opus_decoder;
use std;
use utils;
use vq;

pub const SPREAD_NORMAL: i32 = 2;

pub struct BandCtx {
    i: usize,
    intensity: usize,
    spread: i32,
    tf_change: i32,
    remaining_bits: i32,
    seed: u32,
}

pub struct SplitCtx {
    is_inv: bool,
    imid: i32,
    iside: i32,
    delta: i32,
    itheta: i32,
    qalloc: i32,
}

/// 4.3.4.1. Bits to Pulses
///
/// Although the allocation is performed in 1/8th bit units, the quantization requires an integer number of pulses k. To do this, the encoder searches for the value of k that produces the number of bits nearest to the allocated value (rounding down if exactly halfway between two values), not to exceed the total number of bits available. For efficiency reasons, the search is performed against a precomputed allocation table that only permits some k values for each n. The number of codebook entries can be computed as explained in Section 4.3.4.2. The difference between the number of bits allocated and the number of bits used is accumulated to a "balance" (initialized to zero) that helps adjust the allocation for the next bands. One third of the balance is applied to the bit allocation of each band to help achieve the target allocation. The only exceptions are the band before the last and the last band, for which half the balance and the whole balance are applied, respectively.
fn bits2pulses(cache: &[i32], bits: i32) -> usize {
    let mut lo = 0;
    let mut hi = cache.len() - 1;
    while hi > lo {
        let mid = (lo + hi) / 2;
        if cache[mid] >= bits {
            hi = mid;
        } else {
            lo = mid + 1;
        }
    }
    assert_eq!(lo, hi);
    if bits -
       if lo == 0 {
           0
       } else {
           cache[lo - 1]
       } <= cache[lo] - bits {
        return lo;
    } else {
        return lo + 1;
    }
}

fn quant_nosplit(cache: &[i32],
                 v: &Vec<Vec<Option<u32>>>,
                 ec: &mut entdec::EntropyCoder,
                 ctx: &mut BandCtx,
                 x: &mut [f32],
                 b: i32,
                 stride: usize,
                 lowband: Option<&mut [f32]>,
                 gain: f32,
                 mut fill: u32)
                 -> u32 {
    for q in (0..bits2pulses(cache, b)).rev() {
        let curr_bits = cache[q];
        if curr_bits <= ctx.remaining_bits {
            ctx.remaining_bits -= curr_bits;
            return vq::alg_unquant(x,
                                   bands_utils::get_pulses(q + 1),
                                   ctx.spread as usize,
                                   stride,
                                   ec,
                                   gain,
                                   v);
        }
    }

    fill &= (1 << stride as u32) - 1;
    if fill == 0 {
        for i in x {
            *i = 0.0;
        }
        return 0;
    }

    match lowband {
        None => for i in 0..x.len() {
            ctx.seed = utils::lcg_rand(ctx.seed);
            x[i] = ctx.seed as i32 as f32;
        },
        Some(v) => for i in 0..x.len() {
            ctx.seed = utils::lcg_rand(ctx.seed);
            let tmp = if ctx.seed & 0x8000 != 0 { 1.0 / 256.0 } else { -1.0 / 256.0 };
            x[i] = v[i] + tmp;
        },
    }
    utils::renormalise(x, gain);
    return fill;
}

fn compute_theta(qn: i32,
                 ec: &mut entdec::EntropyCoder,
                 ctx: &mut BandCtx,
                 sctx: &mut SplitCtx,
                 n: usize,
                 b: &mut i32,
                 b0: i32,
                 is_stereo: bool) {
    let tell = ec.tell_frac() as i32;
    sctx.is_inv = false;
    if is_stereo && ctx.i >= ctx.intensity {
        if *b > 16 {
            sctx.is_inv = ec.decode_bit_logp(2) == 1;
        }
        sctx.itheta = 0;
    } else {
        sctx.itheta = bands_utils::get_theta(ec, qn, b0, is_stereo) * 16384 / qn;
    }
    sctx.qalloc = ec.tell_frac() as i32 - tell;
    *b -= sctx.qalloc;

    sctx.imid = bands_utils::bitexact_cos(sctx.itheta as i16) as i32;
    sctx.iside = bands_utils::bitexact_cos(16384 - sctx.itheta as i16) as i32;
    sctx.delta = ((n as i32 - 1) * bands_utils::bitexact_log2tan(sctx.iside, sctx.imid) + (1 << 7)) >> 8;
}

/// 4.3.4.4. Split Decoding
///
/// To avoid the need for multi-precision calculations when decoding PVQ codevectors, the maximum size allowed for codebooks is 32 bits. When larger codebooks are needed, the vector is instead split in two sub-vectors of size n / 2. A quantized gain parameter with precision derived from the current allocation is entropy coded to represent the relative gains of each side of the split, and the entire decoding process is recursively applied. Multiple levels of splitting may be applied up to a limit of lm + 1 splits. The same recursive mechanism is applied for the joint coding of stereo audio.
fn quant_partition(v: &Vec<Vec<Option<u32>>>,
                   ec: &mut entdec::EntropyCoder,
                   ctx: &mut BandCtx,
                   x: &mut [f32],
                   mut b: i32,
                   mut b0: i32,
                   lowband: Option<&mut [f32]>,
                   mut lm: i32,
                   gain: f32,
                   mut fill: u32)
                   -> u32 {
    {
        let cache = consts::BITS_CACHE[consts::BITS_CACHE_INDEX[(lm + 1) as usize][ctx.i]];
        if lm == -1 || x.len() <= 2 || b <= cache[cache.len() - 1] as i32 + 11 {
            return quant_nosplit(cache, v, ec, ctx, x, b, if b0 == 0 { 1 } else { b0 as usize }, lowband, gain, fill);
        }
    }

    let n = x.len() / 2;
    let (x, y) = x.split_at_mut(n);
    lm -= 1;
    fill = (fill & 1) | (fill << 1);

    let mut sctx = SplitCtx {
        is_inv: false,
        imid: 0,
        iside: 0,
        delta: 0,
        itheta: 0,
        qalloc: 0,
    };
    let qn = bands_utils::compute_qn(consts::BAND_WIDTHS[ctx.i], n as i32, b, lm);
    compute_theta(qn, ec, ctx, &mut sctx, n, &mut b, b0, false);

    b0 /= 2;
    if b0 > 0 {
        if sctx.itheta > 8192 {
            sctx.delta -= sctx.delta >> (4 - lm);
        } else {
            sctx.delta = std::cmp::min(0, sctx.delta + (n as i32 >> (2 - lm)));
        }
    }

    if sctx.itheta == 0 {
        fill &= (1 << b0) - 1;
    }
    if sctx.itheta == 16384 {
        fill &= ((1 << b0) - 1) << b0;
    }

    let mut mbits = std::cmp::max(0, std::cmp::min(b, (b - sctx.delta) / 2));
    let mut sbits = b - mbits;
    let mid = sctx.imid as f32 / 32768.0;
    let side = sctx.iside as f32 / 32768.0;
    ctx.remaining_bits -= sctx.qalloc;
    let mut rebalance = ctx.remaining_bits;
    let mut cm;
    match lowband {
        None => {
            if mbits >= sbits {
                cm = quant_partition(v, ec, ctx, x, mbits, b0, None, lm, gain * mid, fill);
                rebalance = mbits - (rebalance - ctx.remaining_bits);
                if rebalance > 3 * 8 && sctx.itheta != 0 {
                    sbits += rebalance - 3 * 8;
                }
                cm |= quant_partition(v, ec, ctx, y, sbits, b0, None, lm, gain * side, fill >> b0) << b0;
            } else {
                cm = quant_partition(v, ec, ctx, y, sbits, b0, None, lm, gain * side, fill >> b0) << b0;
                rebalance = sbits - (rebalance - ctx.remaining_bits);
                if rebalance > 3 * 8 && sctx.itheta != 16384 {
                    mbits += rebalance - 3 * 8;
                }
                cm |= quant_partition(v, ec, ctx, x, mbits, b0, None, lm, gain * mid, fill);
            }
        },
        Some(lb) => {
            // >32-bit split case
            if mbits >= sbits {
                cm = quant_partition(v, ec, ctx, x, mbits, b0, Some(&mut lb[..n]), lm, gain * mid, fill);
                rebalance = mbits - (rebalance - ctx.remaining_bits);
                if rebalance > 3 * 8 && sctx.itheta != 0 {
                    sbits += rebalance - 3 * 8;
                }
                cm |= quant_partition(v, ec, ctx, y, sbits, b0, Some(&mut lb[n..]), lm, gain * side, fill >> b0) << b0;
            } else {
                cm = quant_partition(v, ec, ctx, y, sbits, b0, Some(&mut lb[n..]), lm, gain * side, fill >> b0) << b0;
                rebalance = sbits - (rebalance - ctx.remaining_bits);
                if rebalance > 3 * 8 && sctx.itheta != 16384 {
                    mbits += rebalance - 3 * 8;
                }
                cm |= quant_partition(v, ec, ctx, x, mbits, b0, Some(&mut lb[..n]), lm, gain * mid, fill);
            }
        },
    }
    return cm;
}

fn quant_band_mono(v: &Vec<Vec<Option<u32>>>,
                   ec: &mut entdec::EntropyCoder,
                   ctx: &mut BandCtx,
                   mut x: &mut [f32],
                   b: i32,
                   transient: bool,
                   gain: f32,
                   mut lowband: Option<&mut [f32]>,
                   mut fill: u32)
                   -> u32 {
    let mut scratch = vec![0.0; x.len()];

    let mut is_use_scratch = false;
    if ctx.i != 20 && (ctx.tf_change < 0 || transient) {
        if let Some(ref mut v) = lowband {
            is_use_scratch = true;
            for i in 0..x.len() {
                scratch[i] = v[i];
            }
        }
    }

    let recombine = std::cmp::max(0, ctx.tf_change);
    for i in 0..recombine {
        if let Some(ref mut v) = lowband {
            bands_utils::haar1(v, 1 << i);
        }
    }
    for _ in 0..recombine {
        const BIT_INTERLEAVE_TABLE: [u32; 16] = [0, 1, 1, 1, 2, 3, 3, 3, 2, 3, 3, 3, 2, 3, 3, 3];
        fill = BIT_INTERLEAVE_TABLE[fill as usize & 0xf] |
               BIT_INTERLEAVE_TABLE[fill as usize >> 4] * 4;
    }

    let mut b0 = if transient { 8 } else { 1 };
    b0 >>= recombine;

    let mut time_divide = 0;
    let mut tf_change = ctx.tf_change;
    while ((x.len() / b0) & 1) == 0 && tf_change < 0 {
        if let Some(ref mut v) = lowband {
            bands_utils::haar1(v, b0);
        }
        fill |= fill << b0;
        b0 *= 2;
        time_divide += 1;
        tf_change += 1;
    }

    // frequency order -> time order
    if let Some(ref mut v) = lowband {
        bands_utils::deinterleave_hadamard(v, b0 << recombine, !transient);
    }
    let mut cm;
    match lowband {
        None => {
            cm = quant_partition(v, ec, ctx, x, b, b0 as i32, None, 3, gain, fill);
        },
        Some(lb) => {
            cm = quant_partition(v, ec, ctx, x, b, b0 as i32, Some(lb), 3, gain, fill);
            if is_use_scratch {
                for i in 0..x.len() {
                    lb[i] = scratch[i];
                }
            }
        }
    }
    // time order -> frequency order
    if b0 > 1 {
        bands_utils::interleave_hadamard(&mut x, b0 << recombine, !transient);
    }

    /* Undo time-freq changes that we did earlier */
    for _ in 0..time_divide {
        b0 /= 2;
        cm |= cm >> b0;
        bands_utils::haar1(x, b0);
    }

    for i in 0..recombine {
        const BIT_DEINTERLEAVE_TABLE: [u8; 16] = [0x00, 0x03, 0x0C, 0x0F, 0x30, 0x33, 0x3C, 0x3F,
                                                  0xC0, 0xC3, 0xCC, 0xCF, 0xF0, 0xF3, 0xFC, 0xFF];
        cm = BIT_DEINTERLEAVE_TABLE[cm as usize] as u32;
        bands_utils::haar1(x, 1 << i);
    }

    if !transient {
        cm &= 1;
    }
    return cm;
}

fn quant_band_stereo(v: &Vec<Vec<Option<u32>>>,
                     ec: &mut entdec::EntropyCoder,
                     ctx: &mut BandCtx,
                     x: &mut [f32],
                     y: &mut [f32],
                     mut b: i32,
                     transient: bool,
                     lowband: Option<&mut [f32]>,
                     lowband_out: Option<&mut [f32]>,
                     mut fill: u32)
                     -> u32 {
    let b0 = if transient { 8 } else { 1 };

    let mut sctx = SplitCtx {
        is_inv: false,
        imid: 0,
        iside: 0,
        delta: 0,
        itheta: 0,
        qalloc: 0,
    };
    let qn = bands_utils::compute_qn(consts::BAND_WIDTHS[ctx.i], x.len() as i32, b, 3);
    compute_theta(qn, ec, ctx, &mut sctx, x.len(), &mut b, b0, true);

    if sctx.itheta == 0 {
        fill &= (1 << b0) - 1;
    }
    if sctx.itheta == 16384 {
        fill &= ((1 << b0) - 1) << b0;
    }

    let mut mbits = std::cmp::max(0, std::cmp::min(b, (b - sctx.delta) / 2));
    let mut sbits = b - mbits;
    let side = sctx.iside as f32 / 32768.0;
    ctx.remaining_bits -= sctx.qalloc;
    let mut rebalance = ctx.remaining_bits;
    let mut cm;
    if mbits >= sbits {
        cm = quant_band_mono(v, ec, ctx, x, mbits, transient, 1.0, lowband, fill);
        rebalance = mbits - (rebalance - ctx.remaining_bits);
        if rebalance > 3 * 8 && sctx.itheta != 0 {
            sbits += rebalance - 3 * 8;
        }
        cm |= quant_band_mono(v, ec, ctx, y, sbits, transient, side, None, fill >> b0);
    } else {
        cm = quant_band_mono(v, ec, ctx, y, sbits, transient, side, None, fill >> b0);
        rebalance = sbits - (rebalance - ctx.remaining_bits);
        if rebalance > 3 * 8 {
            mbits += rebalance - 3 * 8;
        }
        cm |= quant_band_mono(v, ec, ctx, x, mbits, transient, 1.0, lowband, fill);
    }
    // Scale output for later folding
    if let Some(v) = lowband_out {
        let f = (x.len() as f32).sqrt();
        for i in 0..x.len() {
            v[i] = f * x[i];
        }
    }

    bands_utils::stereo_merge(x, y, sctx.imid as f32 / 32768.0);
    if sctx.is_inv {
        for i in 0..y.len() {
            y[i] = -y[i];
        }
    }
    return cm;
}

pub fn quant_all_bands(st: &mut opus_decoder::OpusDecoder,
                       x: &mut [f32],
                       y: &mut [f32],
                       collapse_masks: &mut [u8],
                       pulses: &[i32],
                       transient: bool,
                       spread: i32,
                       mut is_dual_stereo: bool,
                       intensity: usize,
                       tf_res: &[i32],
                       total_bits: usize,
                       coded_bands: usize) {
    let mut ctx = BandCtx {
        i: 0,
        intensity: intensity,
        spread: spread,
        tf_change: 0,
        remaining_bits: 0,
        seed: st.range,
    };

    let mut lowband_offset = 0;
    let mut balance = st.ec.tell_frac() as i32;
    let mut is_update_lowband = true;
    let mut norm_x = vec![0.0; 8 * consts::BANDS[20]];
    let mut norm_y = vec![0.0; 8 * consts::BANDS[20]];
    for i in 0..21 {
        let band = 8 * consts::BANDS[i];
        let n = 8 * consts::BAND_WIDTHS[i];
        let tell = st.ec.tell_frac() as i32;
        ctx.i = i;
        ctx.remaining_bits = total_bits as i32 - tell - 1;
        ctx.tf_change = tf_res[i];
        balance -= tell;
        let b = if i < coded_bands {
            pulses[i] + balance / std::cmp::min(3, coded_bands - i) as i32
        } else {
            0
        };
        if is_update_lowband {
            lowband_offset = i;
        }
        let mut x_cm = 0;
        let mut y_cm = 0;
        let mut effective_lowband = 0;
        bands_utils::get_estimate(lowband_offset, collapse_masks, spread, ctx.tf_change, n,
                                  transient, &mut x_cm, &mut y_cm, &mut effective_lowband);

        if is_dual_stereo && i == intensity {
            is_dual_stereo = false;
            for j in 0..band {
                norm_x[j] = (norm_x[j] + norm_y[j]) / 2.0;
            }
        }
        if is_dual_stereo {
            {
                let lowband = if effective_lowband != -1 {
                    Some(&mut norm_x[effective_lowband as usize..effective_lowband as usize + n])
                } else {
                    None
                };
                x_cm = quant_band_mono(&st.mode.v, &mut st.ec, &mut ctx, &mut x[band..band + n], b / 2, transient, 1.0, lowband, x_cm as u32) as u8;
            }
            if i != 20 {
                // Scale output for later folding
                let f = (n as f32).sqrt();
                for j in 0..n {
                    norm_x[band + j] = f * x[band + j];
                }
            }
            {
                let lowband = if effective_lowband != -1 {
                    Some(&mut norm_y[effective_lowband as usize..effective_lowband as usize + n])
                } else {
                    None
                };
                y_cm = quant_band_mono(&st.mode.v, &mut st.ec, &mut ctx, &mut y[band..band + n], b / 2, transient, 1.0, lowband, y_cm as u32) as u8;
            }
            if i != 20 {
                // Scale output for later folding
                let f = (n as f32).sqrt();
                for j in 0..n {
                    norm_y[band + j] = f * y[band + j];
                }
            }
        } else {
            let (lb, lbo) = norm_x.split_at_mut(band);
            let lowband = if effective_lowband != -1 {
                Some(&mut lb[effective_lowband as usize..effective_lowband as usize + n])
            } else {
                None
            };
            let lowband_out = if i != 20 {
                Some(&mut lbo[..n])
            } else {
                None
            };
            x_cm = quant_band_stereo(&st.mode.v, &mut st.ec, &mut ctx, &mut x[band..band + n], &mut y[band..band + n], b, transient, lowband, lowband_out, (x_cm | y_cm) as u32) as u8;
            y_cm = x_cm;
        }
        collapse_masks[i * 2] = x_cm;
        collapse_masks[i * 2 + 1] = y_cm;
        balance += pulses[i] + tell;

        is_update_lowband = b > n as i32 * 8;
    }
    st.range = ctx.seed;
}
