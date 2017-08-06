use consts;
use entdec;
use std;

fn lerp(bits1: &[i32], bits2: &[i32], thresh: &[i32], t: i32) -> i32 {
    let mut sum = 0;
    for i in 0..consts::NUM_BANDS {
        let bits = bits1[i] + (bits2[i] - bits1[i]) * t / 64;
        if bits >= thresh[i] {
            sum += bits;
        } else {
            sum += 8 * consts::NUM_CHANNELS;
        }
    }
    return sum;
}

fn get_lo(bits1: &[i32], bits2: &[i32], thresh: &[i32], total: i32) -> i32 {
    let mut lo = 0;
    let mut hi = 64;
    for _ in 0..6 {
        let mid = (lo + hi) / 2;
        let sum = lerp(bits1, bits2, thresh, mid);
        if sum > total {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    return lo;
}

fn get_bits(ebits: &mut [i32], t: i32, bits1: &[i32], bits2: &[i32], thresh: &[i32]) {
    for i in 0..consts::NUM_BANDS {
        let mut bits = bits1[i] + (bits2[i] - bits1[i]) * t / 64;
        if bits < thresh[i] {
            if bits >= 16 {
                bits = 16;
            } else {
                bits = 0;
            }
        }
        ebits[i] = bits;
    }
}

fn interp_bits2pulses(skip_start: usize,
                      bits1: &[i32],
                      bits2: &[i32],
                      thresh: &[i32],
                      mut total: i32,
                      intensity: &mut usize,
                      mut intensity_rsv: i32,
                      is_dual_stereo: &mut bool,
                      mut dual_stereo_rsv: i32,
                      bits: &mut [i32],
                      ebits: &mut [i32],
                      fine_priority: &mut [u32],
                      ec: &mut entdec::EntropyCoder)
                      -> usize {
    let t = get_lo(bits1, bits2, thresh, total);
    get_bits(bits, t, bits1, bits2, thresh);
    let mut sum = 0;
    for i in bits.iter() {
        sum += *i;
    }

    let mut left;
    let mut percoeff;
    let mut coded_bands = 0;
    for i in (0..consts::NUM_BANDS).rev() {
        if i <= skip_start {
            total += 8;
            coded_bands = i + 1;
            break;
        }
        percoeff = (total - sum) / consts::BANDS[i + 1] as i32;
        left = (total - sum) % consts::BANDS[i + 1] as i32;
        let mut band_bits = bits[i] + percoeff * consts::BAND_WIDTHS[i] as i32 +
                            std::cmp::max(left - consts::BANDS[i] as i32, 0);
        if band_bits >= thresh[i] {
            if ec.decode_bit_logp(1) == 1 {
                coded_bands = i + 1;
                break;
            }
            sum += 8;
            band_bits -= 8;
        }
        sum -= bits[i] + intensity_rsv;
        if intensity_rsv > 0 {
            intensity_rsv = (((i + 1) as f32).log2() * 8.0).ceil() as i32;
        }
        sum += intensity_rsv;
        if band_bits >= 16 {
            sum += 16;
            bits[i] = 16;
        } else {
            bits[i] = 0;
        }
    }

    *intensity = if intensity_rsv > 0 { ec.decode_uint(coded_bands as u32 + 1) as usize } else { 0 };
    if *intensity == 0 {
        total += dual_stereo_rsv;
        dual_stereo_rsv = 0;
    }
    *is_dual_stereo = if dual_stereo_rsv > 0 { ec.decode_bit_logp(1) == 1 } else { false };

    percoeff = (total - sum) / consts::BANDS[coded_bands] as i32;
    left = (total - sum) % consts::BANDS[coded_bands] as i32;
    for i in 0..coded_bands {
        bits[i] += percoeff * consts::BAND_WIDTHS[i] as i32;
    }
    for i in 0..coded_bands {
        let tmp = std::cmp::min(left, consts::BAND_WIDTHS[i] as i32);
        bits[i] += tmp;
        left -= tmp;
    }

    for i in 0..coded_bands {
        let n = consts::BAND_WIDTHS[i] as i32 * 8;
        let den = 2 * n + if *is_dual_stereo { 0 } else { 1 };
        let nc_logn = den * ((n as f32).log2() * 8.0).ceil() as i32;
        let mut offset = nc_logn / 2 - den * 21;
        if bits[i] + offset < den * 2 * 8 {
            offset += nc_logn / 4;
        } else if bits[i] + offset < den * 3 * 8 {
            offset += nc_logn / 8;
        }

        ebits[i] = std::cmp::max(0, bits[i] + offset + den * 4);
        ebits[i] = ebits[i] / den / 8;

        fine_priority[i] = if ebits[i] * den * 8 >= bits[i] + offset {
            1
        } else {
            0
        };
        bits[i] -= 2 * ebits[i] * 8;
    }

    for i in coded_bands..bits.len() {
        ebits[i] = bits[i] / 16;
        bits[i] = 0;
        fine_priority[i] = 0;
    }

    return coded_bands;
}

pub fn compute_allocation(boosts: &[i32],
                          allocation_trim: i32,
                          intensity: &mut usize,
                          is_dual_stereo: &mut bool,
                          length: usize,
                          is_transient: bool,
                          pulses: &mut [i32],
                          ebits: &mut [i32],
                          fine_priority: &mut [u32],
                          ec: &mut entdec::EntropyCoder)
                          -> usize {
    // = log2(consts::FRAME_SIZE / 120)
    let lm = 3;

    // The allocation computation begins by setting up some initial conditions. 'total' is set to the remaining available 8th bits, computed by taking the size of the coded frame times 8 and subtracting ec.tell_frac(). From this value, one (8th bit) is subtracted to ensure that the resulting allocation will be conservative. 'anti_collapse_rsv' is set to 8 (8th bits) if and only if the frame is a transient, LM is greater than 1, and total is greater than or equal to (LM+2) * 8.  Total is then decremented by anti_collapse_rsv and clamped to be equal to or greater than zero.  'skip_rsv' is set to 8 (8th bits) if total is greater than 8, otherwise it is zero.  Total is then decremented by skip_rsv.  This reserves space for the final skipping flag.
    let mut total = length as i32 * 8 * 8 - ec.tell_frac() as i32 - 1;
    let anti_collapse_rsv = if is_transient && lm > 1 && total >= (lm + 2) * 8 { 8 } else { 0 };
    total -= std::cmp::max(anti_collapse_rsv, 0);
    let skip_rsv = if total > 8 { 8 } else { 0 };
    total -= skip_rsv;

    // If the current frame is stereo, intensity_rsv is set to the conservative log2 in 8th bits of the number of coded bands for this frame (given by the table LOG2_FRAC_TABLE in rate.c).  If intensity_rsv is greater than total, then intensity_rsv is set to zero.  Otherwise, total is decremented by intensity_rsv, and if total is still greater than 8, dual_stereo_rsv is set to 8 and total is decremented by dual_stereo_rsv.
    let mut dual_stereo_rsv = 0;
    let mut intensity_rsv = (((consts::NUM_BANDS + 1) as f32).log2() * 8.0).ceil() as i32;
    if intensity_rsv > total {
        intensity_rsv = 0;
    } else {
        dual_stereo_rsv = 8;
        total -= intensity_rsv + 8;
    }

    let mut thresh = [0; consts::NUM_BANDS];
    let mut trim_offsets = [0; consts::NUM_BANDS];
    for i in 0..21 {
        // The allocation process then computes a vector representing the hard minimum amounts allocation any band will receive for shape. This minimum is higher than the technical limit of the PVQ process, but very low rate allocations produce an excessively sparse spectrum and these bands are better served by having no allocation at all. For each coded band, set thresh[band] to 24 times the number of MDCT bins in the band and divide by 16. If 8 times the number of channels is greater, use that instead. This sets the minimum allocation to one bit per channel or 48 128th bits per MDCT bin, whichever is greater. The band-size dependent part of this value is not scaled by the channel count, because at the very low rates where this limit is applicable there will usually be no bits allocated to the side.
        thresh[i] = std::cmp::max(24 * ((consts::BAND_WIDTHS[i] as i32) << lm) / 16, 8 * consts::NUM_CHANNELS);
        // The previously decoded allocation trim is used to derive a vector of per-band adjustments, 'trim_offsets[]'. For each coded band take the alloc_trim and subtract 5 and LM. Then, multiply the result by the number of channels, the number of MDCT bins in the shortest frame size for this mode, the number of remaining bands, 2**LM, and 8. Next, divide this value by 64. Finally, if the number of MDCT bins in the band per channel is only one, 8 times the number of channels is subtracted in order to diminish the allocation by one bit, because width 1 bands receive greater benefit from the coarse energy coding.
        trim_offsets[i] = consts::BAND_WIDTHS[i] as i32 * (allocation_trim - 5 - lm) * consts::NUM_CHANNELS * (consts::NUM_BANDS - 1 - i) as i32;
    }
    // The "static" bit allocation (in 1/8 bits) for a quality q, excluding the minimums, maximums, tilt and boosts, is equal to channels * N * alloc[band][q] << LM >> 2, where alloc[][] is given in Table 57 and LM = log2(frame_size / 120). The allocation is obtained by linearly interpolating between two values of q (in steps of 1/64) to find the highest allocation that does not exceed the number of bits remaining.
    let mut lo = 1;
    let mut hi = consts::NUM_QUALITIES - 1;
    while lo <= hi {
        let mut sum = 0;
        let mid = (lo + hi) / 2;
        for i in 0..consts::NUM_BANDS {
            let mut bits = consts::NUM_CHANNELS * consts::BAND_WIDTHS[i] as i32 *
                           consts::BAND_ALLOCATION[mid][i] << lm >> 2;
            bits = std::cmp::max(0, bits + trim_offsets[i]);
            bits += boosts[i];
            if bits >= thresh[i] {
                sum += bits;
            } else {
                sum += 8 * consts::NUM_CHANNELS;
            }
        }
        if sum > total {
            hi = mid - 1;
        } else {
            lo = mid + 1;
        }
    }

    let mut bits1 = [0; consts::NUM_BANDS];
    let mut bits2 = [0; consts::NUM_BANDS];
    let mut skip_start = 0;
    for i in 0..consts::NUM_BANDS {
        bits1[i] = consts::NUM_CHANNELS * consts::BAND_WIDTHS[i] as i32 *
                   consts::BAND_ALLOCATION[lo - 1][i] << lm >> 2;
        bits2[i] = consts::NUM_CHANNELS * consts::BAND_WIDTHS[i] as i32 *
                   consts::BAND_ALLOCATION[lo][i] << lm >> 2;
        bits1[i] = std::cmp::max(0, bits1[i] + trim_offsets[i]);
        bits2[i] = std::cmp::max(0, bits2[i] + trim_offsets[i]);
        bits1[i] += boosts[i];
        bits2[i] += boosts[i];

        if boosts[i] > 0 {
            skip_start = i;
        }
    }
    return interp_bits2pulses(skip_start,
                              &bits1,
                              &bits2,
                              &thresh,
                              total,
                              intensity,
                              intensity_rsv,
                              is_dual_stereo,
                              dual_stereo_rsv,
                              pulses,
                              ebits,
                              fine_priority,
                              ec);
}
