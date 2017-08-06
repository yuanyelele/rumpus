use anti_collapse;
use bands;
use consts;
use denormalise_bands;
use entdec;
use mdct;
use mode;
use quant_bands;
use rate;
use std;

pub const BUFFER_SIZE: usize = 2048;

const GAINS: [[f32; 3]; 3] = [[0.3066406250, 0.2170410156, 0.1296386719],
                              [0.4638671875, 0.2680664062, 0.0],
                              [0.7998046875, 0.1000976562, 0.0]];

pub struct OpusDecoder<'a> {
    pub mode: mode::CeltMode,
    pub range: u32,
    pub pitch: usize,
    pub gain: f32,
    pub tapset: usize,
    pub preemph_mem: [f32; 2],
    pub decode_mem: [Vec<f32>; 2],
    pub bands: [f32; consts::NUM_BANDS * 6],
    pub ec: entdec::EntropyCoder<'a>,
}

pub fn comb_filter_old(window: &[f32], x: &mut [f32], pitch: usize, gain: f32, tapset: usize) {
    for i in 0..window.len() {
        let f = 1.0 - window[i] * window[i];
        x[i + pitch + 2] += f * gain *
                            (GAINS[tapset][0] * x[i + 2] +
                             GAINS[tapset][1] * (x[i + 3] + x[i + 1]) +
                             GAINS[tapset][2] * (x[i + 4] + x[i]));
    }
}

pub fn comb_filter(window: &[f32], x: &mut [f32], pitch: usize, gain: f32, tapset: usize) {
    for i in 0..window.len() {
        let f = window[i] * window[i];
        x[i + pitch + 2] += f * gain *
                            (GAINS[tapset][0] * x[i + 2] +
                             GAINS[tapset][1] * (x[i + 3] + x[i + 1]) +
                             GAINS[tapset][2] * (x[i + 4] + x[i]));
    }
}

pub fn comb_filter_const(x: &mut [f32], n: usize, pitch: usize, gain: f32, tapset: usize) {
    for i in 0..n {
        x[i + pitch + 2] += gain *
                            (GAINS[tapset][0] * x[i + 2] +
                             GAINS[tapset][1] * (x[i + 3] + x[i + 1]) +
                             GAINS[tapset][2] * (x[i + 4] + x[i]));
    }
}

fn deemphasis(x: &[Vec<f32>], pcm: &mut [f32], mem: &mut [f32]) {
    for c in 0..consts::NUM_CHANNELS as usize {
        for i in 0..consts::FRAME_SIZE {
            mem[c] += x[c][BUFFER_SIZE - consts::FRAME_SIZE + i];
            pcm[2 * i + c] = mem[c] / 32768.0;
            mem[c] *= consts::PRE_EMPHASIS;
        }
    }
}

unsafe fn celt_synthesis(st: &mut OpusDecoder, x: &mut [f32], is_transient: bool) {
    for c in 0..2 {
        denormalise_bands::denormalise_bands(&mut x[960 * c..],
                                             &st.bands[21 * c..21 * (c + 1)]);
        let shift = if is_transient { 3 } else { 0 };
        for b in 0..1 << shift {
            mdct::mdct_backward(&st.mode,
                                &x[960 * c + b..],
                                &mut st.decode_mem[c][BUFFER_SIZE - 960 + 120 * b..
                                                      BUFFER_SIZE - 960 + 120 * b + (960 >> shift) + 120 / 2],
                                shift);
        }
    }
}

pub fn tf_decode(is_transient: bool, tf_res: &mut [i32], ec: &mut entdec::EntropyCoder) {
    const TF_SELECT_TABLE: [[i8; 4]; 2] = [
                /* 20 ms */
                [ 0, -2, 0, -3 ], // !is_transient
                [ 3, 0, 1, -1 ] // is_transient
        ];

    let mut tf_changed = 0;
    let mut logp = if is_transient { 2 } else { 4 };
    tf_res[0] = 0;
    if ec.tell() + logp <= ec.buffer.len() * 8 {
        tf_res[0] = ec.decode_bit_logp(logp as u32) as i32;
        tf_changed = tf_res[0];
    }
    logp = if is_transient { 4 } else { 5 };
    for i in 1..tf_res.len() {
        tf_res[i] = tf_res[i - 1];
        if ec.tell() + logp <= ec.buffer.len() * 8 {
            tf_res[i] ^= ec.decode_bit_logp(logp as u32) as i32;
            tf_changed |= tf_res[i];
        }
    }
    let mut tf_select = 0;
    if TF_SELECT_TABLE[if is_transient { 1 } else { 0 }][tf_changed as usize] !=
       TF_SELECT_TABLE[if is_transient { 1 } else { 0 }][2 + tf_changed as usize] {
        tf_select = ec.decode_bit_logp(1) as usize;
    }
    for i in 0..tf_res.len() {
        tf_res[i] = TF_SELECT_TABLE[if is_transient { 1 } else { 0 }][2 * tf_select + tf_res[i] as usize] as i32;
    }
}

pub fn decode_post_filter_params(total_bits: usize,
                                 pitch: &mut usize,
                                 tapset: &mut usize,
                                 gain: &mut f32,
                                 ec: &mut entdec::EntropyCoder) {
    if ec.tell() + 16 <= total_bits && ec.decode_bit_logp(1) == 1 {
        let octave = ec.decode_uint(6) as usize;
        let pitch_in_octave = ec.decode_bits(4 + octave);
        *pitch = ((16 << octave) + pitch_in_octave) as usize - 1;
        *gain = 3.0 * (ec.decode_bits(3) + 1) as f32 / 32.0;
        const TAPSET_ICDF: [u8; 3] = [2, 1, 0];
        *tapset = ec.decode_icdf(&TAPSET_ICDF, 2) as usize;
    }
}

/// Decode the band boosts
///
/// First, set 'dynalloc_logp' to 6, the initial amount of storage required to signal a boost in bits, 'total_bits' to the size of the frame in 8th bits, 'total_boost' to zero, and 'tell' to the total number of 8th bits decoded so far. For each band from the coding start (0 normally, but 17 in Hybrid mode) to the coding end (which changes depending on the signaled bandwidth), the boost quanta in units of 1/8 bit is calculated as quanta = min(8 * N, max(48, N)). This represents a boost step size of six bits, subject to a lower limit of 1/8th bit/sample and an upper limit of 1 bit/sample. Set 'boost' to zero and 'dynalloc_loop_logp' to dynalloc_logp. While dynalloc_loop_log (the current worst case symbol cost) in 8th bits plus tell is less than total_bits plus total_boost and boost is less than cap[] for this band: Decode a bit from the bitstream with dynalloc_loop_logp as the cost of a one and update tell to reflect the current used capacity. If the decoded value is zero break the loop. Otherwise, add quanta to boost and total_boost, subtract quanta from total_bits, and set dynalloc_loop_log to 1. When the loop finishes 'boost' contains the bit allocation boost for this band. If boost is non-zero and dynalloc_logp is greater than 2, decrease dynalloc_logp. Once this process has been executed on all bands, the band boosts have been decoded.
fn decode_band_boosts(length: usize,
                      ec: &mut entdec::EntropyCoder,
                      boosts: &mut [i32]) -> usize {
    let mut dynalloc_logp = 6;
    let mut total_boost = 0;
    let mut tell = ec.tell_frac();
    for i in 0..boosts.len() {
        let width = 2 * consts::BAND_WIDTHS[i] * 8;
        let quanta = std::cmp::min(8 * width, std::cmp::max(48, width));
        let mut boost = 0;
        let mut dynalloc_loop_logp = dynalloc_logp;
        while dynalloc_loop_logp * 8 + tell < length * 8 * 8 {
            if ec.decode_bit_logp(dynalloc_loop_logp as u32) == 0 {
                break;
            }
            tell = ec.tell_frac();
            boost += quanta;
            total_boost += quanta;
            dynalloc_loop_logp = 1;
        }
        boosts[i] = boost as i32;
        if boost != 0 && dynalloc_logp > 2 {
            dynalloc_logp -= 1;
        }
    }
    return total_boost;
}

#[no_mangle]
pub unsafe extern "C" fn opus_decoder_decode(st: &mut OpusDecoder,
                                             data_ptr: *const u8,
                                             mut length: usize,
                                             pcm_ptr: *mut f32)
                                             -> usize {
    length -= 1;
    st.ec.init(std::slice::from_raw_parts(data_ptr.offset(1), length));

    let is_silence = st.ec.decode_bit_logp(15) == 1;

    let mut pitch: usize = 0;
    let mut tapset: usize = 0;
    let mut gain = 0.0;
    decode_post_filter_params(length * 8, &mut pitch, &mut tapset, &mut gain, &mut st.ec);

    let is_transient = st.ec.tell() + 3 <= length * 8 && st.ec.decode_bit_logp(3) == 1;

    let intra = st.ec.tell() + 3 <= length * 8 && st.ec.decode_bit_logp(3) == 1;

    quant_bands::unquant_coarse_energy(&mut st.bands[..2 * 21], intra, &mut st.ec);

    let mut tf_res = vec![0 as i32; 21];
    tf_decode(is_transient, &mut tf_res, &mut st.ec);

    const SPREAD_ICDF: [u8; 4] = [25, 23, 2, 0];
    let mut spread: i32 = bands::SPREAD_NORMAL;
    if st.ec.tell() + 4 <= length * 8 {
        spread = st.ec.decode_icdf(&SPREAD_ICDF, 5) as i32;
    }

    let mut boosts = vec![0; 21];
    let total_boost = decode_band_boosts(length, &mut st.ec, &mut boosts);

    // The allocation trim is an integer value from 0-10. The default value of 5 indicates no trim. The trim parameter is entropy coded in order to lower the coding cost of less extreme adjustments. Values lower than 5 bias the allocation towards lower frequencies and values above 5 bias it towards higher frequencies. Like other signaled parameters, signaling of the trim is gated so that it is not included if there is insufficient space available in the bitstream. To decode the trim, first set the trim value to 5, then if and only if the count of decoded 8th bits so far (ec.tell_frac) plus 48 (6 bits) is less than or equal to the total frame size in 8th bits minus total_boost (a product of the above band boost procedure), decode the trim value using the PDF in Table 58.
    const TRIM_ICDF: [u8; 11] = [126, 124, 119, 109, 87, 41, 19, 9, 4, 2, 0];
    let mut allocation_trim = 5;
    if st.ec.tell_frac() + 6 * 8 <= length * 8 * 8 - total_boost {
        allocation_trim = st.ec.decode_icdf(&TRIM_ICDF, 7) as i32;
    }

    let mut intensity = 0;
    let mut is_dual_stereo = false;
    let mut pulses = vec![0; 21];
    let mut fine_quant = vec![0; 21];
    let mut fine_priority = vec![0; 21];
    let coded_bands = rate::compute_allocation(&mut boosts,
                                               allocation_trim,
                                               &mut intensity,
                                               &mut is_dual_stereo,
                                               length,
                                               is_transient,
                                               &mut pulses,
                                               &mut fine_quant,
                                               &mut fine_priority,
                                               &mut st.ec);

    quant_bands::unquant_fine_energy(&mut st.bands[..2 * 21], &fine_quant, &mut st.ec);

    for c in 0..2 {
        for i in 0..BUFFER_SIZE - 960 + 120 / 2 {
            st.decode_mem[c][i] = st.decode_mem[c][i + 960];
        }
    }

    let mut collapse_masks = [0; 21 * 2];
    let mut x_y = [0.0; 960 * 2];
    {
        let (x, y) = x_y.split_at_mut(960);
        bands::quant_all_bands(st,
                               x,
                               y,
                               &mut collapse_masks,
                               &pulses,
                               is_transient,
                               spread,
                               is_dual_stereo,
                               intensity,
                               &tf_res,
                               length * 8 * 8 - if is_transient { 8 } else { 0 },
                               coded_bands);
    }

    let is_anti_collapse = is_transient && st.ec.decode_bits(1) == 1;
    quant_bands::unquant_energy_finalise(&mut st.bands[..2 * 21],
                                         &fine_quant,
                                         &fine_priority,
                                         length as i32 * 8 - st.ec.tell() as i32,
                                         &mut st.ec);

    if is_anti_collapse {
        anti_collapse::anti_collapse(&mut x_y,
                                     &collapse_masks,
                                     &st.bands[..2 * 21],
                                     &st.bands[2 * 21..4 * 21],
                                     &st.bands[4 * 21..6 * 21],
                                     &pulses,
                                     st.range);
    }

    if is_silence {
        for i in 0..2 * 21 {
            st.bands[i] = -28.0;
        }
    }

    celt_synthesis(st, &mut x_y, is_transient);

    for c in 0..2 {
        comb_filter_const(&mut st.decode_mem[c][BUFFER_SIZE - 960 - st.pitch - 2..], 120,
                          st.pitch, st.gain, st.tapset);
        comb_filter_old(&st.mode.window,
                        &mut st.decode_mem[c][BUFFER_SIZE - 960 + 120 - st.pitch - 2..],
                        st.pitch, st.gain, st.tapset);
        comb_filter(&st.mode.window,
                    &mut st.decode_mem[c][BUFFER_SIZE - 960 + 120 - pitch - 2..],
                    pitch, gain, tapset);
        comb_filter_const(&mut st.decode_mem[c][BUFFER_SIZE - 960 + 120 * 2 - pitch - 2..], 960 - 120 * 2,
                          pitch, gain, tapset);
    }
    st.pitch = pitch;
    st.gain = gain;
    st.tapset = tapset;

    if !is_transient {
        for i in 0..2 * 21 {
            st.bands[4 * 21 + i] = st.bands[2 * 21 + i];
            st.bands[2 * 21 + i] = st.bands[i];
        }
    }
    st.range = st.ec.range;

    let pcm = std::slice::from_raw_parts_mut(pcm_ptr, 960 * 2);
    deemphasis(&st.decode_mem, pcm, &mut st.preemph_mem);

    return 960;
}
