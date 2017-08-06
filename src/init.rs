extern crate num_complex;

use consts;
use kiss_fft;
use mode;
use opus_decoder;
use std;

#[no_mangle]
pub extern "C" fn opus_decoder_create<'a>() -> *mut opus_decoder::OpusDecoder<'a> {
    let mut mode = mode::CeltMode {
        window: vec![0.0; consts::WINDOW_SIZE],
        fft0: kiss_fft::KissFft::new(0, vec![5, 4, 4, 3, 2]),
        fft3: kiss_fft::KissFft::new(3, vec![5, 4, 3]),
        twiddles: vec![Default::default(); consts::FRAME_SIZE / 2],
        v: vec![vec![None; 176]; 176],
    };

    for i in 0..mode.window.len() {
        let theta = 0.5 * std::f32::consts::PI * (i as f32 + 0.5) / mode.window.len() as f32;
        mode.window[i] = (0.5 * std::f32::consts::PI * theta.sin().powi(2)).sin();
    }

    for i in 0..mode.twiddles.len() {
        let theta = 2.0 * std::f32::consts::PI * (i as f32) / mode.twiddles.len() as f32;
        mode.twiddles[i] = num_complex::Complex::new(theta.cos(), -theta.sin());
    }

    // The number of n-dimensional unit pulse vectors with k pulses.
    //
    // The number of combinations, with replacement, of n items, taken k at a time, when a sign bit is added to each item taken at least once. A table of values for n < 10 and k < 10 looks like:
    //
    //     v[10][10] = {
    //         {1,  0,   0,    0,    0,     0,     0,      0,      0,       0},
    //         {1,  2,   2,    2,    2,     2,     2,      2,      2,       2},
    //         {1,  4,   8,   12,   16,    20,    24,     28,     32,      36},
    //         {1,  6,  18,   38,   66,   102,   146,    198,    258,     326},
    //         {1,  8,  32,   88,  192,   360,   608,    952,   1408,    1992},
    //         {1, 10,  50,  170,  450,  1002,  1970,   3530,   5890,    9290},
    //         {1, 12,  72,  292,  912,  2364,  5336,  10836,  20256,   35436},
    //         {1, 14,  98,  462, 1666,  4942, 12642,  28814,  59906,  115598},
    //         {1, 16, 128,  688, 2816,  9424, 27008,  68464, 157184,  332688},
    //         {1, 18, 162,  978, 4482, 16722, 53154, 148626, 374274,  864146}
    //     };
    mode.v[0][0] = Some(1);
    for k in 1..mode.v[0].len() {
        mode.v[0][k] = Some(0);
    }
    for n in 1..mode.v.len() {
        mode.v[n][0] = Some(1);
        for k in 1..mode.v[n].len() {
	    match (mode.v[n - 1][k], mode.v[n][k - 1], mode.v[n - 1][k - 1]) {
	        (Some(a), Some(b), Some(c)) => {
		    let (temp, overflow) = a.overflowing_add(b);
		    if !overflow {
		        let (temp, overflow) = temp.overflowing_add(c);
			if !overflow {
			    mode.v[n][k] = Some(temp);
			}
		    }
		},
		_ => {},
	    }
        }
    }

    let opus_decoder = opus_decoder::OpusDecoder {
        mode: mode,
        range: 0,
        pitch: 0,
        gain: 0.0,
        tapset: 0,
        preemph_mem: [0.0; 2],
        decode_mem: [
            vec![0.0; opus_decoder::BUFFER_SIZE + 120 / 2],
            vec![0.0; opus_decoder::BUFFER_SIZE + 120 / 2]
        ],
        bands: [0.0; consts::NUM_BANDS * 6],
        ec: Default::default(),
    };
    let b = Box::new(opus_decoder);
    return Box::into_raw(b);
}
