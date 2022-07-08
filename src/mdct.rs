extern crate num_complex;

use kiss_fft;
use mode;
use std;

fn pre_rotate(fft: &kiss_fft::KissFft, x: &[f32], y: &mut [f32]) {
    let stride = 1 << fft.shift;
    let n = y.len();
    for i in 0..n / 2 {
        y[2 * fft.bitrev[i]] =
            x[stride * 2 * i] * fft.trig[i] - x[stride * (n - 1 - 2 * i)] * fft.trig[n / 2 + i];
        y[2 * fft.bitrev[i] + 1] =
            x[stride * (n - 1 - 2 * i)] * fft.trig[i] + x[stride * 2 * i] * fft.trig[n / 2 + i];
    }
}

fn post_rotate(y: &mut [f32]) {
    let n = y.len();
    for i in 0..(n / 2 + 1) / 2 {
        let i0 = 2 * i;
        let i1 = 2 * i + 1;
        let i2 = n - 2 * i - 1;
        let i3 = n - 2 * i - 2;

        let theta = -std::f32::consts::PI * (i as f32 + 1.0 / 8.0) / n as f32;
        let phi = std::f32::consts::PI * (i as f32 + 7.0 / 8.0) / n as f32;
        let n0 = y[i1] * theta.cos() - (-y[i0]) * theta.sin();
        let n2 = y[i1] * theta.sin() + (-y[i0]) * theta.cos();
        let n3 = (-y[i3]) * phi.cos() - (-y[i2]) * phi.sin();
        let n1 = (-y[i3]) * phi.sin() + (-y[i2]) * phi.cos();

        y[i0] = n0;
        y[i1] = n1;
        y[i2] = n2;
        y[i3] = n3;
    }
}

fn mirror(y: &mut [f32], window: &[f32]) {
    let n = y.len();
    for i in 0..n / 2 {
        let y1 = y[i];
        let y2 = y[n - 1 - i];
        y[i] = window[n - 1 - i] * y1 - window[i] * y2;
        y[n - 1 - i] = window[i] * y1 + window[n - 1 - i] * y2;
    }
}

pub unsafe fn mdct_backward(mode: &mode::CeltMode, x: &[f32], y: &mut [f32], shift: usize) {
    let complex = std::slice::from_raw_parts_mut(
        (&mut y[120 / 2..]).as_mut_ptr() as *mut num_complex::Complex<f32>,
        480 >> shift,
    );
    let fft = if shift == 0 { &mode.fft0 } else { &mode.fft3 };
    pre_rotate(fft, x, &mut y[120 / 2..120 / 2 + (960 >> shift)]);
    kiss_fft::opus_fft(fft, complex, &mode.twiddles);
    post_rotate(&mut y[120 / 2..120 / 2 + (960 >> shift)]);
    mirror(&mut y[..120], &mode.window);
}
