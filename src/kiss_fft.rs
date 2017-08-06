extern crate num_complex;

use consts;
use std;

pub struct KissFft {
    pub shift: i32,
    pub factors: Vec<usize>,
    pub bitrev: Vec<usize>,
    pub trig: Vec<f32>,
}

impl KissFft {
    pub fn new(shift: i32, factors: Vec<usize>) -> KissFft {
        let mut fft = KissFft {
            shift: shift,
            factors: factors,
            bitrev: vec![0; consts::FRAME_SIZE / 2 >> shift],
            trig: vec![0.0; consts::FRAME_SIZE >> shift],
        };

        for i in 0..fft.bitrev.len() {
            let mut a = 1;
            for k in 0..fft.factors.len() {
                fft.bitrev[i] = i % (a * fft.factors[k]) / a +
                                fft.factors[k] * fft.bitrev[i];
                a *= fft.factors[k];
            }
        }

        for i in 0..fft.trig.len() {
            let theta = std::f32::consts::PI * (i as f32 + 1.0 / 8.0) / fft.trig.len() as f32;
            fft.trig[i] = theta.cos();
        }

        return fft;
    }
}

fn butterfly(x0: &mut [num_complex::Complex<f32>],
             factor: usize, s: usize, m: usize, n: usize, mm: usize,
             cw: &[num_complex::Complex<f32>]) {
    let mut a = vec![num_complex::Complex::<f32>::new(0.0, 0.0); factor];
    let mut t = vec![num_complex::Complex::<f32>::new(0.0, 0.0); factor];
    for i in 0..n {
        for j in 0..m {
            let x = &mut x0[i * mm + j..];

            // e.g. factor == 5:
            // ┌     ┐     ╭                            ╮┌                 ┐
            // │x[0] │     │0,    0,     0,     0,     0││x[0]             │
            // │x[m] │     │0, 2π/5,  4π/5,  6π/5,  8π/5││x[m]cw(πjs/240)  │
            // │x[2m]│ = cw│0, 4π/5,  8π/5, 12π/5, 16π/5││x[2m]cw(2πjs/240)│
            // │x[3m]│     │0, 6π/5, 12π/5, 18π/5, 24π/5││x[3m]cw(3πjs/240)│
            // │x[4m]│     │0, 8π/5, 16π/5, 24π/5, 32π/5││x[4m]cw(4πjs/240)│
            // └     ┘     ╰                            ╯└                 ┘
            for k in 0..factor {
                a[k] = x[k * m] * cw[k * j * s % 480];
            }
            for k in 0..factor {
                t[k] = num_complex::Complex::<f32>::new(0.0, 0.0);
                for l in 0..factor {
                    t[k] = t[k] + a[l] * cw[480 / factor * k * l % 480];
                }
            }
            for k in 0..factor {
                x[k * m] = t[k];
            }
        }
    }
}

pub fn opus_fft(st: &KissFft,
                x: &mut [num_complex::Complex<f32>],
                cw: &[num_complex::Complex<f32>]) {
    let mut m2 = 1;
    for i in (0..st.factors.len()).rev() {
        let m = m2;
        m2 *= st.factors[i];
        let stride = 480 / m2;
        butterfly(x, st.factors[i], stride, m, stride >> st.shift, m2, cw);
    }
}
