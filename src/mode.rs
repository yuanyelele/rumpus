extern crate num_complex;

use kiss_fft;

pub struct CeltMode {
    pub window: Vec<f32>,
    pub fft0: kiss_fft::KissFft,
    pub fft3: kiss_fft::KissFft,
    pub twiddles: Vec<num_complex::Complex<f32>>,
    pub v: Vec<Vec<Option<u32>>>,
}
