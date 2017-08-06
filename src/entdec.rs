use std;

#[derive(Default)]
pub struct EntropyCoder<'a> {
    pub buffer: &'a [u8],
    pub head: usize,
    pub tail: usize,
    pub bits_buffer: u32,
    pub bits_buffer_size: usize,
    pub nbits_total: usize,
    pub range: u32,
    pub value: u32,
    pub norm_factor: u32,
    pub remainder: u8,
}

impl<'a> EntropyCoder<'a> {
    fn read_byte(&mut self) -> u8 {
        self.nbits_total += 8;
        if self.head >= self.buffer.len() {
            return 0;
        }
        let byte = self.buffer[self.head];
        self.head += 1;
        return byte;
    }

    fn read_byte_from_end(&mut self) -> u8 {
        self.tail += 1;
        return self.buffer[self.buffer.len() - self.tail];
    }

    fn normalize(&mut self) {
        let ec_code_min_range = 1 << 23;
        let ec_code_top = (1 << 31) - 1;
        while self.range <= ec_code_min_range {
            self.range <<= 8;
            let mut sym = self.remainder;
            self.remainder = self.read_byte();
            sym = (((sym as u16) << 8 | self.remainder as u16) >> 1) as u8;
            let rsym = !sym;
            self.value = ((self.value << 8) + rsym as u32) & ec_code_top;
        }
    }

    pub fn init(&mut self, buffer: &'a [u8]) {
        let ec_code_init_range = 1 << 7;
        self.buffer = buffer;
        self.head = 0;
        self.tail = 0;
        self.nbits_total = 1;
        self.remainder = self.read_byte();
        self.range = ec_code_init_range;
        self.value = self.range - (self.remainder >> 1) as u32 - 1;
        self.bits_buffer = 0;
        self.bits_buffer_size = 0;
        self.normalize();
    }

    pub fn tell(&self) -> usize {
        let lg = (self.range as f32).log2() as usize + 1;
        return self.nbits_total - lg;
    }

    pub fn tell_frac(&self) -> usize {
        let mut lg = (self.range as f32).log2() as usize + 1;
        let mut r = self.range;
        for _ in 0..3 {
            r >>= (r as f32).log2() as usize + 1 - 16;
            r = r * r;
            lg = lg * 2 + (r >> 31) as usize;
        }
        return self.nbits_total * 8 - lg;
    }

    pub fn decode(&mut self, ft: u32) -> u32 {
        self.norm_factor = self.range / ft;
        return ft - std::cmp::min(self.value / self.norm_factor + 1, ft);
    }

    pub fn update(&mut self, fl: u32, fh: u32, ft: u32) {
        let temp = self.norm_factor * (ft - fh);
        self.value -= temp;
        self.range = if fl > 0 {
            self.norm_factor * (fh - fl)
        } else {
            self.range - temp
        };
        self.normalize();
    }

    pub fn decode_icdf(&mut self, icdf: &[u8], ftb: u32) -> u8 {
        let ft = 1 << ftb;
        let fs = self.decode(ft);
        let mut k = 0;
        while fs >= (ft - icdf[k] as u32) {
            k += 1;
        }
        let mut fl = 0;
        if k > 0 {
            fl = ft - icdf[k - 1] as u32;
        }
        let fh = ft - icdf[k] as u32;
        self.update(fl, fh, ft);
        return k as u8;
    }

    pub fn decode_bit_logp(&mut self, logp: u32) -> u8 {
        let ft = 1 << logp;
        let fs = self.decode(ft);
        let bit = if fs < ft - 1 {
            self.update(0, ft - 1, ft);
            0
        } else {
            self.update(ft - 1, ft, ft);
            1
        };
        return bit;
    }

    pub fn decode_uint(&mut self, ft: u32) -> u32 {
        let ftb = ((ft - 1) as f32).log2() as usize + 1;
        let mut t: u32;
        if ftb <= 8 {
            t = self.decode(ft);
            self.update(t, t + 1, ft);
        } else {
            let rft = ((ft - 1) >> (ftb - 8)) + 1;
            t = self.decode(rft);
            self.update(t, t + 1, rft);
            t = t << (ftb - 8) | self.decode_bits(ftb - 8);
            if t >= ft {
                t = ft - 1;
            }
        }
        return t;
    }

    pub fn decode_bits(&mut self, bits: usize) -> u32 {
        while self.bits_buffer_size < bits {
            self.bits_buffer |= (self.read_byte_from_end() as u32) << self.bits_buffer_size;
            self.bits_buffer_size += 8;
        }
        let ret = self.bits_buffer & ((1 << bits) - 1);
        self.bits_buffer >>= bits;
        self.bits_buffer_size -= bits;
        self.nbits_total += bits;
        return ret;
    }
}
