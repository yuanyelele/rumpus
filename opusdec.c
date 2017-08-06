#include <limits.h>
#include <math.h>
#include <ogg/ogg.h>
#include <speex/speex_resampler.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "opus.h"

/* 120ms at 48000 */
#define MAX_FRAME_SIZE (960*6)

typedef struct {
	float *b_buf;
	float *a_buf;
} shapestate;

unsigned rngseed = 22222;
inline unsigned fast_rand(void)
{
	rngseed = (rngseed * 96314165) + 907633515;
	return rngseed;
}

inline void shape_dither_toshort(shapestate *_ss, short *_o, float *_i, int _n)
{
	const float fcoef[8] = {
		2.2061, -0.4706, -0.2534, -0.6214, 1.0587, 0.0676, -0.6054, -0.2738,	/* 44.1kHz noise shaping filter sd=2.51 */
	};
	float gain = 32768 - 15;
	float *b_buf = _ss->b_buf;
	float *a_buf = _ss->a_buf;
	/*In order to avoid replacing digital silence with quiet dither noise
	   we mute if the output has been silent for a while */
	for (int i = 0; i < _n; i++) {
		int pos = i * 2;
		for (int c = 0; c < 2; c++) {
			float err = 0;
			float s = _i[pos + c] * gain;
			for (int j = 0; j < 4; j++)
				err += fcoef[j] * b_buf[c * 4 + j] - fcoef[j + 4] * a_buf[c * 4 + j];
			memmove(&a_buf[c * 4 + 1], &a_buf[c * 4], sizeof(float) * 3);
			memmove(&b_buf[c * 4 + 1], &b_buf[c * 4], sizeof(float) * 3);
			a_buf[c * 4] = err;
			s -= err;
			float r = (float)fast_rand() * (1 / (float)UINT_MAX) -
			          (float)fast_rand() * (1 / (float)UINT_MAX);
			/*Clamp in float out of paranoia that the input will be >96 dBFS and wrap if the
			   integer is clamped. */
			int si;
			_o[pos + c] = si = (int)round(s + r);
			/*Including clipping in the noise shaping is generally disastrous:
			   the futile effort to restore the clipped energy results in more clipping.
			   However, small amounts-- at the level which could normally be created by
			   dither and rounding-- are harmless and can even reduce clipping somewhat
			   due to the clipping sometimes reducing the dither+rounding error. */
			b_buf[c * 4] = fmaxf(-1.5f, fminf(si - s, 1.5f));
		}
	}
}

long long audio_write(float *pcm, int frame_size, FILE *fout,
                      SpeexResamplerState *resampler, int *skip,
                      shapestate *shapemem, long long maxout)
{
	short *out = alloca(sizeof(short) * MAX_FRAME_SIZE * 2);
	float *buf = alloca(sizeof(float) * MAX_FRAME_SIZE * 2);
	
	int tmp_skip = 0;
	if (skip) {
		tmp_skip = (*skip > frame_size) ? (int)frame_size : *skip;
		*skip -= tmp_skip;
	}
	unsigned in_len = frame_size - tmp_skip;
	unsigned out_len = 1024 < maxout ? 1024 : maxout;
	speex_resampler_process_interleaved_float(resampler,
						  pcm + 2 * tmp_skip, &in_len,
						  buf, &out_len);
	/*Convert to short and save to output file */
	shape_dither_toshort(shapemem, out, buf, out_len);

	return fwrite((char *)out, 2 * 2, out_len < maxout ? out_len : maxout, fout);
}

void fwrite32(int i32, FILE *file)
{
	fwrite(&i32, 4, 1, file);
}

void fwrite16(int i16, FILE *file)
{
	fwrite(&i16, 2, 1, file);
}

void write_wav_header(FILE *file)
{
	fprintf(file, "RIFF");
	fwrite32(0x7fffffff, file);

	fprintf(file, "WAVEfmt ");
	fwrite32(16, file);
	fwrite16(1, file);
	fwrite16(2, file);
	fwrite32(44100, file);
	fwrite32(2 * 2 * 44100, file);
	fwrite16(2 * 2, file);
	fwrite16(16, file);

	fprintf(file, "data");
	fwrite32(0x7fffffff, file);
}

int main(int argc, char **argv)
{
	(void)argc;
	FILE *fout = NULL;
	struct OpusDecoder *st = opus_decoder_create();
	long long packet_count = 0;
	int total_links = 0;
	int stream_init = 0;
	ogg_int64_t page_granule = 0;
	ogg_int64_t link_out = 0;
	ogg_page og;
	ogg_packet op;
	ogg_stream_state os;
	ogg_int64_t audio_size = 0;
	int wav_format = 16;
	int preskip = 0;
	int gran_offset = 0;
	SpeexResamplerState *resampler = NULL;
	float *output = 0;
	shapestate shapemem;
	shapemem.a_buf = 0;
	shapemem.b_buf = 0;

	FILE *fin = fopen(argv[1], "rb");

	/* .opus files use the Ogg container to provide framing and timekeeping.
	 * http://tools.ietf.org/html/draft-terriberry-oggopus
	 * The easiest way to decode the Ogg container is to use libogg, so
	 *  thats what we do here.
	 * Using libogg is fairly straight forward-- you take your stream of bytes
	 *  and feed them to ogg_sync_ and it periodically returns Ogg pages, you
	 *  check if the pages belong to the stream you're decoding then you give
	 *  them to libogg and it gives you packets. You decode the packets. The
	 *  pages also provide timing information.*/
	ogg_sync_state oy;
	ogg_sync_init(&oy);

	/*Main decoding loop */
	while (!feof(fin)) {
		/*Get the ogg buffer for writing */
		char *data = ogg_sync_buffer(&oy, 200);
		/*Read bitstream from input file */
		int nb_read = fread(data, sizeof(char), 200, fin);
		ogg_sync_wrote(&oy, nb_read);

		/*Loop for all complete pages we got (most likely only one) */
		while (ogg_sync_pageout(&oy, &og) == 1) {
			if (stream_init == 0) {
				ogg_stream_init(&os, ogg_page_serialno(&og));
				stream_init = 1;
			}
			/*Add page to the bitstream */
			ogg_stream_pagein(&os, &og);
			page_granule = ogg_page_granulepos(&og);
			/*Extract all available packets */
			while (ogg_stream_packetout(&os, &op) == 1) {
				if (op.b_o_s) {
					link_out = 0;
					packet_count = 0;
					total_links++;
	                                preskip = op.packet[10] + (op.packet[11] << 8);

					/*Remember how many samples at the front we were told to skip
					   so that we can adjust the timestamp counting. */
					gran_offset = preskip;

					/*Setup the memory for the dithered output */
					shapemem.a_buf = calloc(2, sizeof(float) * 4);
					shapemem.b_buf = calloc(2, sizeof(float) * 4);

					output = malloc(sizeof(float) * MAX_FRAME_SIZE * 2);

					/*Normal players should just play at 48000 or their maximum rate,
					   as described in the OggOpus spec.  But for commandline tools
					   like opusdec it can be desirable to exactly preserve the original
					   sampling rate and duration, so we have a resampler here. */
					resampler = speex_resampler_init(2, 48000, 44100, 5, NULL);
					speex_resampler_skip_zeros(resampler);
					fout = fopen(argv[2], "wb");
					write_wav_header(fout);
				} else if (packet_count != 1) {
					/*Are we simulating loss for this packet? */
					int frame_size = opus_decoder_decode(st, (unsigned char *)op.packet,
					                                     op.bytes, output);

					/*This handles making sure that our output duration respects
					   the final end-trim by not letting the output sample count
					   get ahead of the granpos indicated value. */
					long long maxout = (page_granule - gran_offset) * 44100 / 48000 - link_out;
					long long outsamp = audio_write(output, frame_size, fout,
							                resampler, &preskip, &shapemem,
							                maxout);
					link_out += outsamp;
					audio_size += 2 * outsamp * 2;
				}
				packet_count++;
			}
			if (op.e_o_s) {
				float *zeros = calloc(100 * 2, sizeof(float));
				int drain = speex_resampler_get_input_latency(resampler);
				long long outsamp = audio_write(zeros, drain, fout, resampler, NULL, &shapemem,
						                (page_granule - gran_offset) * 44100 / 48000 - link_out);
				audio_size += 2 * outsamp * 2;
				free(zeros);
				speex_resampler_destroy(resampler);
				free(st);
			}
		}
	}

	fprintf(stderr, "\rDecoding complete.\n");
	fflush(stderr);

	/*If we were writing wav, go set the duration. */
	fseek(fout, 4, SEEK_SET);
	int tmp = audio_size + 20 + wav_format;
	fwrite(&tmp, 4, 1, fout);
	fseek(fout, 16 + wav_format, SEEK_CUR);
	tmp = audio_size;
	fwrite(&tmp, 4, 1, fout);

	ogg_stream_clear(&os);
	ogg_sync_clear(&oy);

	free(shapemem.a_buf);
	free(shapemem.b_buf);
	free(output);
	fclose(fin);
	fclose(fout);

	return 0;
}
