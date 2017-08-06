#ifndef OPUS_H
#define OPUS_H

#include <stddef.h>

/**
 * The Opus codec is designed for interactive speech and audio transmission over the Internet.
 * It is designed by the IETF Codec Working Group and incorporates technology from
 * Skype's SILK codec and Xiph.Org's CELT codec.
 *
 * The Opus codec is designed to handle a wide range of interactive audio applications,
 * including Voice over IP, videoconferencing, in-game chat, and even remote live music
 * performances. It can scale from low bit-rate narrowband speech to very high quality
 * stereo music. Its main features are:

 * @li Sampling rates from 8 to 48 kHz
 * @li Bit-rates from 6 kb/s to 510 kb/s
 * @li Support for both constant bit-rate (CBR) and variable bit-rate (VBR)
 * @li Audio bandwidth from narrowband to full-band
 * @li Support for speech and music
 * @li Support for mono and stereo
 * @li Support for multichannel (up to 255 channels)
 * @li Frame sizes from 2.5 ms to 60 ms
 * @li Good loss robustness and packet loss concealment (PLC)
 * @li Floating point and fixed-point implementation
 */

/**
 * @brief This page describes the process and functions used to decode Opus.
 *
 * The decoding process also starts with creating a decoder state.
 * The decoder state is always continuous in memory and only a shallow copy is sufficient
 * to copy it (e.g. memcpy())
 *
 * To decode a frame, opus_decoder_decode() must be called with a packet of compressed audio data:
 * @code
 * frame_size = opus_decoder_decode(st, packet, len, decoded);
 * @endcode
 * where
 *
 * @li packet is the byte array containing the compressed data
 * @li len is the exact number of bytes contained in the packet
 * @li decoded is the decoded audio data in float
 *
 * opus_decoder_decode() return the number of samples (per channel) decoded from the packet.
 *
 * Opus is a stateful codec with overlapping blocks and as a result Opus
 * packets are not coded independently of each other. Packets must be
 * passed into the decoder serially and in the correct order for a correct
 * decode. Lost packets can be replaced with loss concealment by calling
 * the decoder with a null pointer and zero length for the missing packet.
 *
 * A single codec state may only be accessed from a single thread at
 * a time and any required locking must be performed by the caller. Separate
 * streams must be decoded with separate decoder states and can be decoded
 * in parallel unless the library was compiled with NONTHREADSAFE_PSEUDOSTACK
 * defined.
 */

/**
 * Opus decoder state.
 * This contains the complete state of an Opus decoder.
 * It is position independent and can be freely copied.
 */
struct OpusDecoder;

/**
 * Decode a Opus packet with floating point output.
 * @param st Multistream decoder state.
 * @param data Input payload. Use a NULL pointer to indicate packet loss.
 * @param len Number of bytes in payload.
 * @param pcm Output signal, with interleaved samples.
 * @returns Number of samples decoded
 */
int opus_decoder_decode(struct OpusDecoder *st, const unsigned char *data,
		        size_t len, float *pcm);

struct OpusDecoder* opus_decoder_create();

#endif /* OPUS_H */
