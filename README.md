# rumpus

Opus decoder in rust

## test

    ./clean.sh; ./make.sh; ./opusdec music.opus music.wav; util/wavdiff music.libopus.wav music.wav

## API

    struct OpusDecoder* opus_decoder_create();

    int opus_decoder_decode(struct OpusDecoder *st, const unsigned char *data,
                            size_t len, float *pcm);

