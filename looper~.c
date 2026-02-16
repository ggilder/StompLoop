#include "m_pd.h"
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define DEFAULT_MAX_S 60
#define DECLICK_TIME_S 0.005
#define STATUS_UPDATE_MS 100
#define SPEED_SMOOTH_MS 50
#define OVERSAMPLE_FACTOR 2

typedef enum {
    LOOPER_STOPPED = 0,
    LOOPER_RECORDING,
    LOOPER_PLAYING
} t_looper_state;

typedef struct _looper {
    t_object x_obj;
    t_float f;

    t_sample *bufferL;
    t_sample *bufferR;
    size_t buffer_size;
    t_float read_pos;
    size_t loop_end;
    bool loop_end_set;

    t_float speed;
    t_float target_speed;
    t_float speed_increment;

    t_sample filter_state_l;
    t_sample filter_state_r;

    t_sample downsample_state_l;
    t_sample downsample_state_r;

    t_looper_state state;
    t_looper_state target_state;
    size_t fade_samples;
    size_t fade_pos;
    bool fading;
    bool fade_in; // true if fading in, false if fading out

    t_outlet *status_out;

    int counter;

    t_clock *status_clock;
} t_looper;

static t_class *looper_class;

void report_state(t_looper *x, t_looper_state state) {
    if (state == LOOPER_RECORDING || state == LOOPER_PLAYING) {
        // Compose list: symbol + time remaining
        t_atom out_list[2];
        t_symbol *state_sym = gensym(state == LOOPER_RECORDING ? "recording" : "playing");
        SETSYMBOL(&out_list[0], state_sym);
        t_float time_remaining = (t_float)(x->loop_end - x->read_pos) / (t_float)(sys_getsr() * OVERSAMPLE_FACTOR);
        if (x->speed != 0) time_remaining /= fabsf(x->speed);
        SETFLOAT(&out_list[1], time_remaining);
        outlet_list(x->status_out, &s_list, 2, out_list);
    } else {
        outlet_symbol(x->status_out, gensym("stopped"));
    }
}

void looper_tick(t_looper *x) {
    report_state(x, x->state);
    clock_delay(x->status_clock, STATUS_UPDATE_MS);
}

void *looper_new(t_symbol *s, int argc, t_atom *argv) {
    (void)s; // unused; silence warning
    t_looper *x = (t_looper *)pd_new(looper_class);

    t_float f = 0;
    if (argc > 0 && argv[0].a_type == A_FLOAT) {
        f = atom_getfloat(argv);
    }
    // Allocate buffer at oversampled rate
    x->buffer_size = sys_getsr() * OVERSAMPLE_FACTOR * ((f > 0) ? (size_t)f : DEFAULT_MAX_S);
    logpost(x, PD_DEBUG, "looper~: allocating %.2f seconds at %dx oversampling (%.2f MB per channel)",
            (t_float)x->buffer_size / (sys_getsr() * OVERSAMPLE_FACTOR),
            OVERSAMPLE_FACTOR,
            (t_float)(x->buffer_size * sizeof(t_sample)) / (1024.0 * 1024.0));
    x->bufferL = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    x->bufferR = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));

    if (!x->bufferL || !x->bufferR) {
        pd_error(x, "looper~: failed to allocate buffers");
        return NULL;
    }

    memset(x->bufferL, 0, x->buffer_size * sizeof(t_sample));
    memset(x->bufferR, 0, x->buffer_size * sizeof(t_sample));

    x->read_pos = 0;
    x->loop_end = x->buffer_size; // full buffer loop by default
    x->loop_end_set = false; // hasn't been set by the user

    x->speed = 1.0;
    x->target_speed = 1.0;
    x->speed_increment = 0;

    x->filter_state_l = 0;
    x->filter_state_r = 0;

    x->downsample_state_l = 0;
    x->downsample_state_r = 0;

    x->state = LOOPER_STOPPED;
    x->target_state = LOOPER_STOPPED;
    x->fade_samples = (size_t)(DECLICK_TIME_S * sys_getsr());
    x->fade_pos = 0;
    x->fading = false;
    x->fade_in = false;

    x->counter = 0;

    x->status_clock = clock_new(x, (t_method)looper_tick);
    // schedule first tick
    clock_delay(x->status_clock, STATUS_UPDATE_MS);

    inlet_new(&x->x_obj, &x->x_obj.ob_pd, &s_signal, &s_signal); // right inlet
    outlet_new(&x->x_obj, &s_signal); // left outlet
    outlet_new(&x->x_obj, &s_signal); // right outlet
    x->status_out = outlet_new(&x->x_obj, &s_symbol); // status outlet

    return (void *)x;
}

void looper_free(t_looper *x) {
    if (x->bufferL) freebytes(x->bufferL, x->buffer_size * sizeof(t_sample));
    if (x->bufferR) freebytes(x->bufferR, x->buffer_size * sizeof(t_sample));
    clock_free(x->status_clock);
}

void transition_state(t_looper *x, t_looper_state new_state) {
    if (x->state == new_state) {
        logpost(x, PD_ERROR, "looper: attempted to transition but already in state %d", new_state);
        x->fading = false;
        return; // no change
    }

    if (x->state == LOOPER_STOPPED) {
        if (new_state == LOOPER_PLAYING) {
            // fade in playback
            x->fade_in = true;
        } else if (new_state == LOOPER_RECORDING) {
            // fade in recording
            x->fade_in = true;
        } else {
            pd_error(x, "looper: invalid state transition from STOPPED to %d", new_state);
            return;
        }
    } else if (x->state == LOOPER_RECORDING) {
        // fade out recording
        x->fade_in = false;
    } else if (x->state == LOOPER_PLAYING) {
        if (new_state == LOOPER_RECORDING) {
            // fade in recording
            x->fade_in = true;
        } else if (new_state == LOOPER_STOPPED) {
            // fade out playback
            x->fade_in = false;
        } else {
            pd_error(x, "looper: invalid state transition from PLAYING to %d", new_state);
            return;
        }
    } else {
        pd_error(x, "looper: in unknown state %d", x->state);
        return;
    }

    x->target_state = new_state;
    report_state(x, new_state);
    x->fading = true;
    x->fade_pos = 0;
}

// control messages
void looper_start(t_looper *x) {
    if (x->state == LOOPER_STOPPED) {
        if (x->loop_end_set) {
            // Start playing existing loop
            transition_state(x, LOOPER_PLAYING);
            logpost(x, PD_NORMAL, "looper stopped -> playing");
        } else {
            // Start recording new loop
            transition_state(x, LOOPER_RECORDING);
            logpost(x, PD_NORMAL, "looper stopped -> recording");
        }
    } else if (x->state == LOOPER_RECORDING) {
        if (x->loop_end_set == false) {
            // Recording but haven't set loop end yet - mark end, we will now be overdubbing
            x->loop_end = (size_t)x->read_pos;
            x->loop_end_set = true;
            logpost(x, PD_NORMAL, "looper loop end set at %.2f seconds (%zu samples at %dx)",
                    x->loop_end / (sys_getsr() * OVERSAMPLE_FACTOR), x->loop_end, OVERSAMPLE_FACTOR);
            logpost(x , PD_NORMAL, "looper overdubbing");
        } else {
            // Overdubbing on existing loop, switch to playing
            transition_state(x, LOOPER_PLAYING);
            logpost(x, PD_NORMAL, "looper recording -> playing");
        }
    } else if (x->state == LOOPER_PLAYING) {
        // Already playing, switch to overdubbing
        transition_state(x, LOOPER_RECORDING);
        logpost(x, PD_NORMAL, "looper playing -> recording (overdubbing)");
    } else {
        pd_error(x, "looper: unknown state %d", x->state);
    }
}

void looper_stop(t_looper *x) {
    if (x->state == LOOPER_RECORDING) {
        if (x->loop_end_set == false) {
            // Stopping recording without having set loop end - set it now
            x->loop_end = (size_t)x->read_pos;
            x->loop_end_set = true;
            logpost(x, PD_NORMAL, "looper loop end set via stop at %.2f seconds (%zu samples at %dx)",
                    x->loop_end / (sys_getsr() * OVERSAMPLE_FACTOR), x->loop_end, OVERSAMPLE_FACTOR);
        }
    }
    transition_state(x, LOOPER_STOPPED);
    logpost(x, PD_NORMAL, "looper stopped");
}

void looper_clear(t_looper *x) {
    memset(x->bufferL, 0, x->buffer_size * sizeof(t_sample));
    memset(x->bufferR, 0, x->buffer_size * sizeof(t_sample));
    x->read_pos = 0;
    x->loop_end = x->buffer_size;
    x->loop_end_set = false;
    x->state = LOOPER_STOPPED;
    x->counter = 0;
    x->speed = 1.0;
    x->target_speed = 1.0;
    x->speed_increment = 0;
    x->filter_state_l = 0;
    x->filter_state_r = 0;
    x->downsample_state_l = 0;
    x->downsample_state_r = 0;
    logpost(x, PD_NORMAL, "looper cleared and stopped");
    report_state(x, x->state);
}

void looper_bang(t_looper *x) {
    report_state(x, x->state);
}

void looper_playpause(t_looper *x) {
    if (x->state == LOOPER_PLAYING || x->state == LOOPER_RECORDING) {
        looper_stop(x);
    } else if (x->loop_end_set) {
        looper_start(x);
    } else {
        logpost(x, PD_NORMAL, "looper playpause: no loop recorded yet, cannot start playing");
    }
}

void looper_speed(t_looper *x, t_floatarg f) {
    x->target_speed = f;
    // Calculate increment for smooth transition over SPEED_SMOOTH_MS
    t_float smooth_samples = (SPEED_SMOOTH_MS / 1000.0) * sys_getsr();
    x->speed_increment = (x->target_speed - x->speed) / smooth_samples;
    // Reset filter state on speed changes to prevent stuck filter artifacts
    x->filter_state_l = 0;
    x->filter_state_r = 0;
    logpost(x, PD_NORMAL, "looper speed set to %.2f", f);
}

// 4-point cubic (Hermite) interpolation
static inline t_sample hermite_interp(t_sample y0, t_sample y1, t_sample y2, t_sample y3, t_float frac) {
    t_sample c0 = y1;
    t_sample c1 = 0.5f * (y2 - y0);
    t_sample c2 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
    t_sample c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
    return ((c3 * frac + c2) * frac + c1) * frac + c0;
}

t_int *looper_perform(t_int *w) {
    t_looper *x = (t_looper *)(w[1]);
    t_sample *inL = (t_sample *)(w[2]);
    t_sample *inR = (t_sample *)(w[3]);
    t_sample *outL = (t_sample *)(w[4]);
    t_sample *outR = (t_sample *)(w[5]);
    int n = (int)(w[6]);

    // Process each external sample
    for (int i = 0; i < n; i++) {
        t_sample in_l = inL[i];
        t_sample in_r = inR[i];

        // Accumulator for downsampling
        t_sample accum_l = 0;
        t_sample accum_r = 0;

        // Process OVERSAMPLE_FACTOR internal samples per external sample
        for (int os = 0; os < OVERSAMPLE_FACTOR; os++) {
            // Smooth speed transitions (only on first oversample iteration)
            if (os == 0 && x->speed != x->target_speed) {
                if (fabsf(x->target_speed - x->speed) < fabsf(x->speed_increment)) {
                    x->speed = x->target_speed;
                } else {
                    x->speed += x->speed_increment;
                }
            }

            t_sample play_gain = 1.0;
            t_sample rec_gain = 1.0;

            // Handle state transitions with fading
            if (x->fading) {
                // Fade complete
                if (x->fade_pos >= x->fade_samples) {
                    x->fading = false;
                    x->state = x->target_state;
                    logpost(x, PD_DEBUG, "looper fade complete, new state %d", x->state);
                } else {
                    // Smoother fade using cosine curve
                    t_sample fade_ratio = 0.5 * (1 - cosf((t_float)x->fade_pos / (t_float)x->fade_samples * M_PI));
                    x->fade_pos++;
                    if (x->fade_in) {
                        if (x->target_state == LOOPER_RECORDING) {
                            rec_gain = fade_ratio;
                        } else if (x->target_state == LOOPER_PLAYING) {
                            play_gain = fade_ratio;
                        }
                    } else {
                        // Fade out
                        if (x->state == LOOPER_RECORDING) {
                            rec_gain = 1.0 - fade_ratio;
                            if (x->target_state == LOOPER_STOPPED) {
                                // When fading out from recording to stopped, we also need to fade out playback
                                play_gain = rec_gain;
                            }
                        } else if (x->state == LOOPER_PLAYING) {
                            play_gain = 1.0 - fade_ratio;
                        }
                    }
                }
            }

            if (x->state == LOOPER_STOPPED && !x->fading) {
                // If stopped, just output silence
                accum_l += 0;
                accum_r += 0;
                continue;
            }

            // Simple linear upsampling (hold input value across oversamples)
            t_sample os_in_l = in_l;
            t_sample os_in_r = in_r;

            // Apply anti-aliasing filter when recording at fast speeds
            if (x->state == LOOPER_RECORDING) {
                t_float abs_speed = fabsf(x->speed);
                if (abs_speed > 1.0f) {
                    t_float cutoff = 0.5f / abs_speed;
                    os_in_l = x->filter_state_l = x->filter_state_l + cutoff * (os_in_l - x->filter_state_l);
                    os_in_r = x->filter_state_r = x->filter_state_r + cutoff * (os_in_r - x->filter_state_r);
                }
            }

            // Get interpolated playback samples using cubic interpolation
            size_t pos_int = (size_t)x->read_pos;
            t_float frac = x->read_pos - pos_int;

            size_t i0 = (pos_int > 0) ? pos_int - 1 : x->loop_end - 1;
            size_t i1 = pos_int;
            size_t i2 = (pos_int + 1) % x->loop_end;
            size_t i3 = (pos_int + 2) % x->loop_end;

            t_sample play_l = hermite_interp(x->bufferL[i0], x->bufferL[i1],
                                             x->bufferL[i2], x->bufferL[i3], frac) * play_gain;
            t_sample play_r = hermite_interp(x->bufferR[i0], x->bufferR[i1],
                                             x->bufferR[i2], x->bufferR[i3], frac) * play_gain;

            if (x->state == LOOPER_RECORDING) {
                t_float abs_speed = fabsf(x->speed);
                size_t write_pos = (size_t)(x->read_pos + 0.5f);
                if (write_pos >= x->loop_end) write_pos = 0;
                t_float write_gain = rec_gain * abs_speed;
                x->bufferL[write_pos] += os_in_l * write_gain;
                x->bufferR[write_pos] += os_in_r * write_gain;
            }

            accum_l += play_l;
            accum_r += play_r;

            // Advance read position by speed (now in oversampled domain)
            x->read_pos += x->speed;

            while (x->read_pos >= x->loop_end) {
                x->read_pos -= x->loop_end;
                if (!x->loop_end_set) {
                    x->loop_end_set = true;
                    logpost(x, PD_NORMAL, "looper reached buffer end, loop end set at %.2f seconds (%zu samples at %dx)",
                            x->loop_end / (sys_getsr() * OVERSAMPLE_FACTOR), x->loop_end, OVERSAMPLE_FACTOR);
                }
            }
            while (x->read_pos < 0) {
                x->read_pos += x->loop_end;
            }
        }

        // Downsample with averaging and lowpass
        t_sample raw_l = accum_l / (t_float)OVERSAMPLE_FACTOR;
        t_sample raw_r = accum_r / (t_float)OVERSAMPLE_FACTOR;

        t_float ds_cutoff = 0.45f;
        outL[i] = x->downsample_state_l = x->downsample_state_l + ds_cutoff * (raw_l - x->downsample_state_l);
        outR[i] = x->downsample_state_r = x->downsample_state_r + ds_cutoff * (raw_r - x->downsample_state_r);
    }
    return (w + 7);
}

void looper_dsp(t_looper *x, t_signal **sp) {
    dsp_add(looper_perform, 6,
            x,
            sp[0]->s_vec, // inL
            sp[1]->s_vec, // inR
            sp[2]->s_vec, // outL
            sp[3]->s_vec, // outR
            sp[0]->s_n);
}

void looper_tilde_setup(void) {
    looper_class = class_new(gensym("looper~"),
                             (t_newmethod)looper_new,
                             (t_method)looper_free,
                             sizeof(t_looper),
                             CLASS_DEFAULT,
                             0);

    CLASS_MAINSIGNALIN(looper_class, t_looper, f);
    class_addmethod(looper_class, (t_method)looper_dsp, gensym("dsp"), A_CANT, 0);
    class_addmethod(looper_class, (t_method)looper_start, gensym("start"), 0);
    class_addmethod(looper_class, (t_method)looper_stop, gensym("stop"), 0);
    class_addmethod(looper_class, (t_method)looper_clear, gensym("clear"), 0);
    class_addmethod(looper_class, (t_method)looper_playpause, gensym("playpause"), 0);
    class_addmethod(looper_class, (t_method)looper_speed, gensym("speed"), A_FLOAT, 0);
    class_addbang(looper_class, (t_method)looper_bang);
}

