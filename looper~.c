#include "m_pd.h"
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define DEFAULT_MAX_S 60
#define DECLICK_TIME_S 0.005
#define STATUS_UPDATE_MS 100

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
    size_t pos;
    size_t loop_end;
    bool loop_end_set;

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
        t_float time_remaining = (t_float)(x->loop_end - x->pos) / (t_float)sys_getsr();
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
    x->buffer_size = sys_getsr() * ((f > 0) ? (size_t)f : DEFAULT_MAX_S);
    logpost(x, PD_DEBUG, "looper~: allocating %.2f seconds (%.2f MB per channel) buffer",
            (t_float)x->buffer_size / sys_getsr(),
            (t_float)(x->buffer_size * sizeof(t_sample)) / (1024.0 * 1024.0));
    x->bufferL = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));
    x->bufferR = (t_sample *)getbytes(x->buffer_size * sizeof(t_sample));

    if (!x->bufferL || !x->bufferR) {
        pd_error(x, "looper~: failed to allocate buffers");
        return NULL;
    }

    memset(x->bufferL, 0, x->buffer_size * sizeof(t_sample));
    memset(x->bufferR, 0, x->buffer_size * sizeof(t_sample));

    x->pos = 0;
    x->loop_end = x->buffer_size; // full buffer loop by default
    x->loop_end_set = false; // hasn't been set by the user

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
            x->loop_end = x->pos;
            x->loop_end_set = true;
            logpost(x, PD_NORMAL, "looper loop end set at %.2f seconds (%zu samples)", x->loop_end / sys_getsr(), x->loop_end);
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
            x->loop_end = x->pos;
            x->loop_end_set = true;
            logpost(x, PD_NORMAL, "looper loop end set via stop at %.2f seconds (%zu samples)", x->loop_end / sys_getsr(), x->loop_end);
        }
    }
    transition_state(x, LOOPER_STOPPED);
    logpost(x, PD_NORMAL, "looper stopped");
}

void looper_clear(t_looper *x) {
    memset(x->bufferL, 0, x->buffer_size * sizeof(t_sample));
    memset(x->bufferR, 0, x->buffer_size * sizeof(t_sample));
    x->pos = 0;
    x->loop_end = x->buffer_size;
    x->loop_end_set = false;
    x->state = LOOPER_STOPPED;
    x->counter = 0;
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

t_int *looper_perform(t_int *w) {
    t_looper *x = (t_looper *)(w[1]);
    t_sample *inL = (t_sample *)(w[2]);
    t_sample *inR = (t_sample *)(w[3]);
    t_sample *outL = (t_sample *)(w[4]);
    t_sample *outR = (t_sample *)(w[5]);
    int n = (int)(w[6]);

    x->counter += n;
    // Debug logging every 2 seconds
    /* if (x->counter >= 2 * sys_getsr()) { */
    /*     logpost(x, PD_DEBUG, "Hello from looper, state=%d, pos=%zu, loop_end=%zu, loop_end_set=%d", */
    /*             x->state, x->pos, x->loop_end, x->loop_end_set); */
    /*     if (x->fading) { */
    /*         logpost(x, PD_DEBUG, "  fading %s, fade_pos=%zu/%zu, target_state=%d", */
    /*                 x->fade_in ? "in" : "out", x->fade_pos, x->fade_samples, x->target_state); */
    /*     } */
    /*     x->counter = 0; */
    /* } */

    for (int i = 0; i < n; i++) {
        t_sample play_gain = 1.0;
        t_sample rec_gain = 1.0;

        // Handle state transitions with fading
        if (x->fading) {
            if (x->fade_pos >= x->fade_samples) {
                // Fade complete
                x->fading = false;
                x->state = x->target_state;
                logpost(x, PD_DEBUG, "looper fade complete, new state %d", x->state);
            } else {
                // smoother cosine fade
                t_sample fade_ratio = 0.5 * (1 - cosf((t_float)x->fade_pos / (t_float)x->fade_samples * M_PI));
                x->fade_pos++;
                if (x->fade_in) {
                    if (x->target_state == LOOPER_RECORDING) {
                        rec_gain = fade_ratio;
                    } else if (x->target_state == LOOPER_PLAYING) {
                        play_gain = fade_ratio;
                    }
                } else {
                    // Fading out
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
            outL[i] = 0;
            outR[i] = 0;
            continue;
        }

        t_sample in_l = inL[i];
        t_sample in_r = inR[i];

        t_sample play_l = x->bufferL[x->pos] * play_gain;
        t_sample play_r = x->bufferR[x->pos] * play_gain;

        if (x->state == LOOPER_RECORDING) {
            // Overdub: mix input with existing buffer content
            // TODO: manage levels with compression?
            x->bufferL[x->pos] = play_l + (in_l * rec_gain);
            x->bufferR[x->pos] = play_r + (in_r * rec_gain);
        }

        outL[i] = play_l;
        outR[i] = play_r;

        x->pos = x->pos + 1;
        if (x->pos >= x->loop_end) {
            x->pos = x->pos % x->loop_end; // wrap around to loop start
            if (!x->loop_end_set) {
                // Mark that we've set loop end. This should only happen if we hit the max buffer size
                x->loop_end_set = true;
                logpost(x, PD_NORMAL, "looper reached buffer end, loop end set at %.2f seconds (%zu samples)", x->loop_end / sys_getsr(), x->loop_end);
            }
        }
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
    class_addbang(looper_class, (t_method)looper_bang);
}

