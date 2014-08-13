# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 14:44:19 2013

@author: Ronald L. Sprouse (ronald@berkeley.edu)
"""

cdef extern from '../c/always.h':
    pass

cdef extern from '../c/vtconfig.h':
    pass

cdef extern from '../c/lam_lib.h':
    int NP
    ctypedef struct float2D:
        short x
        short y

    float2D *ivt
    float2D *evt
    float2D *f2d

cdef extern from '../c/vsyn_lib.h':
    pass

cdef extern from '../c/vtt_lib.c':
    pass

cdef extern from '../c/lam_lib.c':
    pass

cdef extern from '../c/vsyn_lib.c':
    pass

cdef extern from '../c/synthesize.c':
    void synth_frame(float *params, short *buffer, short mode)
    int AMloc
    int AMnum
    int AP
    int NPAR
    int F0_LOC
    int TIME
    float FRAME_DUR
    float smpfrq
    float *aa
    float *uw
    float *iy
    float *ey
    float *eh
    float *ah
    float *ao
    float *ow
    float *iw
    float *ew
    float *oe
