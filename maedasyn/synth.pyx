## -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 15:01:24 2013

@author: Ronald L. Sprouse (ronald@berkeley.edu)
"""

import numpy as np
cimport numpy as np
cimport maedasyn.maedasyn as ms


########################################################################
#
#
# C-style API
# This section provides more-or-less direct access to C code internals.
# You are encouraged to use the Python classes instead of direct access.
#
#
#########################################################################
AMloc = ms.AMloc
AMnum = ms.AMnum
AP = ms.AP
NPAR = ms.NPAR
F0_LOC = ms.F0_LOC
TIME = ms.TIME
FRAME_DUR = ms.FRAME_DUR
smpfrq = ms.smpfrq

NP = ms.NP

# Wrapper for C code synth_frame() in synthesize.c.
def synth_frame(
    np.ndarray[float, ndim=1, mode="c"] params not None,
    np.ndarray[short, ndim=1, mode="c"] buff not None,
    mode
):
    ms.synth_frame(&params[0], &buff[0], mode)
    return None
 
###### Access to C arrays #####
#
#
#  We copy values from C arrays. This is a quick and easy way to get the job done.
#  It would be better to import them more directly through their pointers.
#
#
################################

def get_ivt():
    ivt = []
    for idx in np.arange(ms.NP):
        ivt.append( {"x": ms.ivt[idx].x, "y": ms.ivt[idx].y} )
    return ivt

def get_evt():
    evt = []
    for idx in np.arange(ms.NP):
        evt.append( {"x": ms.evt[idx].x, "y": ms.evt[idx].y} )
    return evt

def get_afvt():
    afvt = []
    for idx in np.arange(ms.nss):
        afvt.append( {"A": ms.afvt[idx].A, "x": ms.afvt[idx].x} )
    return afvt

def get_afnt():
    afnt = []
    for idx in np.arange(13):
        afnt.append( {"A": ms.afnt[idx].A, "x": ms.afnt[idx].x} )
    return afnt

def get_alph():
    alph = []
    for idx in np.arange(ms.M4):
        alph.append(ms.alph[idx])
    return alph

def set_alph(alph):
    for idx in np.arange(ms.M4):
        ms.alph[idx] = alph[idx]
    
def get_beta():
    beta = []
    for idx in np.arange(ms.M4):
        beta.append(ms.beta[idx])
    return beta

def set_beta(beta):
    for idx in np.arange(ms.M4):
        ms.beta[idx] = beta[idx]
    
def get_u_wal():
    u_wal = []
    for idx in np.arange(ms.NVRS_WAL):
        u_wal.append(ms.u_wal[idx])
    return u_wal

def set_u_wal(u_wal):
    for idx in np.arange(ms.NVRS_WAL):
        ms.u_wal[idx] = u_wal[idx]

###### End of C array access #####

###### End of C-style API #####

cdef class Synth(object):
    '''The synthesizer object.'''
    cdef public int _bufsize
    cdef np.ndarray _buffer
    property rate:
        def __get__(self):
            return ms.smpfrq
    property buffer:
        def __get__(self):
            return self._buffer
    property ivt:
        def __get__(self):
            return get_ivt()
    property ivt_x:
        def __get__(self):
            return [ms.ivt[idx].x for idx in np.arange(ms.NP)]
    property ivt_y:
        def __get__(self):
            return [ms.ivt[idx].y for idx in np.arange(ms.NP)]
    property evt:
        def __get__(self):
            return get_evt()
    property evt_x:
        def __get__(self):
            return [ms.evt[idx].x for idx in np.arange(ms.NP)]
    property evt_y:
        def __get__(self):
            return [ms.evt[idx].y for idx in np.arange(ms.NP)]
    property afvt:
        def __get__(self):
            return get_afvt()
    property afvt_A:
        def __get__(self):
            return [ms.afvt[idx].A for idx in np.arange(ms.nss)]
    property afvt_x:
        def __get__(self):
            return [ms.afvt[idx].x for idx in np.arange(ms.nss)]
    property afnt:
        def __get__(self):
            return get_afnt()
    property afnt_A:
        def __get__(self):
            return [ms.afnt[idx].A for idx in np.arange(13)]
    property afnt_x:
        def __get__(self):
            return [ms.afnt[idx].x for idx in np.arange(13)]
    property alph:
        def __get__(self):
            return get_alph()
        def __set__(self, vals):
            assert(len(vals) == len(get_alph()))
            set_alph(vals)
    property beta:
        def __get__(self):
            return get_beta()
        def __set__(self, vals):
            assert(len(vals) == len(get_beta()))
            set_beta(vals)
    property u_wal:
        def __get__(self):
            return get_u_wal()
        def __set__(self, vals):
            assert(len(vals) == len(get_u_wal()))
            set_u_wal(vals)

    def __cinit__(self):
        self._bufsize = np.floor(ms.FRAME_DUR * ms.smpfrq)
        self._buffer = np.zeros(self._bufsize, dtype=np.int16)
        self.synthesize(FrameParam(), 1)  # Initialize

    def synthesize(self, params, mode):
        #self._buffer = np.zeros(self._bufsize * mode, dtype=np.int16)
        synth_frame(params.as_ndarray(), self._buffer, mode)

    def time_for_frameidx(self, idx):
        '''Calculate value of time from a frame index.'''
        return (idx * ms.FRAME_DUR) * 1000

        
# You can override one or more defaults identified by name, e.g.
#
# my_uw = FrameParam('uw', lip_aperture=1.2)
#
cdef class FrameParam(object):
    '''FrameParam parameters for Maeda synthesizer.'''
    cdef public float time
    cdef public float f0
    cdef public float glottal_aperture
    cdef public float jaw_position
    cdef public float dorsum_position
    cdef public float dorsum_shape
    cdef public float apex_position
    cdef public float lip_aperture
    cdef public float lip_protrusion
    cdef public float larynx_height
    cdef public float nasal_coupling

    # elements in Maeda model specifications for vowels
    fields = (
        'time',
        'f0',
        'glottal_aperture',
        'jaw_position',
        'dorsum_position',
        'dorsum_shape',
        'apex_position',
        'lip_aperture',
        'lip_protrusion',
        'larynx_height',
        'nasal_coupling'
    )

    # These are fields for which we have defaults for various vowels.
    vowel_fields = fields[3:]

    # Construct a FrameParam. If name is given, use default specifications for various
    # vowels. You can specifiy the value of parameters with keyword arguments, the names
    # of which are the values of FrameParam.fields. These keyword arguments override defaults
    # provided by name.
    #
    # uw = FrameParam('uw')
    # uw = FrameParam('uw', lip_aperture=1.3)   # override 'uw' lip_aperture
    def __init__(self, name=None, **kwargs):

        # Default to 0.0 values
        for idx,articulator in enumerate(FrameParam.fields):
            setattr(self, articulator, 0.0)

        # Overwrite defaults by cloning from another FrameParam or from default specifications.
        if type(name) is FrameParam:
            for idx,field in enumerate(FrameParam.fields):
                setattr(self, field, getattr(name, field))
        elif name in ('aa', 'uw', 'iy', 'ey', 'eh', 'ah', 'ao', 'ow', 'iw', 'ew', 'oe'):
            if name == 'aa':
                aspec = &ms.aa[0]
            elif name == 'uw':
                aspec = &ms.uw[0]
            elif name == 'iy':
                aspec = &ms.iy[0]
            elif name == 'ey':
                aspec = &ms.ey[0]
            elif name == 'eh':
                aspec = &ms.eh[0]
            elif name == 'ah':
                aspec = &ms.ah[0]
            elif name == 'ao':
                aspec = &ms.ao[0]
            elif name == 'ow':
                aspec = &ms.ow[0]
            elif name == 'iw':
                aspec = &ms.iw[0]
            elif name == 'ew':
                aspec = &ms.ew[0]
            elif name == 'oe':
                aspec = &ms.oe[0]
    
            for idx,articulator in enumerate(FrameParam.vowel_fields):
                setattr(self, articulator, aspec[idx])

        # Finally, overwrite existing values with keyword values.
        for k,v in kwargs.iteritems():
            setattr(self, k, v)

    def as_ndarray(self):
        '''Convert parameter attributes to an ndarray that the C code can use.'''
        a = np.ndarray([len(self.fields)], dtype=np.float32)
        for idx,field in enumerate(FrameParam.fields):
            a[idx] = getattr(self, field)
        return a

    # Override addition. If we are adding another FrameParam, add the corresponding
    # field. If we are adding a number, add that number to all the fields.
    def __add__(self, other):
        result = FrameParam(self)
        if type(other) is FrameParam:
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) + getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) + other)
        return result

    # Override subtraction. If we are subtracting another FrameParam, subtract the corresponding
    # field. If we are subtracting a number, subtract that number from all the fields.
    def __sub__(self, other):
        result = FrameParam(self)
        if type(other) is FrameParam:
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) - getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) - other)
        return result

    # Override multiplication. If we are multiplying another FrameParam, multiply the corresponding
    # field. If we are multiplying by a number, multiply all the fields by that number.
    def __mul__(self, other):
        result = FrameParam(self)
        if type(other) is FrameParam:
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) * getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) * other)
        return result

    # Override division. If we are dividing another FrameParam, divide the corresponding
    # field. If we are dividing by a number, divide all the fields by that number.
    def __div__(self, other):
        result = FrameParam(self)
        if type(other) is FrameParam:
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) / getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FrameParam.fields):
                setattr(result, field, getattr(self, field) / other)
        return result


