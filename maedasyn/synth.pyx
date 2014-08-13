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
    for idx in np.arange(ms.NPAR):
        ivt.append( {"x": ms.ivt[idx].x, "y": ms.ivt[idx].y} )
    return ivt

def get_evt():
    evt = []
    for idx in np.arange(ms.NPAR):
        evt.append( {"x": ms.ivt[idx].x, "y": ms.ivt[idx].y} )
    return evt

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
    property evt:
        def __get__(self):
            return get_evt()

    def __cinit__(self):
        self._bufsize = np.floor(ms.FRAME_DUR * ms.smpfrq)
        self._buffer = np.zeros(self._bufsize, dtype=np.int16)
        self.synthesize(FParam(), 1)  # Initialize

    def synthesize(self, params, mode):
        synth_frame(params.as_ndarray(), self._buffer, mode)

    def time_for_frameidx(self, idx):
        '''Calculate value of time from a frame index.'''
        return (idx * ms.FRAME_DUR) * 1000

        
# You can override one or more defaults identified by name, e.g.
#
# my_uw = FParam('uw', lip_aperture=1.2)
#
cdef class FParam(object):
    '''FParam parameters for Maeda synthesizer.'''
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

    # Construct a FParam. If name is given, use default specifications for various
    # vowels. You can specifiy the value of parameters with keyword arguments, the names
    # of which are the values of FParam.fields. These keyword arguments override defaults
    # provided by name.
    #
    # uw = FParam('uw')
    # uw = FParam('uw', lip_aperture=1.3)   # override 'uw' lip_aperture
    def __init__(self, name=None, **kwargs):

        # Default to 0.0 values
        for idx,articulator in enumerate(FParam.fields):
            setattr(self, articulator, 0.0)

        # Overwrite defaults by cloning from another FParam or from default specifications.
        if type(name) is FParam:
            for idx,field in enumerate(FParam.fields):
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
    
            for idx,articulator in enumerate(FParam.vowel_fields):
                setattr(self, articulator, aspec[idx])

        # Finally, overwrite existing values with keyword values.
        for k,v in kwargs.iteritems():
            setattr(self, k, v)

    def as_ndarray(self):
        '''Convert parameter attributes to an ndarray that the C code can use.'''
        a = np.ndarray([len(self.fields)], dtype=np.float32)
        for idx,field in enumerate(FParam.fields):
            a[idx] = getattr(self, field)
        return a

    # Override addition. If we are adding another FParam, add the corresponding
    # field. If we are adding a number, add that number to all the fields.
    def __add__(self, other):
        result = FParam(self)
        if type(other) is FParam:
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) + getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) + other)
        return result

    # Override subtraction. If we are subtracting another FParam, subtract the corresponding
    # field. If we are subtracting a number, subtract that number from all the fields.
    def __sub__(self, other):
        result = FParam(self)
        if type(other) is FParam:
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) - getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) - other)
        return result

    # Override multiplication. If we are multiplying another FParam, multiply the corresponding
    # field. If we are multiplying by a number, multiply all the fields by that number.
    def __mul__(self, other):
        result = FParam(self)
        if type(other) is FParam:
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) * getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) * other)
        return result

    # Override division. If we are dividing another FParam, divide the corresponding
    # field. If we are dividing by a number, divide all the fields by that number.
    def __div__(self, other):
        result = FParam(self)
        if type(other) is FParam:
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) / getattr(other, field))
        else:    # should be a number
            for idx,field in enumerate(FParam.fields):
                setattr(result, field, getattr(self, field) / other)
        return result



def synth_to_file(fname):
    bufsize = np.floor(ms.FRAME_DUR * ms.smpfrq)
    final_buffer = np.zeros(bufsize*10, dtype=np.int16)
    mybuffer = np.zeros(bufsize, dtype=np.int16)
    parameters = np.zeros(ms.NPAR, dtype=np.float32)
    target1 = np.zeros(ms.NPAR, dtype=np.float32)
    target2 = np.zeros(ms.NPAR, dtype=np.float32)
    d = np.zeros(ms.NPAR, dtype=np.float32)

    target1[ms.F0_LOC] = 130
    target1[ms.AP] = 0.2
    target2[ms.F0_LOC] = 100
    target2[ms.AP] = 0.2
    for j in np.arange(ms.AMnum):
        target1[ms.AMloc+j] = ms.uw[j]
        target2[ms.AMloc+j] = ms.iy[j]

    n = 80-20
    for j in np.arange(ms.NPAR):
        d[j] = (target2[j] - target1[j])/n
        parameters[j] = target1[j]

    # Initialize the synthesizer.
    mode = 1
    synth_frame(parameters, mybuffer, mode)

    mode = 2
    with open(fname, 'wb') as soundfile:
    #    /* construct a par array for testing */
        for i in np.arange(100):
            parameters[ms.TIME] = (i*ms.FRAME_DUR)*1000
            if i >= 80:
                for j in np.arange(ms.NPAR):
                    parameters[j] = target2[j]
            elif i <= 20:
                for j in np.arange(ms.NPAR):
                    parameters[j] = target1[j]
            else:
                for j in np.arange(ms.NPAR):
                    parameters[j] = target1[j] + d[j]*(i-20.0)
            synth_frame(parameters, mybuffer, mode)
            soundfile.write(mybuffer);
        synth_frame(parameters, final_buffer, bufsize*10)
        soundfile.write(final_buffer)

#####UNDO#####bufsize = np.floor(ms.FRAME_DUR * ms.smpfrq)
#####UNDO######mybuffer = np.zeros(bufsize, dtype=SHORT_t)
#####UNDO#####cdef np.ndarray final_buffer = np.zeros(bufsize*10, dtype=np.int)
#####UNDO#####mybuffer = np.zeros(bufsize, dtype=np.int)
#####UNDO#####parameters = np.zeros(ms.NPAR, dtype=np.float)
#####UNDO#####cdef np.ndarray target1 = np.zeros(ms.NPAR, dtype=np.float)
#####UNDO#####cdef np.ndarray target2 = np.zeros(ms.NPAR, dtype=np.float)
#####UNDO#####cdef short num_written
#####UNDO#####cdef int i,j,n
#####UNDO#####cdef np.ndarray d = np.zeros(ms.NPAR, dtype=np.float)
#####UNDO#####
#####UNDO#####target1[ms.F0_LOC] = 130
#####UNDO#####target1[ms.AP] = 0.2
#####UNDO#####target2[ms.F0_LOC] = 100
#####UNDO#####target2[ms.AP] = 0.2
#####UNDO#####for j in np.arange(ms.AMnum):
#####UNDO#####    target1[ms.AMloc+j] = ms.uw[j]
#####UNDO#####    target2[ms.AMloc+j] = ms.iy[j]
#####UNDO#####
#####UNDO#####n = 80-20
#####UNDO#####for j in np.arange(ms.NPAR):
#####UNDO#####    d[j] = (target2[j] - target1[j])/n
#####UNDO#####    parameters[j] = target1[j]
#####UNDO#####
#####UNDO###### Initialize the synthesizer.
#####UNDO#####cdef int mode = 1
#####UNDO#####synth_frame(parameters, mybuffer, mode)
#####UNDO#####
#####UNDO#####mode = 2
#####UNDO#####with open('../resources/mtest2.raw', 'wb') as soundfile:
#####UNDO######    /* construct a par array for testing */
#####UNDO#####    for i in np.arange(100):
#####UNDO#####        parameters[ms.TIME] = (i*ms.FRAME_DUR)*1000
#####UNDO#####        if i >= 80:
#####UNDO#####            for j in np.arange(ms.NPAR):
#####UNDO#####                parameters[j] = target2[j]
#####UNDO#####        elif i <= 20:
#####UNDO#####            for j in np.arange(ms.NPAR):
#####UNDO#####                parameters[j] = target1[j]
#####UNDO#####        else:
#####UNDO#####            for j in np.arange(ms.NPAR):
#####UNDO#####                parameters[j] = target1[j] + d[j]*(i-20.0)
#####UNDO#####        synth_frame(parameters, mybuffer, mode);
#        
#        /* at this point in the calling code you can refer to evt[i].x and evt[i].y to 
#            draw the "passive" surface of the vocal tract, and ivt[i].x and ivt[i].y to
#            draw the "active" surface.
#         
#         
#        for (j=0;j<NP;j++) {
#            printf("%d\t%3.4f\t%3.4f\n",j,evt[j].x,evt[j].y);
#            
#        }
#        */
#         
#        if ((num_written = fwrite( buffer, sizeof(short), bufsize, soundfile )) != bufsize) {
#            printf("%s\n","write to file failed");
#        }
#    }
#    synth_frame(parameters,final_buffer,bufsize*10);
#    if ((num_written = fwrite( final_buffer, sizeof(short), bufsize*10, soundfile )) != bufsize*10) {
#        printf("%s\n","write to file failed");
#    }

#int main() {
#    
#    short bufsize = (short)(FRAME_DUR*smpfrq);
#    short buffer[bufsize];
#    short final_buffer[bufsize*10];
#    float parameters[NPAR],target1[NPAR],target2[NPAR];
#    short num_written;
#    int i,j,n;
#    float d[NPAR];
#    FILE *soundfile;
#    
#    
#    if( (soundfile = fopen("/Users/transferredkj/Google Drive/maeda/mtest2.raw", "w+b")) == NULL ) {
#        fprintf(stderr,"Can't open file.");
#        exit(-1);
#	}
#    
#    target1[F0_LOC] = 130;
#    target2[F0_LOC] = 100;
#    target1[AP]=0.2;
#    target2[AP]=0.2;
#    for (j=0; j<AMnum; j++) target1[AMloc+j] = uw[j];
#    for (j=0; j<AMnum; j++) target2[AMloc+j] = iy[j];
#
#    n= 80-20;
#    for (j=1; j<NPAR; j++) {
#        d[j] = (target2[j] - target1[j])/n;  // delta values for interpolation
#        parameters[j] = target1[j];  // starting params for initialization
#    }
#    
#    synth_frame(parameters, buffer, 1);  // mode 1 = initialize the synthesizer
#
#    
#    /* construct a par array for testing */
#    for (i=0; i<100; i++) {
#        parameters[TIME] = (i*FRAME_DUR)*1000;
#        if (i>=80) {for (j=0; j<NPAR; j++) parameters[j] = target2[j];}
#        if (i<=20) {for (j=0; j<NPAR; j++) parameters[j] = target1[j];}
#        if (i>20 && i<80) {
#            for (j=0; j<NPAR; j++) parameters[j] = target1[j] + d[j]*(i-20.0);
#        }
#        synth_frame(parameters, buffer, 2);
#        
#        /* at this point in the calling code you can refer to evt[i].x and evt[i].y to 
#            draw the "passive" surface of the vocal tract, and ivt[i].x and ivt[i].y to
#            draw the "active" surface.
#         
#         
#        for (j=0;j<NP;j++) {
#            printf("%d\t%3.4f\t%3.4f\n",j,evt[j].x,evt[j].y);
#            
#        }
#        */
#         
#        if ((num_written = fwrite( buffer, sizeof(short), bufsize, soundfile )) != bufsize) {
#            printf("%s\n","write to file failed");
#        }
#    }
#    synth_frame(parameters,final_buffer,bufsize*10);
#    if ((num_written = fwrite( final_buffer, sizeof(short), bufsize*10, soundfile )) != bufsize*10) {
#        printf("%s\n","write to file failed");
#    }
#    
#    fclose(soundfile);
#
#    
#}

## These map param names to index locations in defval. They are defined in the
## symb1 and symb2 arrays in the C code. This is much easier to read.
#params_map = {
#    'sr': 0,  # C
#    'nf': 1,  # C
#    'du': 2,  # C
#    'ss': 3,  # C
#    'ui': 4,  # C
#    'rs': 5,  # C
#    'f0': 6,
#    'av': 7,
#    'F1': 8,
#    'b1': 9,
#    'F2': 10,
#    'b2': 11,
#    'F3': 12,
#    'b3': 13,
#    'F4': 14,
#    'b4': 15,
#    'F5': 16,
#    'b5': 17,
#    'f6': 18,
#    'b6': 19,
#    'fz': 20,
#    'bz': 21,
#    'fp': 22,
#    'bp': 23,
#    'ah': 24,
#    'oq': 25,
#    'at': 26,
#    'tl': 27,
#    'af': 28,
#    'sk': 29,
#    'a1': 30,
#    'p1': 31,
#    'a2': 32,
#    'p2': 33,
#    'a3': 34,
#    'p3': 35,
#    'a4': 36,
#    'p4': 37,
#    'a5': 38,
#    'p5': 39,
#    'a6': 40,
#    'p6': 41,
#    'an': 42,
#    'ab': 43,
#    'ap': 44,
#    'os': 45,  # C
#    'g0': 46,
#    'dF': 47,
#    'db': 48
#}
#
## These are other parameters that are not in defval and not in Klatt's original
## .doc file format. They are present in the .klp format.
#extra_params = (
#    'agc',  # automatic gain control mode (formerly the -g command line switch)
#            # comma is necessary to maintain this as a single-valued tuple
#)
#
#class synthesizer(object):
#    def __init__(self):
#        ''' Initialize global variables as done in C code's main() function.'''
#        # This section is from main().
#        for idx in range(klatt_defs.NPAR):
#            klatt_defs.defval[idx] = klatt_defs.cdefval[idx]
#        #self.update_constant_params_from_defval()  # seems to be okay without this
#        klatt_defs.initpars()
#
#    def update_constant_params_from_defval(self, verbose=0):
#        ''' Set C variables for constant params from current values in defval[]. '''
#        klatt_defs.OUTSELECT = klatt_defs.defval[params_map['os']]
#        klatt_defs.RANSEED = klatt_defs.defval[params_map['rs']]
#        klatt_defs.NFCASC = klatt_defs.defval[params_map['nf']]
#        klatt_defs.SOURCE_SELECT = klatt_defs.defval[params_map['ss']]
#        klatt_defs.setlimits(verbose)
#
#    def get_defval(self, param):
#        ''' Get a parameter value from klatt_defs.defval[]. '''
#        return klatt_defs.defval[params_map[param]]
#
#    def get_constant_params(self):
#        ''' Return constant params in a dict. '''
#        params = {}
#        for idx in range(klatt_defs.NPAR):
#            param = params_map.keys()[params_map.values().index(idx)]
#            if klatt_defs.cv[idx] != klatt_defs.VARRIED:
#                params[param] = klatt_defs.defval[idx]
#        for param in extra_params:
#            if param == 'agc':
#                params[param] = klatt_defs.gain_control
#        return params
#
#    def get_varied_params(self):
#        ''' Return varied params in a dict. '''
#        params = {}
#        for idx in range(klatt_defs.NPAR):
#            param = params_map.keys()[params_map.values().index(idx)]
#            if klatt_defs.cv[idx] == klatt_defs.VARRIED:
#                val = np.zeros([klatt_defs.nframes], dtype=np.int16)
#                for nf in range(len(val)):
#                    val[nf] = klatt_defs.pdata[idx][nf]
#                params[param] = val
#        return params
#
#    def get_params(self):
#        ''' Return all params in a dict. '''
#        return dict(
#                   self.get_constant_params().items() +
#                   self.get_varied_params().items()
#        )
#
#    def set_params(self, params=None):
#        ''' Set C variables based on params. If a param has only a single value, it is treated as a constant variable. If its value is list-like it is treated as a varied parameter. '''
#
#        for idx in range(klatt_defs.NPAR):
#            klatt_defs.clearpar(idx)
#        klatt_defs.ipsw = klatt_defs.MYTRUE
#
#        num_varied = 0
#        for (param, val) in params.items():
#            if param in extra_params:
#                if param == 'agc':
#                    klatt_defs.gain_control = int(val)
#            else:
#                # kl_np corresponds to np in C code. We don't use 'np' here
#                # because it conflicts with the usual use of 'np' for numpy.
#                kl_np = params_map[param]
#                if np.isscalar(val):
#                    klatt_defs.defval[kl_np] = int(val)
#                    self.update_constant_params_from_defval()
#                else:
#                    for (nf, v) in enumerate(val):
#                        klatt_defs.cv[kl_np] = klatt_defs.VARRIED
#                        klatt_defs.pdata[kl_np][nf] = int(v)
#                    num_varied += 1
#        klatt_defs.nvar = num_varied
#        return self.get_params()
#
#    def get_ms_per_frame(self):
#        ''' Get the milliseconds-per-frame value. '''
#        return klatt_defs.ms_frame
#
#    def synthesize(self):
#        ''' Synthesize waveform based on state of global variables in C code and return the waveform data and sample rate. '''
#        klatt_defs.actonrequest('s', 1) # Assume synthesis and batch mode.
#        # TODO: don't assume one channel, or disable multichannel synthesis in C code
#        pywav = np.zeros([klatt_defs.nsamtot], dtype=np.int16)
#        # TODO: figure out how to access klatt_defs.iwave without copying
#        for idx in range(klatt_defs.nsamtot):
#            pywav[idx] = klatt_defs.iwave[idx]
#        return (pywav, self.get_defval('sr'))
#
