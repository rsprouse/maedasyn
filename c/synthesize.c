//
//  synthesize.c
//  maeda
//
//  Created by Keith Johnson on 1/22/14.
//  Copyright (c) 2014 Keith Johnson. All rights reserved.
//

#include	"always.h"
#include	"lam_lib.h"
#include	"vtconfig.h"
#include	"vsyn_lib.h"

/* parameter matrix - a sequence of frames  */
#define NPAR 10     /* number of model parameters per frame */
#define AMloc 3     /* starting location of the 7 Maeda model parameters */
#define AMnum 7     /* how many Maeda model parameters are there */
#define TIME 0  /* index of the time value of the frame */
#define F0_LOC 1 /* index of the f0 value */
#define AP 2 /* target amplitude of glottal opening */
#define FRAME_DUR 0.005  /* duration (seconds) of a frame */

/* NOTE BY RLS
  The values of NPAR and AMnum are incorrect as given. They should be 11 and 8,
  i.e. they are really 'last index', not 'number of'.
  This probably hasn't mattered since the last parameter (nasal coupling) is
  is always 0 in this code.
*/

/*   The columns in <par>:
 0 time
 1 f0
 2 max glottal opening
 3 Jaw position             <- Artic model params start here
 4 Tongue dorsum position
 5 Tongue dorsum shape
 6 Tongue apex position
 7 Lip height (aperture)
 8 Lip protrusion
 9 Larynx height
 10 Nasal coupling (cm2)  *** not used? ***
  */

/* Specific to time domain calculations */

float	simfrq = 30000.0f;	/* simulation frequency in Hz	*/
float	smpfrq = 10000.0f;	/* sampling freq. in Hz		*/

/*( Used in the time and frequency domain calculations )*/

float	Psub = 8.f;		/* subglottal air pressure in cmH2O	*/
float	Ag = 0.0f;		/* glottis area in cm**2		*/
float	xg = 0.3f;		/* glottis thickness in cm (=0.3 cm)	*/
float	lg = 1.2f;		/* fold length in cm			*/
float	Kc = 1.42f;		/* Ishizaka's for turbulent flow	*/
/* =.875, Van der Berg's constant	*/

short	nbu = 8;		/* # of sections in the bucal tube	*/
short	nph = 9;		/* # of sections in the phryngeal tube	*/
short	nvt;			/* = nbu + nph, # of VT sections	*/
area_function	*afvt;		/* AF from the glottis to the lips	*/

float	anc = 0.f;		/* nasal coupling area in cm2		*/
short	nna = 13;			/* # of sections in the nasal tract	*/
/* AF from nostrils to coupling point	*/
area_function afnt_array[13] = {
    {0.25, 1.},{6.0, 1.}, {6.0, 1.}, {6.0, 1.}, {6.0, 1.}, {6.0, 1.}, {6.0, 1.}, {6.0, 1.},
    {6.0, 1.}, {6.0, 1.}, {4.0,1.}, {2.0,1.}, {0.1, 1.}};
area_function	*afnt = afnt_array;


short	nss;			/* total # of sections in VT system	*/

/************************( simulation options )**************************/

short	nasal_tract  = OFF;		/* or ON			*/
short	wall         = YIELDING;	/* or RIGID			*/
short	rad_boundary = BESSEL_FUNCTION;	/* RL_CIRCUIT SHORT_CIRCUIT, or BESSEL_FUN	*/
short	glt_boundary = CLOSE;		/* or OPEN			*/

short	source_loc = 0;		/* source location in VT section number	*/
short	source_typ = FLOW;		/* or PRESSURE			*/

/* specific to time_domain calculations */

short	vocal_tract  = TIME_VARYING;	/* STATIONARY or TIME_VARYING		*/
short	dynamic_term = OFF;		/* or OFF			*/

/************( an extra heat loss factor for the nasal tract )***********/

float	extra_loss_factor = 50.f;

/************************( physical constants )**************************/

float	ro    = 1.14e-3f;	/* air density, gm/cm**3		*/
float	c     = 3.5e+4f;		/* sound velocity, cm/s			*/
float	eta   = 1.4f;		/* adiabatic constant			*/
float	cp    = 0.24f;		/* specific heat, cal/gm.degree		*/
float	lamda = 5.5e-5f;		/* heat conduction, cal/cm.sec.degree	*/
float	mu    = 1.86e-4f;	/* viscosity coef, dyne.sec/cm**2	*/
float	wall_resi = 1600.f;	/* wall resistance, gm/s/cm2		*/
float	wall_mass = 1.5f;	/* wall mass per area, gm/cm2		*/
float	wall_comp = 3.0e+5f;	/* wall compliance			*/
float	H2O_bar= 980.39f;



long update_VT(float **par, long buf_count, int time_steps){
    
    short	ns0 = 29;
	static	area_function	*af0;
    int target_time, i;
    float AMpar[7];  /* articulatory model parameters */
    
    /* Initialization */
    if (buf_count==0L) {
        
        af0 = (area_function *) calloc( ns0, sizeof(area_function) );
        nss = nbu + nph;
        afvt  = (area_function *) calloc( nss, sizeof(area_function) );
        convert_scale();
        semi_polar();
    
        for (i=0;i<AMnum;i++) AMpar[i] = par[0][i+AMloc];  /* first vocal tract shape */

        lam( AMpar );				/* compute VT sagittal section */
        printf("time = %0.3f\n",(buf_count/smpfrq));
        //print_lam();
        sagittal_to_area( &ns0, af0 );		/* compute area function from sagittal section */
        appro_area_function( ns0, af0, nss, afvt);  /* make tube lengths equal */
        print_af(nss,afvt);
        vtt_ini();
        return (smpfrq*FRAME_DUR);  /* next update is due */
    }
    target_time = buf_count/smpfrq * 1000;  /* get parameters for frame at target_time */
    do time_steps--; while (par[time_steps][TIME] > target_time);
    for (i=0;i<AMnum;i++) AMpar[i] = par[time_steps][i+AMloc];
    
    lam( AMpar );				/* compute VT sagittal section */
    printf("time = %0.3f\n",(buf_count/smpfrq));
    //print_lam();
    sagittal_to_area( &ns0, af0 );		/* compute area function from sagittal section */
    appro_area_function( ns0, af0, nss, afvt);  /* make tube lengths equal */
    print_af(nss,afvt);
   
    if( nasal_tract == ON ){  /* add nose */
        anc = (float) min( anc, afvt[nph].A );
        afvt[nph].A -= anc;
    }
    return (buf_count + smpfrq*FRAME_DUR);  /* 0.005 = 5 ms */
    
}

short update_pitch(float **par,long buf_count, int time_steps, float *Ap) {
    int target_time;
    
    target_time = buf_count/smpfrq * 1000;  /* get parameters for frame at target_time */
    do time_steps--; while (par[time_steps][TIME] > target_time);
    
    *Ap = par[time_steps][AP];
    return (short) ( 0.5+ smpfrq/par[time_steps][F0_LOC]);

}


/* synth_frame 
 input: 
 par: array of parameters
 buffer:  a pointer to a buffer in which one frame (5ms) of speech samples will be saved
            size of buffer must be at least FRAME_DUR*smpfrq
 mode: 1 = initialize, 2 = normal, >2 = fade-out mode for "mode" samples
 */
void synth_frame(float *params, short *buffer, short mode) {
    short i;
    short	ns0 = 29;
	static	area_function	*af0;
    static short period,t0;
    float AMpar[7];
    float Ap = 0.2;
    short nsamp = FRAME_DUR*smpfrq;
    
    for (i=0;i<AMnum;i++) AMpar[i]=params[i+AMloc];
    
    if (mode==1) {  // initialize
        
        af0 = (area_function *) calloc( ns0, sizeof(area_function) );
        nss = nbu + nph;
        afvt  = (area_function *) calloc( nss, sizeof(area_function) );
        convert_scale();
        semi_polar();
        
        lam(AMpar);				/* compute VT sagittal section */
        sagittal_to_area( &ns0, af0 );		/* compute area function from sagittal section */
        appro_area_function( ns0, af0, nss, afvt);  /* make tube lengths equal */
        vtt_ini();
        
        Ap = params[AP];
        t0 = period = (short)(0.5 + smpfrq/params[F0_LOC]);
        
    }

    if (mode==2) {  // normal
        lam(AMpar);				/* compute VT sagittal section */
        sagittal_to_area( &ns0, af0 );		/* compute area function from sagittal section */
        appro_area_function( ns0, af0, nss, afvt);  /* make tube lengths equal */
        for (i=0;i<nsamp;i++) {
            Ag = glottal_area( 'F', 'o', Ap, &t0 );  /* voice source */
            buffer[i] = (short) (DACscale * vtt_sim());  /* synthesize next sample */
            period--;
            if (period <= 0) {
                t0 = period = (short)(0.5 + smpfrq/params[F0_LOC]);
                Ap = params[AP];
            }
        }
    }
    if (mode > 2) {  //fade-out mode - buffer must be long enough to accommodate mode ms of samples
        t0 = period = (short)(0.5 + smpfrq/params[F0_LOC]);
        Ap = params[AP];
        for (i=0;i<mode;i++) {
            Ag = glottal_area( 'F', 't', Ap, &t0 );  /* voice source  'transition' */
            buffer[i] = (short) (DACscale * vtt_sim());  /* synthesize next sample */

        }
        vtt_term();

    }
    
}

/*  ----------------------------synthesize -----------------------------
 Inputs:
 
 par: a two dimensional array of control parameters - column 1 is time value for each frame (in samples)
 sig_buf: a pointer - memory allocation will take place inside the routine and sound wave will be saved here
 
 Output:
 returns number of samples in sig_buf
 */


long synthesize(float **par, short **sig_buf, int time_steps) {
    
    long buf_length, nextVTupdate, nextPitchUpdate, buf_count = 0L;
    short t0, temp;  /* number of samples in a pitch period */
    float	Ap = 0.2;
    
    buf_length = (par[time_steps-1][TIME]/1000)*smpfrq;  /* duration in sec */
    
    if ((*sig_buf = (short *) calloc( buf_length + smpfrq*0.06, sizeof(short) ))==NULL) {
        fprintf(stderr,"%s/n","error allocating memory for wave");
    }
    
    nextVTupdate = update_VT(par,buf_count,time_steps);   /* set initial VT area function */
    
    
    t0 = update_pitch(par, buf_count, time_steps, &Ap);  /* set initial pitch period */
    nextPitchUpdate = t0;
    
    while (buf_count < buf_length) {
        
        Ag = glottal_area( 'F', 'o', Ap, &t0 );  /* voice source */
        temp = (short) (DACscale * vtt_sim());  /* synthesize next sample */
        //printf("%ld\t%0.3f\t%hd\n",buf_count,Ag,temp);
        (*sig_buf)[buf_count++] = temp; /* add the sample to the sound buffer */
        
        if (buf_count >= nextVTupdate) {
            nextVTupdate = update_VT(par,buf_count,time_steps);   /* move mouth every 5 ms */
        }
        if (buf_count >= nextPitchUpdate) {
            t0 = update_pitch(par, buf_count, time_steps, &Ap);  /* change pitch every cycle */
            nextPitchUpdate += t0;
        }
    }
    t0=update_pitch(par,buf_count,time_steps, &Ap);
    while (buf_count < buf_length + smpfrq*0.06)  {  /* 'transition' glottal vibration - 60 ms */
        Ag = glottal_area( 'F', 't', Ap, &t0 );
        temp = (short) (DACscale * vtt_sim());  /* synthesize next sample */
        //printf("%ld\t%0.3f\t%hd\n",buf_count,Ag,temp);
        (*sig_buf)[buf_count++] = temp;
	}
    
    return buf_count;
}


/* some Maeda model specifications for vowels
0 Jaw position  1 Tongue dorsum position    2 Tongue dorsum shape   3 Tongue apex position
4 Lip height    5 Lip protrusion            6 Larynx height         7 Nasal coupling (cm2) */
float aa[7]=  {   -1.5,  2.0, 0.0, -0.5,  0.5, -0.5, 0.0};    /* aa */
float uw[7]=  {   0.5, -1.0, 1.0, -2.0, -0.5,  1.0, 0.0 };    /* uw */
float iy[7]=  {   0.5, -2.0, 1.0, -2.0,  1.0, -1.0, 0.0 };    /* iy */
float ey[7]=  {   0.0, -1.0, 1.0, -2.0,  1.0, -1.0, 0.0 };    /* ey */
float eh[7]=  {   -1.0,  0.0, 1.0, -2.0,  1.0, -0.5, 0.0 };   /* eh */
float ah[7]=  {   -1.5,  0.5, 0.0, -0.5,  0.5, -0.5, 0.0 };   /* ah */
float ao[7]=  {   -0.4,  3.0, 1.5,  0.0, -0.3,  0.0, 0.0 };   /* ao */
float ow[7]=  {   -.7,  3.0, 1.5,  0.0, -0.6,  0.0, 0.0 };    /* ow */
float iw[7]=  {   0.5, -1.0, 1.0, -2.0, -0.5,  1.0, 0.0 };    /* iw */
float ew[7]=  {   0.0, -0.2, 1.0, -1.5, -0.25, 0.5, 0.0 };    /* ew */
float oe[7]=  {   -1.0, -0.5, 0.5, -2.0,  0.2, -0.5, 0.0 };   /* oe */

int main() {
    
    short bufsize = (short)(FRAME_DUR*smpfrq);
    short buffer[bufsize];
    short final_buffer[bufsize*10];
    float parameters[NPAR],target1[NPAR],target2[NPAR];
    short num_written;
    int i,j,n;
    float d[NPAR];
    FILE *soundfile;
    
    
    if( (soundfile = fopen("../resources/mtest2.raw", "w+b")) == NULL ) {
        fprintf(stderr,"Can't open file.");
        exit(-1);
	}
    
    target1[F0_LOC] = 130;
    target2[F0_LOC] = 100;
    target1[AP]=0.2;
    target2[AP]=0.2;
    for (j=0; j<AMnum; j++) target1[AMloc+j] = uw[j];
    for (j=0; j<AMnum; j++) target2[AMloc+j] = iy[j];

    n= 80-20;
    for (j=1; j<NPAR; j++) {
        d[j] = (target2[j] - target1[j])/n;  // delta values for interpolation
        parameters[j] = target1[j];  // starting params for initialization
    }
    
    synth_frame(parameters, buffer, 1);  // mode 1 = initialize the synthesizer

    
    /* construct a par array for testing */
    for (i=0; i<100; i++) {
        parameters[TIME] = (i*FRAME_DUR)*1000;
        if (i>=80) {for (j=0; j<NPAR; j++) parameters[j] = target2[j];}
        if (i<=20) {for (j=0; j<NPAR; j++) parameters[j] = target1[j];}
        if (i>20 && i<80) {
            for (j=0; j<NPAR; j++) parameters[j] = target1[j] + d[j]*(i-20.0);
        }
        synth_frame(parameters, buffer, 2);
        
        /* at this point in the calling code you can refer to evt[i].x and evt[i].y to 
            draw the "passive" surface of the vocal tract, and ivt[i].x and ivt[i].y to
            draw the "active" surface.
         
         
        for (j=0;j<NP;j++) {
            printf("%d\t%3.4f\t%3.4f\n",j,evt[j].x,evt[j].y);
            
        }
        */
         
        if ((num_written = fwrite( buffer, sizeof(short), bufsize, soundfile )) != bufsize) {
            printf("%s\n","write to file failed");
        }
    }
    synth_frame(parameters,final_buffer,bufsize*10);
    if ((num_written = fwrite( final_buffer, sizeof(short), bufsize*10, soundfile )) != bufsize*10) {
        printf("%s\n","write to file failed");
    }
    
    fclose(soundfile);

    
}
/*
int main() {
    
    int time_steps=100;  // 500 ms at 5 ms per step
    short *sound_buffer = NULL;
    long buf_length,num_written;
    int i,j,n;
    float d[AMnum];
    FILE *soundfile;
    
    
	float **par;
	par = malloc(time_steps * sizeof(float *));
	if(par == NULL) {
		fprintf(stderr, "out of memory\n");
		exit(1);
    }
	for(i = 0; i < time_steps; i++)   {
		par[i] = malloc(NPAR * sizeof(float));
		if(par[i] == NULL)
        {
			fprintf(stderr, "out of memory\n");
			exit(1);
        }
    }
    
    // construct a par array for testing
    for (i=0; i<100; i++) {
        par[i][TIME] = (i*FRAME_DUR)*1000;
        par[i][F0_LOC] = 100;
        par[i][2] = 0.2;
    }
    for (i=0; i<20; i++) {for (j=0; j<AMnum; j++) par[i][AMloc+j] = iy[j];}
    for (i=80;i<100; i++) {for (j=0; j<AMnum; j++) par[i][AMloc+j] = uw[j];}
    
    for (i=20; i<80; i++) {
        n= 80-20;
        for (j=0; j<AMnum; j++) d[j] = (par[80][j+AMloc] - par[19][j+AMloc])/n;
        // interpolate
        for (j=0; j<AMnum; j++) par[i][j+AMloc] = par[i-1][j+AMloc] + d[j];
    }
    for (i=0; i<100; i++) {
        for (j=0; j<NPAR; j++) {printf("%0.3f\t",par[i][j]); }
        printf("\n");
    }
    
    buf_length = synthesize(par, &sound_buffer, time_steps);
    
    
	vtt_term();
    
    if( (soundfile = fopen("/Users/transferredkj/Google Drive/maeda/mtest2.raw", "w+b")) == NULL ) {
        fprintf(stderr,"Can't open temporary SIG file.");
        exit(-1);
	}
    
    if ((num_written = fwrite( sound_buffer, sizeof(short), buf_length, soundfile )) != buf_length) {
        printf("%s\n","write to file failed");
    }
    fclose(soundfile);
    
    
}
*/
