/***************************************************************************
*	File :	vsyn_lib						   *
*	Note :	This file contains functions for synthesis of vowels.	   *
***************************************************************************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include    "vtconfig.h"
#include	"vsyn_lib.h"

/****************( externally defined global variables )*******************/

extern	float	smpfrq;
extern	float	Ag;

/****************************( functions )*********************************/
/*****
*	Function : glottal_area
*	Note :	To calculate glottal area (Ag cm2) at sample point n during
*		a fundamental period or a transition quotient (t0 samples).
*		In the oscilation mode (mode = 1), Ag is specified either
*		Fant's model or Maeda's model.  In the transition mode
*		(mode = 2), the Ag variation from the initial to the target
*		area is specified by a raised cosine curve.  The calculated
*		area is returned by value to the calling program.
*****/

float	glottal_area (
	char	model,		/* 'F' for Fant, 'M' for Maeda model	  */
	char	mode,		/* 'o' for oscilating, 't' for transition */
	float	Ap,		/* peak glottal area (cm2) with mode=1,	  */
				/* target area (cm2) with mode=2	  */
	short	*t0 )		/* fundamental period or transion quotient*/
				/* in samples. The non-zero value signals */
				/* a new glottal cycle or transion.	  */
				/* t0 will be zeroed			  */
{
	float	oqF, cqF;	/* opening & closing quotients for Fant's */
				/* model.				  */
	float	wp, oqM, cqM;	/* time warping coefficients for Maeda's  */
				/* model.				  */
	static	short	n;
	static	float	amp;
	static	short	t1, t2, t3;
	static	float	a, b, A, A0;
	float	t;

/*** oscilation mode ***/
	if(mode == 'o' ) {
        if( *t0 > 0 ) {		/* set a new glottal cycle */
            A = (float)0.5*Ap;
            if( model == 'F'){
                oqF = 0.36f;
                cqF = 0.26f;
                t1 = (short)(oqF * (*t0));
                t2 = (short)(cqF * (*t0));
                t3 = t1 + t2;
                a  = 3.141593f/t1;
                b  = (float)(1./(1. - cos(a*t2)));
            }
            if( model == 'M') {
                wp = 2.f;
                oqM = 0.5f;
                cqM = 0.2f;
                t1 = (short)(oqM * (*t0));
                t2 = (short)(cqM * (*t0));
                t3 = t1 + t2;
                a  = (float)(3.141593f/t1);
                b  = (float)( (t1 - t2)/pow( t2, wp ));
            }
            *t0 = 0;
            n   = 0;
        }

        if( n < t1 ) amp = (float)(A*(1.0 - cos(a*n)));	/* opening */
        if( n >= t1 && n < t3 )	{		/* closing */
            t = (float)(n - t1);
            if( model == 'F' ) amp = (float)(Ap*(1. - b + b*cos(a*t)));
            if( model == 'M' ) amp = (float)(A*(1. + cos(a*(t+pow(t,wp)))));
        }
        if( n >= t3 ) amp = 0.0;			/* closed */
        n++;
        return( amp );
	}
/*** transion mode ***/
	if( mode == 't' ) {
        if( *t0 > 0 ) {
            t1 = *t0;
            A0 = amp;
            A  = (float)(0.5*( Ap - A0 ));
            a  = (float)(3.141593/t1);
            n  = 0;
            *t0 = 0;
        }

        if( n < t1 ) amp = (float)(A0 + A*(1. - cos(a*n)));
        else         amp = Ap;
        n++;
        return( amp );
	}
	return 0;
}

/*****
*	Function : vowel_synthesis
*	Note :	Synthesis of a stationary vowels with varying F0.
*****/

void	vowel_synthesis ( FILE *sig_file ) {
	float	Ap = 0.2f;
	long    buf_count;
	short	*sig_buf;

    float	f[3] = { 120., 130., 100. };	/* F0 pattern, Hz */
    float	d[3] = { 0., 100., 200. };	/* duration, ms */

	float	p_trgt[3], d_trgt[3];
    long    num_samples;   /* how long is this sound in samples? */
    long    num_written;                                                  /* we add room for three cycles of "transition" at 50 Hz*/
	float	a, p, n;
	short	np, t0;
	short	i, j;
    short   temp;
    float   duration = 0;
    
    /* define F0 contour */
    for(i=0; i<3; i++){
        p_trgt[i] = smpfrq/f[i];		/* target pitch (samples) */
        d_trgt[i] = (float)(smpfrq*d[i]/1000.);	/* duration in samples    */
        duration += d[i];
	}
    num_samples = (smpfrq * duration/1000.) + (3*(smpfrq/50));
    /* allocate buffer */
	if ((sig_buf = (short *) calloc( num_samples, sizeof(short) ))==NULL) {
        fprintf(stderr,"%s/n","error allocating memory for wave");
    }
    /* initialize simulator */
	Ag = 0;
	vtt_ini(afvt);		/* discard returned value, constant delay */

    /* pitch loop */
	n = 0;
	buf_count = 0;
	for(j=0; j<2; j++)  {  /* looping here over a sequence of target values */
        a = (p_trgt[j+1] - p_trgt[j])/d_trgt[j+1];  /* slope of pitch contour for interpolation */
        n = n - d_trgt[j];   /* time location in pitch contour */
        p = 0;
        do  {              /* we're doing a pitch-synchronous update of f0 target */
            n = n + p;    /* sample count in the output buffer */
            p = a*n + p_trgt[j];   /* updated pitch period in samples */
            t0 = np = (short)(p + 0.5);   /* pitch period as an integer */
            //printf("%d\t%ld\t%0.2f\t%0.2f\n",np,buf_count,n,d_trgt[2]);
            /* time vary vocal tract parameters here */

            for(i=0; i<np; i++){		/* single period */
                Ag = glottal_area( 'F', 'o', Ap, &t0 );

                temp = (short) (DACscale * vtt_sim());  /* synthesize next sample */
                printf("%ld\t%0.3f\t%hd\n",buf_count,Ag,temp);
                sig_buf[buf_count++] = temp;
            }
        } while( n < d_trgt[j+1] );  /* keep synthesizing until the number of samples is reached */
	}
	t0 = np;			/* open glottis to decay signal */
	for(j=0; j<3; j++)  {  /* add three "transition" glottal pulses */
        for(i=0; i<np; i++)  {
           Ag = glottal_area( 'F', 't', Ap, &t0 );
            temp = (short) (DACscale * vtt_sim(afvt));  /* synthesize next sample */
            printf("%ld\t%0.3f\t%hd\n",buf_count,Ag,temp);
            sig_buf[buf_count++] = temp;
        }
	}
    if ((num_written = fwrite( sig_buf, sizeof(short), num_samples, sig_file )) != num_samples) {
        printf("%s\n","write to file failed");
    }

	vtt_term();
}
