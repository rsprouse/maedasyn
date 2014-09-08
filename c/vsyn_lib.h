#ifndef VSYN_LIB_H
#define VSYN_LIB_H

/*****
*	File :	vsyn_lib.h
*	Note :	Prototypes etc..
*****/

#define	DACscale 10000.0		/* scale-up signal to fit 16 bits integer */
#define gltDACscale 100.*DACscale	/* sacle up for glottal source */

float	glottal_area( char model, char mode, float Ap, short *t0 );
void	vowel_synthesis( FILE *sig_file );

short	vtt_ini( );
float	vtt_sim( );
void	vtt_term( void );

#endif
