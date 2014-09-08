#ifndef LAM_LIB_H
#define LAM_LIB_H

/********
*	File :	lam-lib.h
*	Note :	The definition of structures and prototypes.
********/

#include "vtconfig.h"
#define NP	 29	/* = np        */

typedef	struct{ short x, y;} int2D;
typedef	struct{ float x, y;} float2D;

extern float2D	ivt[NP];		/* VT inside contours             */
extern float2D	evt[NP];		/* VT exterior contours           */


void	read_model_spec( void );
void	convert_scale( void );
void	semi_polar( void );
void	lam( float *para);
void	sagittal_to_area( short *ns, area_function *af );
void	appro_area_function (short ns1, area_function *A1, short ns2, area_function *A2);
void	print_lam ( void );
void    print_af (short ns, area_function *af);
void	plot_semi_polar( short nvp );
void	plot_lam( short nvp );
void	deplot_lam( void );

#endif
