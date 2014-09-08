#ifndef VTCONFIG_H
#define VTCONFIG_H

/******
*	File :	vtconfig.h
*	Note :	Acoustic and geometrical configurations for the vocal
*		tract calculations.  It should be declared in a "main".
******/

#define	OFF		0x00
#define	ON		0x01
#define	RIGID		0
#define	YIELDING	1
#define	FLOW		0
#define	PRESSURE	1
#define CLOSE		0
#define OPEN		1
#define	SHORT_CIRCUIT	0
#define RL_CIRCUIT	1
#define	BESSEL_FUNCTION	2

#define	TIME_VARYING	1
#define STATIONARY	0

/*************************( typedefs )***********************************/
typedef struct{ float A, x; } area_function;


/******************( Global variables defined in "main" )****************/


	/* Specific to time domain calculations */

extern float	simfrq;	/* simulation frequency in Hz	*/
extern float	smpfrq;	/* sampling freq. in Hz		*/

	/*( Used in the time and frequency domain calculations )*/

extern float	Psub;		/* subglottal air pressure in cmH2O	*/
extern float	Ag;		/* glottis area in cm**2		*/
extern float	xg;		/* glottis thickness in cm (=0.3 cm)	*/
extern float	lg;		/* fold length in cm			*/
extern float	Kc;		/* Ishizaka's for turbulent flow	*/
				/* =.875, Van der Berg's constant	*/

extern short	nbu;		/* # of sections in the bucal tube	*/
extern short	nph;		/* # of sections in the phryngeal tube	*/
extern short	nvt;			/* = nbu + nph, # of VT sections	*/
extern area_function	*afvt;		/* AF from the glottis to the lips	*/

extern float	anc;		/* nasal coupling area in cm2		*/
extern short	nna;			/* # of sections in the nasal tract	*/
/* AF from nostrils to coupling point	*/
extern area_function	*afnt;

extern short	nss;			/* total # of sections in VT system	*/

/************************( simulation options )**************************/

extern short	nasal_tract;		/* or ON			*/
extern short	wall;	/* or RIGID			*/
extern short	rad_boundary;	/* SHORT_CIRCUIT, or BESSEL_FUN	*/
extern short	glt_boundary;		/* or OPEN			*/

extern short	source_loc;		/* source location in VT section number	*/
extern short	source_typ;		/* or PRESSURE			*/

	/* specific to time_domain calculations */

extern short	vocal_tract;	/* or TIME_VARYING		*/
extern short	dynamic_term;		/* or OFF			*/

/************( an extra heat loss factor for the nasal tract )***********/

extern float	extra_loss_factor;

/************************( physical constants )**************************/

extern float	ro;	/* air density, gm/cm**3		*/
extern float	c;		/* sound velocity, cm/s			*/
extern float	eta;		/* adiabatic constant			*/
extern float	cp;		/* specific heat, cal/gm.degree		*/
extern float	lamda;		/* heat conduction, cal/cm.sec.degree	*/
extern float	mu;	/* viscosity coef, dyne.sec/cm**2	*/
extern float	wall_resi;	/* wall resistance, gm/s/cm2		*/
extern float	wall_mass;	/* wall mass per area, gm/cm2		*/
extern float	wall_comp;	/* wall compliance			*/
extern float	H2O_bar;

#endif

