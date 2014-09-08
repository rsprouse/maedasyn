/***************************************************************************
*                                                                          *
*	File : lam_lib.c						   *
*	Note : functions for linear articulatory model (LAM).              *
*                                                                          *
***************************************************************************/

#include	"always.h"
#include	"vtconfig.h"
#include	"lam_lib.h"

/************************( global )****************************/

/* Number of articulatory parameters */

#define	JAW	1	/* number of jaw parameter */
#define	LIP	2	/* number of intrinsic lip parameters (HT and PR) */
#define	TNG	3	/* number of intrinsic tongue parameters */
#define	LRX	1	/* number of intrinsic larynx parameter (height) */

/* Number of variables (for arrays) */

#define	M4	 31	/* =m1 + m2 + m3, total number of semi-polar grids */
#define	NVRS_LIP  4	/* = nvrs_lip  */
#define	NVRS_TNG 26	/* = nvrs_tng  */
#define	NVRS_LRX  5	/* = nvrs_lrx  */
#define NVRS_WAL 25	/* = nvrs_wal  */

const	float	pi = 3.14159265f;

/***************( global variables specified by data file )*************/

/* Semi-polar specification */
short	m1 = 14;   /* # of grids in the 3 regions          */
short m2 = 11;
short m3 = 6;
float dl = 0.5;  /* grid spacing in cm and degree        */
float omega = -11.25;
float theta = 11.25;
short	ix0 = 3000;
short iy0 = 1850;	  /* the origin on 4095*4095 TEK space    */

float TEKvt = 188.679245;
float TEKlip = 0.00;	  /* from cm to TEK unit map (points/cm)  */
/* alpha-beta sagittal-area coefficients */

static float alph[M4] = {1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.8,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.7,1.8,1.8,1.9,2.0,2.6};
static float beta[M4] = {1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.2,1.3,1.4,1.4,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5,1.5};

/* Lip-tube specification */
short	nvrs_lip = 4;		/* number of variables            */
short	jaw_lip = 1;		/* number of jaw variables(0 or 1)*/
/* labals for factors (parameters)*/
static	char	flab_lip[JAW+LIP][3] = { "JW", "HT", "P1"};
/* mean values in TEK unit        */
static	float	u_lip[NVRS_LIP] = { 104.271675, 122.812141, 135.938339, 460.440857};
/* standard deviations in TEK unit*/
static	float	s_lip[NVRS_LIP] = { 27.674635, 33.068081, 99.392258, 213.996170};
/* mean maxillary incisor position*/
float inci_x = 2212.354492;
float inci_y = 1999.574219;
/* factor patterns (loadings) */
static	float	A_lip[NVRS_LIP][JAW+LIP] =
    {{1.000000,   0.000000,   0.000000},
    {0.178244,  -0.395733,   0.888897},
    {-0.154638,   0.987971,   0.000000},
    {-0.217332,   0.825187,  -0.303429}};

/* tongue-contour specification */
short	nvrs_tng =26;		/* number of variables            */
short	jaw_tng = 1;		/* number of jaw variables(0 or 1)*/
short	iniva_tng = 7;	/* coordinate number of initial point     */
short	lstva_tng = 31;	/* coorninate number of last point        */
/* labals for factors (parameters)*/
static	char	flab_tng[JAW+TNG][3]= {"JW", "P1", "P2", "P3"};
/* mean values in TEK unit        */
static	float	u_tng[NVRS_TNG] = { 104.271675,     443.988434,     450.481689,   399.942200,     348.603088,    351.181122,     365.404633,     370.290955,     356.202301,     341.890167,    332.117523,     326.826599,     326.512512,     331.631989,     343.175323,    361.265900,     385.231201,     411.826599,     435.691711,     455.040466,    462.736023,     453.025055,     432.250488,     407.358368,     384.551056,    363.836212};
/* standard deviations in TEK unit*/
static	float	s_tng[NVRS_TNG] = { 27.674635,      29.947931,      44.694466,      99.310226,      96.871323,    84.140404,      78.357513,      73.387718,      72.926758,      71.453232,    69.288765,      66.615509,      63.603722,      59.964859,      56.695446,    56.415058,      62.016468,      73.235176,      84.008438,      91.488312,    94.124176,      95.246323,      93.516365,      93.000343,     100.934669,    106.512482 };
/* factor patterns (loadings) */
static	float	A_tng[NVRS_TNG][JAW+TNG] = {
    {1.000000,   0.000000,   0.000000,   0.000000},
    {-0.464047,   0.098776,  -0.251690,   0.228351},
    {-0.328015,   0.337579,  -0.283667,   0.568234},
    {-0.213039,   0.485565,  -0.283533,   0.653696},
    {-0.302565,   0.705432,  -0.379044,   0.392917},
    {-0.327806,   0.786897,  -0.388116,   0.245703},
    {-0.325065,   0.852409,  -0.285125,   0.176843},
    {-0.325739,   0.904725,  -0.142602,   0.138558},
    {-0.313741,   0.926339,   0.021042,   0.122976},
    {-0.288138,   0.924019,  0.131949,   0.116762},
    {-0.249008,   0.909585,   0.250320,   0.112433},
    {-0.196936,   0.882236,   0.369083,   0.112396},
    {-0.128884,   0.830243,   0.499894,   0.115700},
    {-0.040825,   0.730520,   0.651662,   0.112048},
    {0.073420,   0.543080,  0.807947,   0.126204},
    {0.202726,   0.230555,   0.919065,   0.163735},
    {0.298853,  -0.162541,   0.899074,   0.213884},
    {0.332785,  -0.491647,   0.748869,   0.243163},
    {0.349955,  -0.681313,   0.567615,   0.245295},
    {0.377277,  -0.771200,   0.410502,  0.249425},
    {0.422713,  -0.804874,   0.270513,   0.274015},
    {0.474635,  -0.797704,   0.129324,   0.314454},
    {0.526087,  -0.746938,  -0.026201,   0.366149},
    {0.549466,  -0.643572,  -0.190005,   0.422848},
    {0.494200,  -0.504012,  -0.350434,   0.488056},
    {0.448797,  -0.417352,  -0.445410,   0.500909}};
/* Larynx-tube specification */
	short	nvrs_lrx = 5;		/* number of variables            */
	short	jaw_lrx = 1;		/* number of jaw variables(0 or 1)*/
	short	iniva_lrx = 7;	/* coordinate number of initial point     */
	short	lstva_lrx = 6;	/* coordinate number of last point        */
/* labals for factors (parameters)*/
static	char	flab_lrx[JAW+LRX][3] = {"JW","Y1"};
/* mean values in TEK unit        */
static	float	u_lrx[NVRS_LRX] ={104.271675, 143.138733, -948.229309, 404.678223, -962.936401};
/* standard deviations in TEK unit*/
static	float	s_lrx[NVRS_LRX] = {27.674635, 41.593315, 65.562340, 44.372742, 66.147499};
/* factor patterns (loadings) */
static	float	A_lrx[NVRS_LRX][JAW+LRX]  =
    {{ 1.00000,   0.00000},
    {-0.208338,   0.262446},
    {0.127814,   0.991798},
    {-0.131840,   0.300784},
    {0.097688,   0.934267}};

/* Wall-contour specification */
	short	nvrs_wal = 25;		/* number of variables            */
	short	jaw_wal = 0;		/* number of jaw variables(0 or 1)*/
	short	iniva_wal = 7;	/* coordinate number of initial point     */
	short	lstva_wal = 31;	/* coorninate number of last point        */
/* mean values in TEK unit        */
static	float	u_wal[NVRS_WAL] = {550.196533,     604.878601,    674.127197,     678.776489,     665.905579,    653.312134,     643.223511,    633.836243,     636.994202,     668.834290,    703.098267,     657.815002,    649.919067,     565.194580,     529.824646,    573.250488,     603.023132,    621.433533,     643.055847,     650.136780,    630.809265,     589.867065,    556.134888,     541.551086,     525.210022};

/*************************( global variables )***************************/

/* Semi-polar coordinate and VT contours */
float	vp_map = 1.0f;			/* viewport map coef. (cm/point)  */
const	float	size_correction = 1.10f; /* 10% size increase              */
const	float	vp_width_cm = 10.0f;	/* viewport height in cm          */
static	float2D	igd[M4], egd[M4];	/* Semi-polar coordinate grids    */
static	float2D	vtos[M4];		/* vector to semi-polar map       */
short	np;			/* number of points               */
float2D	ivt[NP];		/* VT inside contours             */
float2D	evt[NP];		/* VT exterior contours           */
const	float	inci_lip = 0.8f;	/* (cm) dist. btwn. incisor and upper lip */
float	inci_lip_vp;		/* inci_lip in viewport unit      */
float	lip_w, lip_h;		/* lip-tube width and height      */

/***************************( functions )**********************************/


/*****
*	Function : semi_polar
*		The semi-ploar is specified by x-y values of the two ends
*		of each grid line.
*
*****/

void	semi_polar ( void )
{
	float	r  = 5.0;		/* grid length (cm) */

	float	r_vp, dl_vp, ome, the, gam;
	float	dx_i, dy_i, dx_e, dy_e;
	float	p, q, s;
	short	i;

	r_vp  = r/vp_map;
	dl_vp = dl/vp_map;
	ome    = pi*omega/180.0f;
	the    = pi*theta/180.0f;

/* linear coordinate in the pharynx region */
	dx_i  = dl_vp*(float)cos(ome - pi/2.);
	dy_i  = dl_vp*(float)sin(ome - pi/2.);
	dx_e  = r_vp*(float)cos(ome);
	dy_e  = r_vp*(float)sin(ome);

	for(i=0; i<m1; i++)
	{  igd[i].x = dx_i*(m1 - (i + 1)) + ix0;
	   igd[i].y = dy_i*(m1 - (i + 1)) + iy0;
	   egd[i].x = dx_e + igd[i].x;
	   egd[i].y = dy_e + igd[i].y;
	}
/* polar coordinate in the velar region */
	for(i=m1; i<m1+m2; i++)
	{  gam = the*(i + 1 - m1) + ome;
	   igd[i].x = (float)ix0;
	   igd[i].y = (float)iy0;
	   egd[i].x = r_vp*(float)cos(gam) + ix0;
	   egd[i].y = r_vp*(float)sin(gam) + iy0;
	}
/* linear coordinate in the palato-dental region */
	dx_i  = dl_vp*(float)cos(gam + pi/2.0f);
	dy_i  = dl_vp*(float)sin(gam + pi/2.0f);
	dx_e  = r_vp*(float)cos(gam);
	dy_e  = r_vp*(float)sin(gam);

	for(i=m1+m2; i<m1+m2+m3; i++)
	{  igd[i].x = dx_i*(i + 1 - m1 - m2) + ix0;
	   igd[i].y = dy_i*(i + 1 - m1 - m2) + iy0;
	   egd[i].x = dx_e + igd[i].x;
	   egd[i].y = dy_e + igd[i].y;
	}
/* Mapping coefficients from vector to semi-polar */
	for(i=0; i<m1+m2+m3; i++)
	{  p = egd[i].x - igd[i].x;
	   q = egd[i].y - igd[i].y;
	   s = (float)sqrt( p*p + q*q );
	   vtos[i].x = p/s;
	   vtos[i].y = q/s;
	}
}

/*****
*	Function : lam
*	Note :	Compute vector representetion of the articulator positions
*		for given parameter values and project them on the semi-
*		polar coordinate to generate a VT profile.
*
*		Articulatory parameters are defined as follows :
*			para[0] : jaw
*			para[1] : tongue-body position
*			para[2] : tongue-body shape
*			para[3] : tongue-tip position
*			para[4] : lip height (aperture)
*			para[5] : lip protrusion
*			para[6] : larynx height
*
*		Note :	In order to avoid the crossover of contours, the lip
*			tube dimensions are blocked at zero and the tract
*			inside contour (tongue) is also blocked at the
*			exterior walls (10/08/92).
*
****/
void	lam ( float *pa )		/* a set of JAW+TNG+LIP+LRX */
{                                       /* articulatory parameter   */
	float	p[JAW+TNG];
	float	v_lip[NVRS_LIP], v_tng[NVRS_TNG], v_lrx[NVRS_LRX];
	float	v, x1, y1, x2, y2;
	short	i, j;

/*** copy parameter values and compute vectors ****/

	p[0] = pa[0];		    /* p[0] is always jaw parameter value */

/* tongue */
	for(i=1; i<=TNG; i++) p[i] = pa[i];
	for(i=0; i<nvrs_tng; i++)
	{  v = 0;
	   for(j=0; j<JAW+TNG; j++) v = v + A_tng[i][j]*p[j];
	   v_tng[i] = s_tng[i]*v + u_tng[i];
	}
/* lip */
	for(i=1; i<=LIP; i++) p[i] = pa[i+TNG];/* copy intrinsic lip pars.*/
	for(i=0; i<nvrs_lip; i++)
	{  v = 0;
	   for(j=0; j<JAW+LIP; j++) v = v + A_lip[i][j]*p[j];
	   v_lip[i] = s_lip[i]*v + u_lip[i];
	   if( v_lip[i] < 0. ) v_lip[i] = 0.;		/** block at zero **/
	}
/* larnx */
	for(i=1; i<=LRX; i++) p[i] = pa[i+TNG+LIP];
	for(i=0; i<nvrs_lrx; i++)
	{  v = 0;
	   for(j=0; j<JAW+LRX; j++) v = v + A_lrx[i][j]*p[j];
	   v_lrx[i] = s_lrx[i]*v + u_lrx[i];
	}

/*** Projection of vectors on the semi-polar coordinate ***/

/* larynx back edge */
	np=0;
	ivt[np].x = v_lrx[JAW]   + ix0;		/* front edge */
	ivt[np].y = v_lrx[JAW+1] + iy0;
	evt[np].x = v_lrx[JAW+2] + ix0;		/* rear edge */
	evt[np].y = v_lrx[JAW+3] + iy0;

/* larynx, pharynx and buccal */
	for(i=iniva_tng; i<lstva_tng; i++)  // this was i<=lstva_tng which left slice 26 at (0,0)
	{  j = i - iniva_tng;
	   /** block tongue contour at walls **/
	   v = (float)min( v_tng[j+JAW], u_wal[j] );
	   x1 = vtos[i].x * v + igd[i].x;		/* inside */
	   y1 = vtos[i].y * v + igd[i].y;
	   x2 = vtos[i].x * u_wal[j] + igd[i].x;	/* outside */
	   y2 = vtos[i].y * u_wal[j] + igd[i].y;

	   if(i == iniva_tng)			/* add an extra point */
	   {  ivt[++np].x = (ivt[0].x + x1)/2;
	      ivt[  np].y = (ivt[0].y + y1)/2;
	      evt[  np].x = (evt[0].x + x2)/2;
	      evt[  np].y = (evt[0].y + y2)/2;
	   }
	   ivt[++np].x = x1;
	   ivt[  np].y = y1;
	   evt[  np].x = x2;
	   evt[  np].y = y2;
	}
/* lips */
	evt[++np].x = inci_x + ix0;     /* pos. of inner edge of upper lip */
	evt[  np].y = inci_y + inci_lip_vp + iy0;
	ivt[  np].x = evt[np].x;
	ivt[  np].y = evt[np].y - v_lip[2];	/* lower lip */

    np++;
    
	evt[  np].x = evt[np-1].x - v_lip[1];
	evt[  np].y = evt[np-1].y;
	ivt[  np].x = evt[np].x;
	ivt[  np].y = ivt[np-1].y;

/* number of points */
	++np;

/*** Lip frontal shape (ellips) ***/
	lip_h = (float)(v_lip[2]/2.);
	lip_w = (float)(v_lip[3]/2.);
}

/*****
*	Function : amo
*	Note : returns the distance between two points.
*****/
float	amo (float2D p, float2D q)
{
	return( (float)sqrt((p.x-q.x)*(p.x-q.x) + (p.y-q.y)*(p.y-q.y)) );
}

/******
*	Function : sagittal_to_area
*	Note :	To calculate an area function of the vocal tract (VT)
*		from a quadrilateral representation of the profile.
*		The conversion employes an "alpha-beta" sagittal-to-
*		area relationships.  The alpha and beta values are
*		specified in relation with the semi-polar coordinate
*		grids.
*****/

void	sagittal_to_area (
	short	*ns,		/* number of sections */
	area_function	*af)	/* af.A = cross-sectional area (cm**2) */
				/* af.x = section length (cm) */
{
	float	p, q, r, s, t, a1, a2, s1, s2, x1, y1, d, w;
	float	c, cc;
	short	i, j;

/* vt_unit to cm conversion coef. with size_correction */
	c = size_correction*vp_map;
	cc = c*c;

/* from larynx to buccal */
	for(i=1; i<np-1; i++)
	{  p  = amo(ivt[i],   ivt[i-1]);
	   q  = amo(evt[i],   evt[i-1]);
	   r  = amo(ivt[i-1], evt[i-1]);
	   s  = amo(evt[i],   ivt[i]  );
	   t  = amo(evt[i],   ivt[i-1]);
	   a1 = (float)(0.5*(p + s + t));
	   a2 = (float)(0.5*(q + r + t));
	   s1 = (float)sqrt(a1*(a1 - p)*(a1 - s)*(a1 - t));
	   s2 = (float)sqrt(a2*(a2 - q)*(a2 - r)*(a2 - t));
	   x1 = ivt[i-1].x + evt[i-1].x - ivt[i].x - evt[i].x;
	   y1 = ivt[i-1].y + evt[i-1].y - ivt[i].y - evt[i].y;
	   d  = 0.5f*(float)sqrt(x1*x1 + y1*y1);
	   w  = c*(s1 + s2)/d;
	   af[i-1].x = c*d;
	   j  = i + iniva_tng - 3;
	   af[i-1].A = (float)(1.4*alph[j]*pow(w, beta[j])); /* 40% ad hoc increase */
	}
/* lips (2 sections with the equel length) */
	af[np-2].A = af[np-1].A = pi * lip_h * lip_w * cc;
	af[np-2].x = af[np-1].x = (float)(0.5 * (ivt[np-2].x - ivt[np-1].x) * c);

/* number of sections */
	*ns = np;

/* Check areas */
	for(i=0; i<*ns; i++)
	{  if(af[i].A <= 0.0) af[i].A = 0.0001f;
	   if(af[i].x <= 0.0) af[i].x = 0.01f;
	}
}

/*****
*	Function : appro_area_function
*	Note :	Approximate the area function from LAM, in which the section
*		length varies, by an area function with M tubes with a fixed
*		section length.  In the approximation, the fixed section
*		lenght is determined as dx = total length/M, so that the
*		total length remains invariant.
*****/

void	appro_area_function (
	short		ns1,	/* number of sections */
	area_function	*af1,	/* input area function with ns1 sections */
	short		ns2,	/* number of fixed sections */
	area_function	*af2 )	/* output areas function */
{
	float	dx;
	float	x, z1, z2, s1, s2;
	short	i, j;

/* compute dx */
	for(x=0., i=0; i<ns1; i++) x = x + af1[i].x;
	for(dx=x/ns2,i=0; i<ns2; i++) af2[i].x = dx;

/* approximation of areas */
	x = z2 = s2 = 0.;
	i = j = 0;
	while( i < ns2 )
	{  x += dx;
	   while( z2 <= x )
	   {  z1 = z2; z2 += af1[j].x;
	      s1 = s2; s2 += af1[j].x*af1[j].A;
	      if( ++j == ns1 ) break;
	   }
	   af2[i++].A = (s1 + (x - z1)*af1[j-1].A)/dx;
	   while( x+dx <= z2 )
	   {  af2[i++].A = af1[j-1].A;
	      x += dx;
	   }
	   s2 = (z2 - x)*af1[j-1].A;
	}

/*	Let the approximated lip area, af2[ns2-1].A, equal to the lip area*/
/*	of the original area function, af1[ns1-1].A.  This operation is   */
/*	necessary, since the averaged area, af2.A, tends to be too large  */
/*	in the case of the vowel such as [u], resulting in too high F2    */
/*	value.                                                            */

	af2[ns2-1].A = af1[ns1-1].A;
}

void print_af (short ns, area_function *af) {
    short i;
    
    printf("area_function = c(");
    for (i=0;i<ns;i++) { printf("%0.3f,",af[i].A); }
    printf(");\n");
}
/*****************( functions for reading the data file )******************/

/*****
*	Function : skiplines
*	Note :	Skip a specified number of lines in a text file.
*****/

void	skiplines( FILE *in, short nlines )
{
	short	i;
	for(i=0; i<nlines; i++) while( fgetc(in) != '\n' );
}


/*****
*	Function : read_model_spec
*	Note :	Reads specifications for a linear articulatory model from
*		a text file (PB1_spec.dat).  Only items needed for the
*		model are read.
*****/
void	read_model_spec( void )
{
	FILE	*in;
	short	nfs, nafs;
	short	i, j, dummy;
	char	vlab_dummy[3];

/* open file */
	if((in = fopen( "pb1_spec.dat", "rt")) == NULL)
	{  printf("Impossible to open the input file.\n");
	//	   getch();
	   exit(1);
	}

/* semi-polar coordinate specs. */
	skiplines( in, 9 );
	fscanf(in, "%hd %hd %hd %f %f %f %hd %hd\n",
		    &m1, &m2, &m3, &dl, &omega, &theta, &ix0, &iy0 );

	skiplines( in, 1 );
	fscanf(in, "%f %f\n", &TEKvt, &TEKlip);

	skiplines( in, 1 );
	for(i=0; i<m1+m2+m3; i++)
	   fscanf(in, "%hd %f %f\n", &dummy, &alph[i], &beta[i]);

/* Lip specifications */
	skiplines( in, 2 );
	fscanf(in, "%hd %hd %hd %hd %hd %hd\n",
	&nvrs_lip, &jaw_lip, &dummy, &dummy, &nfs, &nafs);
	if( jaw_lip != JAW || nfs < JAW+LIP )
	{  printf("Not enough factors in the lip spec..");
	   exit(1);
	}

	skiplines( in, nafs + 2 );
	for(i=0; i<nvrs_lip; i++) fscanf(in, "%s\n", vlab_dummy);

	skiplines( in, 1 );
	for(i=0; i<JAW+LIP; i++)  fscanf(in, "%s\n", flab_lip[i]);

	skiplines( in, 2 );
	for(i=0; i<nvrs_lip; i++) fscanf(in, "%f\n", &u_lip[i]);

	skiplines( in, 1 );
	for(i=0; i<nvrs_lip; i++) fscanf(in, "%f\n", &s_lip[i]);

	skiplines( in, 1 );
	fscanf(in, "%f %f\n", &inci_x, &inci_y);

	skiplines( in, 3 );
	for(i=0; i<nvrs_lip; i++)
	{  for(j=0; j<JAW+LIP; j++) fscanf(in, "%f\n", &A_lip[i][j]);
	   skiplines( in, 1 );
	}

	skiplines( in, 1 + nvrs_lip );

/* Tongue */
	skiplines( in, 2 );
	fscanf(in, "%hd %hd %hd %hd %hd %hd\n",
	&nvrs_tng, &jaw_tng, &iniva_tng, &lstva_tng, &nfs, &nafs);
	if( jaw_tng != JAW || nfs < JAW+TNG )
	{  printf("Not enough factors in the tongue spec..");
	   exit(1);
	}
	iniva_tng--;	/* coordinate address for C, now same as lable */
	lstva_tng--;

	skiplines( in, nafs + 2 );
	for(i=0; i<nvrs_tng; i++) fscanf(in, "%s\n", vlab_dummy);

	skiplines( in, 1 );
	for(i=0; i<JAW+TNG; i++)  fscanf(in, "%s\n", flab_tng[i]);

	skiplines( in, 2 );
	for(i=0; i<nvrs_tng; i++) fscanf(in, "%f\n", &u_tng[i]);

	skiplines( in, 1 );
	for(i=0; i<nvrs_tng; i++) fscanf(in, "%f\n", &s_tng[i]);

	skiplines( in, 1 );
	for(i=0; i<nvrs_tng; i++)
	{  for(j=0; j<JAW+TNG; j++) fscanf(in, "%f\n", &A_tng[i][j]);
	   skiplines( in, 1 );
	}

	skiplines( in, 1 + nvrs_tng );

/* Larynx */
	skiplines( in, 2 );
	fscanf(in, "%hd %hd %hd %hd %hd %hd\n",
	&nvrs_lrx, &jaw_lrx, &iniva_lrx, &lstva_lrx, &nfs, &nafs);
	if( jaw_lrx != JAW || nfs < JAW+LRX )
	{  printf("No enough factors in the larynx spec..");
	   exit(1);
	}

	skiplines( in, nafs + 2 );
	for(i=0; i<nvrs_lrx; i++) fscanf(in, "%s\n", vlab_dummy);

	skiplines( in, 1 );
	for(i=0; i<JAW+LRX; i++)  fscanf(in, "%s\n", flab_lrx[i]);

	skiplines( in, 2 );
	for(i=0; i<nvrs_lrx; i++) fscanf(in, "%f\n", &u_lrx[i]);

	skiplines( in, 1 );
	for(i=0; i<nvrs_lrx; i++) fscanf(in, "%f\n", &s_lrx[i]);

	skiplines( in, 1 );
	for(i=0; i<nvrs_lrx; i++)
	{  for(j=0; j<JAW+LRX; j++) fscanf(in, "%f\n", &A_lrx[i][j]);
	   skiplines( in, 1 );
	}

	skiplines( in, 1 + nvrs_lrx );

/* Wall */
	skiplines( in, 2 );
	fscanf(in, "%hd %hd %hd %hd %hd %hd\n",
	&nvrs_wal, &jaw_wal, &iniva_wal, &lstva_wal, &nfs, &nafs);
	iniva_wal--;	/* coordinate address for C, now same as lable */
	lstva_wal--;

	skiplines( in, nafs + 2 );
	for(i=0; i<nvrs_wal; i++) fscanf(in, "%s\n", vlab_dummy);

	skiplines( in, 3 );
	for(i=0; i<nvrs_wal; i++) fscanf(in, "%f\n", &u_wal[i]);

	fclose(in);
}

/*************************( plotting functions )**************************/
/*****
*	Function : convert_scale
*	Note :	Convert the cm-TEK unit mapping coefs., TEKvt and TEKlip,
*		and change the center coordinate (ix0, iy0), in order to
*		plot the genetated VT contours inside the specified
*		viewport.  Mean and standard deviation in TEK unit are
*		converted into the viewport unit. A viewport mapping coef.,
*		vp_map, is defined as its width corresponds to vp_width_cm
*		cm.
*****/
#define DWIDTH	10*20		/* display width in pixels */
#define DHEIGHT	10*20		/* display height in pixels */

void	convert_scale ( void )
{
	short	i;
	short done_flag=0;

/* viewport-to-cm scale factors (mapping coefficient) */
	vp_map = vp_width_cm/(DWIDTH/10);

/* Modify TEK-to-cm to TEK-to-viewport scale factor */
	TEKvt  = TEKvt*vp_map;
	if( TEKlip == 0. ) TEKlip = TEKvt;
	else               TEKlip = TEKlip*vp_map;

/* convert absolute mean maxillary incisor position to that relative
   to the semi-polar coordinate center, and to in viewport unit */
	inci_x = (inci_x - ix0)/TEKvt;
	inci_y = (inci_y - iy0)/TEKvt;

/* convert din in cm to in vp-unit */
	inci_lip_vp = inci_lip/vp_map;

/* convert TEK-unit to the viewport unit */
	for(i=0; i<2; i++)              /* jaw and lip protrusion */
	{  u_lip[i] = u_lip[i]/TEKvt;
	   s_lip[i] = s_lip[i]/TEKvt;
	}
	for(i=2; i<nvrs_lip; i++)	/* front lip height and width */
	{  u_lip[i] = u_lip[i]/TEKlip;
	   s_lip[i] = s_lip[i]/TEKlip;
	}
	for(i=0; i<nvrs_tng; i++)	/* tongue profile */
	{  u_tng[i] = u_tng[i]/TEKvt;
	   s_tng[i] = s_tng[i]/TEKvt;
	}
	for(i=0; i<nvrs_lrx; i++)	/* larynx profile */
	{  u_lrx[i] = u_lrx[i]/TEKvt;
	   s_lrx[i] = s_lrx[i]/TEKvt;
	}

	if (done_flag == 0) {
	for(i=0; i<nvrs_wal; i++)	/* vt rear wall profile */
	   u_wal[i] = u_wal[i]/TEKvt;
	done_flag=1;
	}

/* new coordinate center in the viewport */
	ix0 = 0.6*(DWIDTH/10);
	iy0 = 0.6*(DHEIGHT/10);
}



void	print_lam ( void )
{
	float	theta,x,y;
	short	i;
    
    printf("ivt_x = c(");
    for (i=0;i<np;i++) printf("%0.3f,",ivt[i].x);
    printf(");\n");
    printf("ivt_y = c(");
    for (i=0;i<np;i++) printf("%0.3f,",ivt[i].y);
    printf(");\n");
    printf("evt_x = c(");
    for (i=0;i<np;i++) printf("%0.3f,",evt[i].x);
    printf(");\n");
    printf("evt_y = c(");
    for (i=0;i<np;i++) printf("%0.3f,",evt[i].y);
    printf(");\n");
    
    /* lip frontal shape */
    theta = 2.*pi/20.;
	
    printf("lip_x = c(");
	for(i=0; i<=20; i++){
        x = lip_w*cos(i*theta);
        printf("%0.3f,",x);
    }
    printf(");\n");
    
    printf("lip_y = c(");
	for(i=0; i<=20; i++){
        y = lip_h*sin(i*theta);
        printf("%0.3f,",y);
    }
    printf(");\n");
}


#ifdef PLOTME

/*****
*	Function : plot_semi_polar
*	Note :	Plot the semi_polar coordinate in a current viewport.
*****/
void	plot_semi_polar ( short nvp )
{
	extern	viewport  vp[NVPS];
	float	vp_height;
	short	i;

	vp_height = pixtolo_y(vp[nvp].h);
	for(i=0; i<m1+m2+m3; i++)
	{  mov( igd[i].x, vp_height - igd[i].y );
	   drw( egd[i].x, vp_height - egd[i].y );
	}
}


/*****
*	Function : plot_lam
*	Note :	Plot a VT profile generated by "lam".  In order to erase
*		the plotted contours, the points are stored in arrays.
*		This operation is necessary only to create an animation
*		effect.
*****/

void	plot_lam ( short nvp )
{
	extern	viewport  vp[NVPS];
	float	vp_height;
	short	lip_x0, lip_y0;		/* fronttal-lip position in pixels */
	float	theta;
	short	i;

/* inside contour */

	for(i=0; i<np; i++)
	{  ivt_pix[i].x = lotopix_x(ivt[i].x);
	   ivt_pix[i].y = vp[nvp].h - lotopix_y(ivt[i].y);
	   if(i == 0) moveto(ivt_pix[i].x, ivt_pix[i].y);
	   else       lineto(ivt_pix[i].x, ivt_pix[i].y);
	}

/* outside contour */

	for(i=0; i<np; i++)
	{  evt_pix[i].x = lotopix_x(evt[i].x);
	   evt_pix[i].y = vp[nvp].h - lotopix_y(evt[i].y);
	   if(i == 0) moveto(evt_pix[i].x, evt_pix[i].y);
	   else       lineto(evt_pix[i].x, evt_pix[i].y);
	}

/* lip frontal shape */

	lip_x0 = .25*vp[nvp].w;
	lip_y0 = .50*vp[nvp].h;
	theta = 2.*pi/20.;
	for(i=0; i<=20; i++)
	{  lip_pix[i].x = lotopix_x(lip_w*cos(i*theta))+lip_x0;
	   lip_pix[i].y = lotopix_y(lip_h*sin(i*theta))+lip_y0;
	   if(i == 0) moveto(lip_pix[i].x, lip_pix[i].y);
	   else       lineto(lip_pix[i].x, lip_pix[i].y);
	}
}

/*****
*	Function : deplot_lam
*	Note :	Erase a current VT profile generated by "lam".
*****/

void	deplot_lam ( void )
{
	short	oldcolor;
	short	i;

	oldcolor = getcolor();

	setcolor( getbkcolor() );

	moveto(ivt_pix[0].x, ivt_pix[0].y);
	for(i=1; i<np; i++) lineto(ivt_pix[i].x, ivt_pix[i].y);

	moveto(evt_pix[0].x, evt_pix[0].y);
	for(i=1; i<np; i++) lineto(evt_pix[i].x, evt_pix[i].y);

	moveto( lip_pix[0].x, lip_pix[0].y);
	for(i=1; i<=20; i++)
	   lineto( lip_pix[i].x, lip_pix[i].y);

	setcolor(oldcolor);
}

#endif
