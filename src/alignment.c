/***************************************************************************
 *   Copyright (C) 2007-2018 by Georg Hennig                               *
 *   georg.hennig@web.de                                                   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

/*
gimp alignment plug-in
(C) Georg Hennig <georg.hennig@web.de>
Take one ore two selections, and aligns the layers using these selections.
Alignment criteria are 2D cross correlation (FFT or straight forward),
gauss fit (for stars) or center of brightness.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multiroots.h>

#include <fftw3.h>

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-align-layers"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

/*
parameters
*/
typedef struct
{
	gint32 alignment_method;
	gint32 search_radius;
	gint32 cross_correlation;

	gint32 scale_percent;

	gboolean visible_only;
	gboolean invisible_remove;

	gboolean trim_image;
	gboolean two_star;
	gboolean subpixel_resolution;
} tparameter;

enum ALIGNMENT_METHOD
{
	GAUSS_FIT = 0,
	CENTER_OF_BRIGHTNESS,
	CROSS_CORRELATION,
	CROSS_CORRELATION_FFT
};

enum QUALITY_METHOD
{
	NONE = 0,
	CONTRAST,
	FIDELITY,
	ALIGNMENT_RELATED
};

static tparameter parameters =
{
	CROSS_CORRELATION_FFT,
	3,
	400,
	100,
	TRUE,
	FALSE,
	TRUE,
	TRUE,
	TRUE
};

const int NUMBER_COEFFICIENTS = 9;

/*
Prototypes
*/

/*
communication to the GIMP
*/
static void query( void );

static void run( const gchar *name, gint nparams, const GimpParam *param,
	gint *nreturn_vals, GimpParam **return_vals );

/*
center methods
*/
void get_center( gint32 layer, gint alignment_method, gint sel_pos_x, gint sel_pos_y, gint sel_width, gint sel_height,
	gint radius, guchar *reference_data, const gint fit_bpp, const gint cross_correlation,
	gint *move_x, gint *move_y, gdouble *quality );

/*
quality
*/
void get_quality( gint32 layer, gint quality_method, gint sel_pos_x, gint sel_pos_y, gint sel_width, gint sel_height,
	gdouble *quality );

gint get_dynamic_range( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gdouble *dynamic_range );
gint get_fidelity( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gdouble *fidelity );

/*gauss fit*/
gint gauss_f( const gsl_vector *p, void *data, gsl_vector *f );
gint gauss_df( const gsl_vector *p, void *data, gsl_matrix *J );
gint gauss_fdf( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J );
gint get_gauss_fit( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	gdouble *pos_x, gdouble *pos_y, gdouble *sigma_x, gdouble *sigma_y );

/*center of brightness*/
gint get_center_of_brightness( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	gdouble *pos_x, gdouble *pos_y, gdouble *maximum );

/*cross correlation*/
gint get_cross_correlation( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	const guchar *reference_data, const gint fit_x, const gint fit_y, const gint fit_bpp, const gint cross_correlation,
	gdouble *pos_x, gdouble *pos_y, const gboolean subpixel_resolution, gdouble *maximum );

/*cross correlation using fftw*/
gint get_cross_correlation_fft( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	const guchar *reference_data, const gint fit_x, const gint fit_y, const gint fit_bpp, const gint cross_correlation,
	gdouble *pos_x, gdouble *pos_y, const gboolean subpixel_resolution, gdouble *maximum );

/*
Align layers
*/
static void align_layers( gint32 image_id );

/*
user interface
*/
static gint dialog();


/*
Variables
*/

/*
PLUG_IN_INFO
*/
GimpPlugInInfo PLUG_IN_INFO =
{
	NULL, /* init_proc  */
	NULL, /* quit_proc  */
	query,/* query_proc */
	run   /* run_proc   */
};

/*
access to levels
*/
guchar *level_buffer;
GtkWidget *level_preview;
gdouble level_log_max;


/*
procedures
*/
MAIN()

/*
communication to the GIMP
*/
static void query( void )
{
	static GimpParamDef params[] =
	{
		{ GIMP_PDB_INT32, "run_mode", "Interactive,non-interactive" },
		{ GIMP_PDB_IMAGE, "image_id", "Input image" },
		{ GIMP_PDB_DRAWABLE, "drawable", "Input drawable" },
		{ GIMP_PDB_INT32, "alignment_method", "Layer alignment Method (0,1,2,3)" },
		{ GIMP_PDB_INT32, "search_radius", "Search radius in px (1-100)" },
		{ GIMP_PDB_INT32, "cross_correlation", "Search area in % of selected area width/height (100-1000)" },
		{ GIMP_PDB_INT32, "scale_percent", "Scale before processing in % (100-400)" },
		{ GIMP_PDB_INT32, "visible_only", "Process visible layers only" },
		{ GIMP_PDB_INT32, "invisible_remove", "Remove invisible layers" },
		{ GIMP_PDB_INT32, "trim_image", "Trim image to overlap area of layers" },
		{ GIMP_PDB_INT32, "two_star", "Two star alignment (de-rotate)" },
		{ GIMP_PDB_INT32, "subpixel_resolution", "Subpixel resolution when using cross correlation" }
	};

	/*  Initialize i18n support  */
	bindtextdomain( GETTEXT_PACKAGE, gimp_locale_directory() );
#ifdef HAVE_BIND_TEXTDOMAIN_CODESET
	bind_textdomain_codeset( GETTEXT_PACKAGE, "UTF-8" );
#endif
	textdomain( GETTEXT_PACKAGE );

	static GimpParamDef *return_vals  = NULL;
	static gint nparams = sizeof( params )/sizeof( params[0] );
	static gint nreturn_vals = 0;

	gimp_install_procedure( PLUG_IN_NAME,
		_("Move and rotate layers so that they overlap correctly"),
		_("This plug-in translates and rotates layers so they fit best."),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Align Layers"),
		"RGB*,GRAY*",
		GIMP_PLUGIN,
		nparams,
		nreturn_vals,
		params,
		return_vals );

	gimp_plugin_menu_register( PLUG_IN_NAME, _("<Image>/Filters/Astronomy") );
}

static void run( const gchar *name, gint nparams, const GimpParam  *param,
	gint *nreturn_vals, GimpParam **return_vals )
{
	static GimpParam values[1];
	GimpPDBStatusType status;
	GimpRunMode run_mode;

	status = GIMP_PDB_SUCCESS;
	run_mode = param[0].data.d_int32;
	*nreturn_vals = 1;
	*return_vals = values;

	/*  Initialize i18n support  */
	bindtextdomain( GETTEXT_PACKAGE, gimp_locale_directory() );
#ifdef HAVE_BIND_TEXTDOMAIN_CODESET
	bind_textdomain_codeset( GETTEXT_PACKAGE, "UTF-8" );
#endif
	textdomain( GETTEXT_PACKAGE );

	values[0].type = GIMP_PDB_STATUS;
	values[0].data.d_status = status;

	switch(run_mode)
	{
		case GIMP_RUN_INTERACTIVE:
			gimp_get_data( PLUG_IN_NAME, &parameters );

			if ( !dialog() ) return;

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 12 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				parameters.alignment_method = param[3].data.d_int32;
				parameters.search_radius = param[4].data.d_int32;
				parameters.cross_correlation = param[5].data.d_int32;
				parameters.scale_percent = param[6].data.d_int32;
				parameters.visible_only = param[7].data.d_int32;
				parameters.invisible_remove = param[8].data.d_int32;
				parameters.trim_image = param[9].data.d_int32;
				parameters.two_star = param[10].data.d_int32;
				parameters.subpixel_resolution = param[11].data.d_int32;
			}
			break;
		case GIMP_RUN_WITH_LAST_VALS:
			gimp_get_data( PLUG_IN_NAME, &parameters );
			break;
		default:
			status = GIMP_PDB_CALLING_ERROR;
			break;
	}

	if( status == GIMP_PDB_SUCCESS )
	{
		align_layers( param[1].data.d_image );
	}

	values[0].data.d_status = status;
}

typedef struct
{
	size_t n;
	gdouble *x;
	gdouble *y;
	gdouble *z;
	gdouble *sigma;
} fit_data;

gint gauss_f( const gsl_vector *p, void *data, gsl_vector *f )
{
	size_t n = ((fit_data *)data)->n;
	gdouble *x = ((fit_data *) data)->x;
	gdouble *y = ((fit_data *) data)->y;
	gdouble *z = ((fit_data *) data)->z;
	gdouble *sigma = ((fit_data *) data)->sigma;

	gdouble A = gsl_vector_get( p, 0 );
	gdouble x0 = gsl_vector_get( p, 1 );
	gdouble sigma_x = gsl_vector_get( p, 2 );
	gdouble y0 = gsl_vector_get( p, 3 );
	gdouble sigma_y = gsl_vector_get( p, 4 );
	gdouble b = gsl_vector_get( p, 5 );

	size_t i;

	for ( i=0; i<n; i++ )
	{
		/* Zi = A * exp(-0.5*(x-x0)^2/sigma_x^2) * exp(-0.5*(y-y0)^2/sigma_y^2) + b */
		gdouble Zi = A * exp( -0.5*(x[i]-x0)*(x[i]-x0)/(sigma_x*sigma_x) ) * exp( -0.5*(y[i]-y0)*(y[i]-y0)/(sigma_y*sigma_y) ) + b;
		gsl_vector_set( f, i, (Zi - z[i])/sigma[i] );
	}

	return GSL_SUCCESS;
}

gint gauss_df( const gsl_vector *p, void *data, gsl_matrix *J )
{
	size_t n = ((fit_data *) data)->n;
	gdouble *x = ((fit_data *) data)->x;
	gdouble *y = ((fit_data *) data)->y;
	gdouble *sigma = ((fit_data *) data)->sigma;

	gdouble A = gsl_vector_get( p, 0 );
	gdouble x0 = gsl_vector_get( p, 1 );
	gdouble sigma_x = gsl_vector_get( p, 2 );
	gdouble y0 = gsl_vector_get( p, 3 );
	gdouble sigma_y = gsl_vector_get( p, 4 );

	size_t i;

	for ( i=0; i<n; i++ )
	{
		/* Jacobian matrix J(i,j) = dfi / dxj, */
		/* where fi = (Yi - yi)/sigma[i],      */
		/*       Yi = A * exp(-lambda * i) + b  */
		/* and the xj are the parameters (A,lambda,b) */
		gdouble e_x = exp( -0.5*(x[i]-x0)*(x[i]-x0)/(sigma_x*sigma_x) );
		gdouble e_y = exp( -0.5*(y[i]-y0)*(y[i]-y0)/(sigma_y*sigma_y) );
		gsl_matrix_set( J, i, 0, e_x * e_y / sigma[i] );
		gsl_matrix_set( J, i, 1, A * e_x * e_y * ( -(x[i]-x0)/(sigma_x*sigma_x) ) / sigma[i] );
		gsl_matrix_set( J, i, 2, A * e_x * e_y * ( (x[i]-x0)*(x[i]-x0)/(sigma_x*sigma_x*sigma_x) ) / sigma[i] );
		gsl_matrix_set( J, i, 3, A * e_x * e_y * ( -(y[i]-y0)/(sigma_y*sigma_y) ) / sigma[i] );
		gsl_matrix_set( J, i, 4, A * e_x * e_y * ( (y[i]-y0)*(y[i]-y0)/(sigma_y*sigma_y*sigma_y) ) / sigma[i] );
		gsl_matrix_set( J, i, 5, 1./sigma[i] );
	}

	return GSL_SUCCESS;
}

gint gauss_fdf( const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J )
{
	gauss_f( x, data, f );
	gauss_df( x, data, J );

	return GSL_SUCCESS;
}

gint get_gauss_fit( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	gdouble *pos_x, gdouble *pos_y, gdouble *sigma_x, gdouble *sigma_y )
{
	const gsl_multifit_fdfsolver_type *T;
	gsl_multifit_fdfsolver *s;

	const gint n = X*Y;
	const gint p = 6;

	gsl_matrix *covar = gsl_matrix_alloc( p, p );
	gdouble x[n], y[n], z[n], sigma[n];
	fit_data d = { n, x, y, z, sigma };
	gsl_multifit_function_fdf f;

	/* This is the data to be fitted */
	gint i, j;
	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			x[i*X+j] = j;
			y[i*X+j] = i;
			z[i*X+j] = data[bpp*(i*X+j)+bit];
			sigma[i*X+j] = sqrt( z[i*X+j]+1 ); /* Photon noise is sqrt(N) */
		}
	}

	/* Set initial values of: A, x0, sigma_x, y0, sigma_y, b */
	/* Remember: Zi = A * exp(-0.5*(x-x0)^2/sigma_x^2) * exp(-0.5*(y-y0)^2/sigma_y^2) + b */
	/* Getting a good (?) guess for x0, y0, A+b */
	gdouble maximum_x = 0., maximum_y = 0.;
	gdouble maximum = 0.;
	get_center_of_brightness( data, bit, bpp, X, Y, radius, &maximum_x, &maximum_y, &maximum );
	/* The mean value of the corners is taken as a guess for b */
	gdouble b_guessed = (gdouble)((data[bit]+data[bpp*(X-1)+bit]+data[bpp*((Y-1)*X)+bit]+data[bpp*(X*Y-1)+bit])/4);
	/* Trying to guess the standard deviations */
	gdouble A_over_e = (maximum-b_guessed)/exp(1.);
	gdouble sigma_x_guessed = 1., sigma_y_guessed = 1.;
	for ( i=0; i<Y; i++ )
	{
		if ( maximum_y - i >= 0 && maximum_y + i < Y )
		{
			if ( 0.5*( data[bpp*(ROUND(maximum_y-i)*X+ROUND(maximum_x))+bit] + data[bpp*(ROUND(maximum_y+i)*X+ROUND(maximum_x))+bit] ) < A_over_e )
			{
				sigma_y_guessed = i;
				break;
			}
		}
	}
	for ( j=0; j<X; j++ )
	{
		if ( maximum_x - j >= 0 && maximum_x + j < X )
		{
			if ( 0.5*( data[bpp*(ROUND(maximum_y)*X+ROUND(maximum_x-j))+bit] + data[bpp*(ROUND(maximum_y)*X+ROUND(maximum_x+j))+bit] ) < A_over_e )
			{
				sigma_x_guessed = j;
				break;
			}
		}
	}
	gdouble x_init[6] = { maximum-b_guessed, maximum_x, sigma_x_guessed, maximum_y, sigma_y_guessed, b_guessed };
	printf( "alignment:\n  gauss fit inital values: A=%.5f, x0=%.5f, sigma_x=%.5f, y0=%.5f, sigma_y=%.5f, b=%.5f\n",
		x_init[0], x_init[1], x_init[2], x_init[3], x_init[4], x_init[5] );

	gsl_vector_view view = gsl_vector_view_array( x_init, p );

	f.f = &gauss_f;
	f.df = &gauss_df;
	f.fdf = &gauss_fdf;
	f.n = n;
	f.p = p;
	f.params = &d;

	T = gsl_multifit_fdfsolver_lmsder;
	s = gsl_multifit_fdfsolver_alloc( T, n, p );
	gsl_multifit_fdfsolver_set( s, &f, &view.vector );

	gint status, iter = 0;
	do
	{
		iter++;
		status = gsl_multifit_fdfsolver_iterate( s );

		printf( "alignment:\n  gauss fit status = %s\n", gsl_strerror( status ) );

		if ( status ) break;

		status = gsl_multifit_test_delta( s->dx, s->x, 1e-8, 1e-8 );
	} while (status == GSL_CONTINUE && iter < 5000);

	gsl_matrix* J = gsl_matrix_alloc(n, p);
 	gsl_multifit_fdfsolver_jac(s, J);
 	gsl_multifit_covar(J, 0.0, covar);
 	gsl_matrix_free(J);

	gdouble chi = gsl_blas_dnrm2( s->f );
	gdouble dof = n - p;
	gdouble c = GSL_MAX_DBL( 1, chi/sqrt( dof ) );

	*pos_x = (gdouble)gsl_vector_get( s->x, 1 );
	*pos_y = (gdouble)gsl_vector_get( s->x, 3 );
	*sigma_x = (gdouble)gsl_vector_get( s->x, 2 );
	*sigma_y = (gdouble)gsl_vector_get( s->x, 4 );

	printf( "alignment:\n  chisq/dof = %g\n",  pow( chi, 2.0 )/dof );

	printf( "alignment:\n  A       = %.5f +/- %.5f\n", gsl_vector_get( s->x, 0 ), c*sqrt( gsl_matrix_get( covar, 0, 0 ) ) );
	printf( "alignment:\n  x0      = %.5f +/- %.5f\n", gsl_vector_get( s->x, 1 ), c*sqrt( gsl_matrix_get( covar, 1, 1 ) ) );
	printf( "alignment:\n  sigma_x = %.5f +/- %.5f\n", gsl_vector_get( s->x, 2 ), c*sqrt( gsl_matrix_get( covar, 2, 2 ) ) );
	printf( "alignment:\n  y0      = %.5f +/- %.5f\n", gsl_vector_get( s->x, 3 ), c*sqrt( gsl_matrix_get( covar, 3, 3 ) ) );
	printf( "alignment:\n  sigma_y = %.5f +/- %.5f\n", gsl_vector_get( s->x, 4 ), c*sqrt( gsl_matrix_get( covar, 4, 4 ) ) );
	printf( "alignment:\n  b       = %.5f +/- %.5f\n", gsl_vector_get( s->x, 5 ), c*sqrt( gsl_matrix_get( covar, 5, 5 ) ) );

	printf( "alignment:\n  gauss fit status = %s\n", gsl_strerror( status ) );

	gsl_multifit_fdfsolver_free( s );
	gsl_matrix_free( covar );

	return 0;
}

gint get_center_of_brightness( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	gdouble *pos_x, gdouble *pos_y, gdouble *maximum )
{
	if ( radius < 0 ) radius = 0;
	if ( radius > X ) radius = X;
	if ( radius > Y ) radius = Y;

	/* Calculate the mean value of X*Y boxes with the size mean_diameter*mean_diameter */
	gint maximum_x = 0, maximum_y = 0;
	gdouble maximum_value = 0.;
	gdouble mean_value = 0.;
	gint mean_value_counter;
	/* This is the data to be fitted */
	gint i, j, k, l;
	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			mean_value = 0.;
			mean_value_counter = 0;
			for ( k=i-radius; k<i+radius; k++ )
			{
				for( l=j-radius; l<j+radius; l++ )
				{
					if ( k >= 0 && k < Y && l >= 0 && l < X )
					{
						mean_value += data[bpp*(k*X+l)+bit];
						mean_value_counter++;
					}
				}
			}

			if ( maximum_value < mean_value/mean_value_counter )
			{
				maximum_value = mean_value/mean_value_counter;
				maximum_x = j;
				maximum_y = i;
			}
			
		}
	}

	*pos_x = (gdouble)maximum_x;
	*pos_y = (gdouble)maximum_y;
	*maximum = (gdouble)maximum_value;

	return 0;
}

void fit_polynom( const double *x, const double *y, const double *z, const double *err, int start, int n,
	double *coeffs, int start_coeffs )
{
	gsl_matrix *X = NULL, *cov = NULL;
	gsl_vector *Y = NULL, *W = NULL, *c = NULL;

	int i;
	double chisq;

	X = gsl_matrix_alloc( n, NUMBER_COEFFICIENTS );
	Y = gsl_vector_alloc( n );
	W = gsl_vector_alloc( n );

	c = gsl_vector_alloc( NUMBER_COEFFICIENTS );
	cov = gsl_matrix_alloc( NUMBER_COEFFICIENTS, NUMBER_COEFFICIENTS );

	for ( i=0; i<n; i++ )
	{
		gsl_matrix_set( X, i, 0, 1.0 );
		gsl_matrix_set( X, i, 1, x[i+start] );
		gsl_matrix_set( X, i, 2, x[i+start]*x[i+start] );
		gsl_matrix_set( X, i, 3, y[i+start] );
		gsl_matrix_set( X, i, 4, y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 5, x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 6, x[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 7, x[i+start]*x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 8, x[i+start]*x[i+start]*y[i+start]*y[i+start] );

		gsl_vector_set( Y, i, z[i+start] );
		gsl_vector_set( W, i, 1.0/(err[i+start]*err[i+start]) );
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc( n, NUMBER_COEFFICIENTS );
	gsl_multifit_wlinear( X, W, Y, c, cov, &chisq, work );
	gsl_multifit_linear_free( work );

/*	printf(
		"alignment: subpixel resolution polynom: f(x,y)=%g %+g*x %+g*x**2 %+g*y %+g*y**2 %+g*x*y %+g*x*y**2 %+g*x**2*y %+g*x**2*y**2\n",
		gsl_vector_get( c, 0 ), gsl_vector_get( c, 1 ), gsl_vector_get( c, 2 ), gsl_vector_get( c, 3 ),
		gsl_vector_get( c, 4 ), gsl_vector_get( c, 5 ), gsl_vector_get( c, 6 ), gsl_vector_get( c, 7 ),
		gsl_vector_get( c, 8 ) );

	printf ( "background_gradient: covariance matrix:\n" );
	printf ( "  [ %+.5e, %+.5e, %+.5e  \n", gsl_matrix_get( cov, 0, 0 ),
		gsl_matrix_get( cov, 0, 1 ), gsl_matrix_get( cov, 0, 2 ) );
	printf ( "    %+.5e, %+.5e, %+.5e  \n", gsl_matrix_get( cov, 1, 0 ),
		gsl_matrix_get( cov, 1, 1 ), gsl_matrix_get( cov, 1, 2 ) );
	printf ( "    %+.5e, %+.5e, %+.5e ]\n", gsl_matrix_get( cov, 2, 0 ),
		gsl_matrix_get( cov, 2, 1 ), gsl_matrix_get( cov, 2, 2 ) );
	printf ( "background_gradient: chisq = %g\n", chisq );*/

	for ( i=0; i<NUMBER_COEFFICIENTS; i++ )
	{
		coeffs[start_coeffs+i] = gsl_vector_get( c, i );
	}

	gsl_matrix_free( X );
	gsl_vector_free( Y );
	gsl_vector_free( W );
	gsl_vector_free( c );
	gsl_matrix_free( cov );
}

inline double polynom_value( const double x, const double y, double *coeffs, int start_coeffs )
{
	return coeffs[start_coeffs+0] + coeffs[start_coeffs+1]*x + coeffs[start_coeffs+2]*x*x +
		coeffs[start_coeffs+3]*y + coeffs[start_coeffs+4]*y*y +
		coeffs[start_coeffs+5]*x*y + coeffs[start_coeffs+6]*x*y*y +
		coeffs[start_coeffs+7]*x*x*y + coeffs[start_coeffs+8]*x*x*y*y;
}

gint get_cross_correlation( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	const guchar *reference_data, const gint fit_x, const gint fit_y, const gint fit_bpp, const gint cross_correlation,
	gdouble *pos_x, gdouble *pos_y, const gboolean subpixel_resolution, gdouble *maximum )
{
	gint i, j, k, l;
	gdouble cross_correlation_data = 0., cross_correlation_maximum = 0.;
	gdouble fit_square = 0., data_square = 0.;
	guchar data_current, fit_current;
	gint cross_correlation_x = 0, cross_correlation_y = 0;

	for ( i=fit_y/2; i<=(Y-fit_y/2); i++ )
	{
		for ( j=fit_x/2; j<=(X-fit_x/2); j++ ) /* calculate the cross correlation at each point inside this rectangle */
		{
			fit_square = 0.;
			data_square = 0.;
			cross_correlation_data = 0.;
			for ( k=-fit_y/2; k<fit_y/2; k++ )
			{
				for ( l=-fit_x/2; l<fit_x/2; l++ )
				{
					data_current = data[bpp*((i+k)*X+j+l)+bit];
					fit_current = reference_data[fit_bpp*((k+fit_y/2)*fit_x+l+fit_x/2)+bit];

					cross_correlation_data += fit_current*data_current;
				}
			}

			if ( cross_correlation_maximum < cross_correlation_data )
			{
				/* coordinate system of the fit layer */
				cross_correlation_maximum = cross_correlation_data;
				cross_correlation_x = j-X/2;
				cross_correlation_y = i-Y/2;
			}
		}
	}

	*pos_x = (gdouble)cross_correlation_x;
	*pos_y = (gdouble)cross_correlation_y;
	*maximum = (gdouble)cross_correlation_maximum;

	return 0;
}

gint get_cross_correlation_fft( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gint radius,
	const guchar *reference_data, const gint fit_x, const gint fit_y, const gint fit_bpp, const gint cross_correlation,
	gdouble *pos_x, gdouble *pos_y, const gboolean subpixel_resolution, gdouble *maximum )
{
	gint i, j;
	gdouble cross_correlation_maximum = 0.;
	gdouble cross_correlation_x = 0., cross_correlation_y = 0.;

	/* Fill fftw_real_input with the normalized feature, and transform it. */
	fftw_plan p;
	double *fftw_real_input;
	fftw_real_input = g_malloc( sizeof(double)*X*Y );

	fftw_complex *fftw_complex_feature;
	fftw_complex_feature = fftw_malloc( sizeof(fftw_complex)*X*Y );

	/* Swapping X and Y, because fftw has a different definition*/
	p = fftw_plan_dft_r2c_2d( Y, X, fftw_real_input, fftw_complex_feature, FFTW_ESTIMATE );

	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			if ( i<fit_y && j<fit_x ) fftw_real_input[i*X+j] = (double)reference_data[fit_bpp*(i*fit_x+j)+bit];
			else fftw_real_input[i*X+j] = 0.;
		}
	}

	fftw_execute( p );

	fftw_complex *fftw_complex_output;
	fftw_complex_output = fftw_malloc( sizeof(fftw_complex)*X*Y );

	p = fftw_plan_dft_r2c_2d( Y, X, fftw_real_input, fftw_complex_output, FFTW_ESTIMATE );

	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			fftw_real_input[i*X+j] = (double)data[bpp*(i*X+j)+bit];
		}
	}

	fftw_execute( p );

	p = fftw_plan_dft_c2r_2d( Y, X, fftw_complex_output, fftw_real_input, FFTW_ESTIMATE );

	fftw_complex product;
	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			product[0] = fftw_complex_output[i*X+j][0]*fftw_complex_feature[i*X+j][0] - fftw_complex_output[i*X+j][1]*(-fftw_complex_feature[i*X+j][1]);
			product[1] = fftw_complex_output[i*X+j][0]*(-fftw_complex_feature[i*X+j][1]) + fftw_complex_output[i*X+j][1]*fftw_complex_feature[i*X+j][0];

			fftw_complex_output[i*X+j][0] = product[0];
			fftw_complex_output[i*X+j][1] = product[1];
		}
	}

	fftw_execute( p );

	for ( i=0; i<Y; i++ )
	{
		for ( j=0; j<X; j++ )
		{
			if ( fftw_real_input[i*X+j] > cross_correlation_maximum )
			{
				cross_correlation_maximum = fftw_real_input[i*X+j];
				cross_correlation_x = j - X/2 + fit_x/2;
				cross_correlation_y = i - Y/2 + fit_y/2;
			}
		}
	}

	if ( subpixel_resolution && fit_x/2 >= 1 && fit_y/2 >=1 )
	{
		printf( "alignment: cross correlation (fft): Before subpixel resolution: x=%.5f, y=%.5f, max=%.5f\n",
			cross_correlation_x, cross_correlation_y, cross_correlation_maximum );

		double x[9] = { -1., 0., 1., -1., 0., 1., -1., 0., 1. };
		double y[9] = { 1., 1., 1., 0., 0., 0., -1., -1., -1. };
		double z[9] = {
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2+1)*X+((int)cross_correlation_x+X/2-fit_x/2)-1],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2+1)*X+((int)cross_correlation_x+X/2-fit_x/2)],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2+1)*X+((int)cross_correlation_x+X/2-fit_x/2)+1],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2)*X+((int)cross_correlation_x+X/2-fit_x/2)-1],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2)*X+((int)cross_correlation_x+X/2-fit_x/2)],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2)*X+((int)cross_correlation_x+X/2-fit_x/2)+1],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2-1)*X+((int)cross_correlation_x+X/2-fit_x/2)-1],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2-1)*X+((int)cross_correlation_x+X/2-fit_x/2)],
				fftw_real_input[((int)cross_correlation_y+Y/2-fit_y/2-1)*X+((int)cross_correlation_x+X/2-fit_x/2)+1]
			};
		double err[9] = { 1., 1., 1., 1., 1., 1., 1., 1., 1. };
		double coeffs[9];

		fit_polynom( x, y, z, err, 0, NUMBER_COEFFICIENTS, coeffs, 0 );

		double search = 1.;
		double positions[18] =
			{
				-1.,  1.,
				 0.,  1.,
				 1.,  1.,
				-1.,  0.,
				 0.,  0.,
				 1.,  0.,
				-1., -1.,
				 0., -1.,
				 1., -1.
			};
		int max = 4; /* position of maximum inside array positions */
		double maximum_value = polynom_value( 0., 0., coeffs, 0 );
		double max_x, max_y;
		for ( i=0; i<50; i++ ) /* 2^-50 = 10^-15 = double precision */
		{
			search *= 0.5;
			max_x = positions[2*max];
			max_y = positions[2*max+1];
			/* initialize new positions around maximum */
			positions[0]  = max_x-search;	positions[1]  = max_y+search;
			positions[2]  = max_x;				positions[3]  = max_y+search;
			positions[4]  = max_x+search;	positions[5]  = max_y+search;
			positions[6]  = max_x-search;	positions[7]  = max_y;
			positions[8]  = max_x;				positions[9]  = max_y;
			positions[10] = max_x+search;	positions[11] = max_y;
			positions[12] = max_x-search;	positions[13] = max_y-search;
			positions[14] = max_x;				positions[15] = max_y-search;
			positions[16] = max_x+search;	positions[17] = max_y-search;

			maximum_value = 0.;
			for ( j=0; j<9; j++ )
			{
				if ( polynom_value( positions[2*j], positions[2*j+1], coeffs, 0 ) > maximum_value )
				{
					maximum_value = polynom_value( positions[2*j], positions[2*j+1], coeffs, 0 );
					max = j;
				}
			}
		}
		printf( "alignment: cross correlation (fft): Subpixel resolution found x=%.5f, y=%.5f, max=%.5f\n", positions[2*max], positions[2*max+1],
			polynom_value( positions[2*max], positions[2*max+1], coeffs, 0 ) );

		if ( cross_correlation_maximum < polynom_value( positions[2*max], positions[2*max+1], coeffs, 0 ) )
		{
			cross_correlation_maximum = polynom_value( positions[2*max], positions[2*max+1], coeffs, 0 );
			cross_correlation_x += positions[2*max];
			cross_correlation_y += positions[2*max+1];
		}
	}

	*pos_x = cross_correlation_x;
	*pos_y = cross_correlation_y;
	*maximum = cross_correlation_maximum;

	return 0;
}

void get_center( gint32 layer, gint alignment_method, gint sel_pos_x, gint sel_pos_y, gint sel_width, gint sel_height,
	gint radius, guchar *reference_data, const gint fit_bpp, const gint cross_correlation,
	gint *move_x, gint *move_y, gdouble *quality )
{
	*move_x = 0;
	*move_y = 0;
	*quality = 0;

	gint fit_width = 0, fit_height = 0;
	if ( alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
	{
		fit_width = sel_width;
		fit_height = sel_height;
		sel_pos_x = (double)sel_pos_x - ((double)(cross_correlation-100)*sel_width)/200. + 0.5;
		sel_pos_y = (double)sel_pos_y - ((double)(cross_correlation-100)*sel_height)/200. + 0.5;

		sel_width = (double)(cross_correlation*sel_width)/100. + 0.5;
		sel_height = (double)(cross_correlation*sel_height)/100. + 0.5;
		if ( sel_width%2 ) sel_width--;
		if ( sel_height%2 ) sel_height--;
	}

	if ( sel_pos_x < 0 ) sel_pos_x = 0;
	if ( sel_pos_y < 0 ) sel_pos_y = 0;
	if ( sel_pos_x + sel_width > gimp_drawable_width( layer ) )  sel_width = gimp_drawable_width( layer ) - sel_pos_x;
	if ( sel_pos_y + sel_height > gimp_drawable_height( layer ) ) sel_height = gimp_drawable_height( layer ) - sel_pos_y;

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( layer ) )
	{
		case GIMP_RGB_IMAGE:
		case GIMP_RGBA_IMAGE:
			image_type = GIMP_RGB;
			layer_type = GIMP_RGB_IMAGE;
			bpp = 3;
			break;
		case GIMP_GRAY_IMAGE:
		case GIMP_GRAYA_IMAGE:
			image_type = GIMP_GRAY;
			layer_type = GIMP_GRAY_IMAGE;
			bpp = 1;
			break;
		case GIMP_INDEXED_IMAGE:
		case GIMP_INDEXEDA_IMAGE:
			return; /* Not allowed */
			break;
		default:
			g_message( _("Invalid color type!") );
			return;
	}

	guchar *data;
	data = malloc( gimp_drawable_bpp( layer )*sel_width*sel_height*sizeof(guchar) );

	GimpPixelRgn region_source;
	gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layer ),
		0, 0, gimp_drawable_width( layer ), gimp_drawable_height( layer ), FALSE, FALSE );

	gimp_pixel_rgn_get_rect( &region_source, data, sel_pos_x, sel_pos_y, sel_width, sel_height );

	gint i;
	gint ret;
	gdouble pos_x, pos_y, parameter_x, parameter_y;
	gdouble pos_x_mean, pos_y_mean, parameter_x_mean, parameter_y_mean;
	switch( alignment_method )
	{
		case GAUSS_FIT:
			pos_x_mean = 0.;
			pos_y_mean = 0.;
			parameter_x_mean = 0.;
			parameter_y_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				printf( "alignment: gauss fit of color %d\n", i );
				ret = get_gauss_fit( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, radius,
					&pos_x, &pos_y, &parameter_x, &parameter_y );
				pos_x_mean += pos_x;
				pos_y_mean += pos_y;
				parameter_x_mean += parameter_x;
				parameter_y_mean += parameter_y;
			}
			pos_x_mean /= bpp;
			pos_y_mean /= bpp;
			parameter_x_mean /= bpp;
			parameter_y_mean /= bpp;

			printf( "alignment: gauss fit results: x0=%.5f, y0=%.5f, sigma_x=%.5f, sigma_y=%.5f\n",
				pos_x_mean, pos_y_mean, parameter_x_mean, parameter_y_mean );

			*quality = (gdouble)( parameter_x_mean + parameter_y_mean +
				sqrt( ( parameter_x_mean - parameter_y_mean )*( parameter_x_mean - parameter_y_mean ) ) );
			*move_x = (gint)ROUND( pos_x_mean );
			*move_y = (gint)ROUND( pos_y_mean );

			break;

		case CENTER_OF_BRIGHTNESS:
			pos_x_mean = 0.;
			pos_y_mean = 0.;
			parameter_x_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				printf( "alignment: center of brightness of color %d\n", i );
				ret = get_center_of_brightness( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, radius,
					&pos_x, &pos_y, &parameter_x );
				pos_x_mean += pos_x;
				pos_y_mean += pos_y;
				parameter_x_mean += parameter_x;
			}
			pos_x_mean /= bpp;
			pos_y_mean /= bpp;
			parameter_x_mean /= bpp;

			printf( "alignment: center of brightness results: x0=%.5f, y0=%.5f, maximum=%.5f\n",
				pos_x_mean, pos_y_mean, parameter_x_mean );

			*quality = (gdouble)parameter_x_mean;
			*move_x = (gint)ROUND( pos_x_mean );
			*move_y = (gint)ROUND( pos_y_mean );

			break;

		case CROSS_CORRELATION:
			pos_x_mean = 0.;
			pos_y_mean = 0.;
			parameter_x_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				printf( "alignment: cross correlation of color %d: ", i );
				ret = get_cross_correlation( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, radius,
					reference_data, fit_width, fit_height, fit_bpp, cross_correlation,
					&pos_x, &pos_y, parameters.subpixel_resolution, &parameter_x );
				printf( "x0=%.5f, y0=%.5f, parameter=%.5f\n", pos_x, pos_y, parameter_x );
				pos_x_mean += pos_x;
				pos_y_mean += pos_y;
				parameter_x_mean += parameter_x;
			}
			pos_x_mean /= bpp;
			pos_y_mean /= bpp;
			parameter_x_mean /= bpp;

			printf( "alignment: cross correlation: x0=%.5f, y0=%.5f, maximum=%.5f\n",
				pos_x_mean, pos_y_mean, parameter_x_mean );

			*quality = (gdouble)parameter_x_mean;
			*move_x = (gint)ROUND( pos_x_mean );
			*move_y = (gint)ROUND( pos_y_mean );

			break;

		case CROSS_CORRELATION_FFT:
			pos_x_mean = 0.;
			pos_y_mean = 0.;
			parameter_x_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				printf( "alignment: cross correlation (fft) of color %d: ", i );
				ret = get_cross_correlation_fft( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, radius,
					reference_data, fit_width, fit_height, fit_bpp, cross_correlation,
					&pos_x, &pos_y, parameters.subpixel_resolution, &parameter_x );
				printf( "x0=%.5f, y0=%.5f, parameter=%.5f\n", pos_x, pos_y, parameter_x );
				pos_x_mean += pos_x;
				pos_y_mean += pos_y;
				parameter_x_mean += parameter_x;
			}
			pos_x_mean /= bpp;
			pos_y_mean /= bpp;
			parameter_x_mean /= bpp;

			printf( "alignment: cross correlation (fft): x0=%.5f, y0=%.5f, maximum=%.5f\n",
				pos_x_mean, pos_y_mean, parameter_x_mean );

			*quality = (gdouble)parameter_x_mean;
			*move_x = (gint)ROUND( pos_x_mean );
			*move_y = (gint)ROUND( pos_y_mean );

			break;

		default:
			break;
	}

	free( data );
}

gint get_dynamic_range( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gdouble *dynamic_range )
{
// 	gdouble data_mean = 0.;
// 	gint i, j;
// 	for ( i=0; i<Y; i++ )
// 	{
// 		for ( j=0; j<X; j++ )
// 		{
// 			data_mean += data[bpp*(i*X+j)+bit];
// 		}
// 	}
// 	data_mean /= X*Y;

	return 0;
}

gint get_fidelity( const guchar *data, const gint bit, const gint bpp, const gint X, const gint Y, gdouble *fidelity )
{
	return 0;
}

void get_quality( gint32 layer, gint quality_method, gint sel_pos_x, gint sel_pos_y, gint sel_width, gint sel_height,
	gdouble *quality )
{
	if ( quality_method == NONE || quality_method == ALIGNMENT_RELATED ) return;

	if ( sel_pos_x < 0 ) sel_pos_x = 0;
	if ( sel_pos_y < 0 ) sel_pos_y = 0;
	if ( sel_pos_x + sel_width > gimp_drawable_width( layer ) )  sel_width = gimp_drawable_width( layer ) - sel_pos_x;
	if ( sel_pos_y + sel_height > gimp_drawable_height( layer ) ) sel_height = gimp_drawable_height( layer ) - sel_pos_y;

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( layer ) )
	{
		case GIMP_RGB_IMAGE:
		case GIMP_RGBA_IMAGE:
			image_type = GIMP_RGB;
			layer_type = GIMP_RGB_IMAGE;
			bpp = 3;
			break;
		case GIMP_GRAY_IMAGE:
		case GIMP_GRAYA_IMAGE:
			image_type = GIMP_GRAY;
			layer_type = GIMP_GRAY_IMAGE;
			bpp = 1;
			break;
		case GIMP_INDEXED_IMAGE:
		case GIMP_INDEXEDA_IMAGE:
			return; /* Not allowed */
			break;
		default:
			g_message( _("Invalid color type!") );
			return;
	}

	guchar *data;
	data = malloc( gimp_drawable_bpp( layer )*sel_width*sel_height*sizeof(guchar) );

	GimpPixelRgn region_source;
	gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layer ),
		0, 0, gimp_drawable_width( layer ), gimp_drawable_height( layer ), FALSE, FALSE );

	gimp_pixel_rgn_get_rect( &region_source, data, sel_pos_x, sel_pos_y, sel_width, sel_height );

	gdouble parameter, parameter_mean;
	int i;

	switch( quality_method )
	{
		case CONTRAST:
			parameter_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				parameter = 0.;
				get_dynamic_range( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, &parameter );
				parameter_mean += parameter;
			}
			parameter_mean /= bpp;
			break;
		case FIDELITY:
			parameter_mean = 0.;
			for ( i=0; i<bpp; i++ )
			{
				parameter = 0.;
				get_fidelity( data, i, gimp_drawable_bpp( layer ), sel_width, sel_height, &parameter );
				parameter_mean += parameter;
			}
			parameter_mean /= bpp;
			break;
		default:
			parameter_mean = 0.;
			break;
	}

	*quality = parameter_mean;

	free( data );
}

/*
align layers
*/
static void align_layers( gint32 image_id )
{
	gimp_image_undo_group_start( image_id );

	gint32 layers_number, number;
	gint32 *layers;

	layers = gimp_image_get_layers( image_id, &layers_number );

	if ( layers_number < 2 )
	{
		gimp_image_undo_group_end( image_id );

		gimp_displays_flush();
		return;
	}

	gint layer_move_x[layers_number];
	gint layer_move_y[layers_number];
	gdouble layer_quality[layers_number];
	gint move_x = 0, move_y = 0;
	gdouble quality = 0.;

	gimp_progress_init( _("Align Layers...") );
	gint sel_pos_x, sel_pos_y, sel_width, sel_height;

	gimp_image_resize_to_layers( image_id );
	for( number = 0; number < layers_number; number++ )
	{
		gimp_layer_resize_to_image_size( layers[number] );
	}

/* FIXME: first visible layer or selected layer? */
	gint32 active_layer = /*gimp_image_get_active_layer( image_id );
	if ( active_layer == -1 ) active_layer =*/ 0;

	if ( parameters.scale_percent > 100 )
	{
		gint new_size_x = (double)(parameters.scale_percent*gimp_drawable_width( layers[active_layer] ))/100. + 0.5;
		gint new_size_y = (double)(parameters.scale_percent*gimp_drawable_height( layers[active_layer] ))/100. + 0.5;

		gimp_image_scale( image_id, new_size_x, new_size_y ); /* This also scales the selection! */
	}

	if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
	{
		sel_width = sel_width - sel_pos_x;
		sel_height = sel_height - sel_pos_y;
	}
	else
	{
		sel_pos_x = 0;
		sel_pos_y = 0;
		sel_width = gimp_drawable_width( layers[active_layer] );
		sel_height = gimp_drawable_height( layers[active_layer] );
	}

	guchar *reference_data = NULL;
	if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
	{
		if ( sel_width%2 ) sel_width--;
		if ( sel_height%2 ) sel_height--;

		reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );

		GimpPixelRgn region_source;
		gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
			0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );

		gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
	}

	for( number=layers_number - 1; number >= 0; number-- )
	{
		if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
		{
			if ( number == active_layer && parameters.alignment_method == CROSS_CORRELATION ) /* Don't skip this for CROSS_CORRELATION_FFT ! */
			{
				move_x = 0;
				move_y = 0;
				quality = 1.;
			}
			else
			{
				get_center( layers[number], parameters.alignment_method, sel_pos_x, sel_pos_y, sel_width, sel_height,
					parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
					&move_x, &move_y, &quality );
			}

			layer_move_x[number] = move_x;
			layer_move_y[number] = move_y;
			layer_quality[number] = quality;
		}
		else
		{
			layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
			layer_move_y[number] = 0;
			layer_quality[number] = 0;
		}

		gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
	}
	free( reference_data );

	gimp_progress_init( _("Moving layers by estimated values...") );
	char buffer[1024];
/* Calculate the "quality" always in the same way. Contrast? */
	sprintf( buffer, _("Quality: %.5f"), layer_quality[0] );
	gimp_drawable_set_name( layers[0], buffer );
	for( number=layers_number - 1; number > 0; number-- )
	{
		if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
		{
			if ( parameters.invisible_remove )
			{
				gimp_image_remove_layer( image_id, layers[number] );
			}
		}
		else
		{
			/* Sanity check of calculated move values */
			if ( abs( layer_move_x[0]-layer_move_x[number] ) <= sel_width*parameters.cross_correlation/100 &&
				abs( layer_move_y[0]-layer_move_y[number] ) <= sel_height*parameters.cross_correlation/100 )
			{
				gimp_layer_translate( layers[number], layer_move_x[0]-layer_move_x[number], layer_move_y[0]-layer_move_y[number] );
				sprintf( buffer, _("Quality: %.5f"), layer_quality[number] );
				gimp_drawable_set_name( layers[number], buffer );
			}
			else
			{
				g_message( _("Requested invalid move values. Not moving layer!") );
			}
		}

		gimp_progress_update( ((gdouble)layers_number-number) / (layers_number-1) );
	}

	gimp_image_resize_to_layers( image_id ); /* Resized and moved layers require a bigger image */

	/* Position of selection might have changed during resize */
	gint center_x = sel_pos_x + sel_width/2;
	gint center_y = sel_pos_y + sel_height/2;

	gint x_min = 0;
	gint y_min = 0;
	gint x_max = gimp_image_width( image_id );
	gint y_max = gimp_image_height( image_id );
	gint x_off, y_off;

	if ( parameters.trim_image && !parameters.two_star )
	{
		for( number=layers_number - 1; number >= 0; number-- )
		{
			if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
			{
				gimp_drawable_offsets( layers[number], &x_off, &y_off );
				if ( x_min < x_off ) x_min = x_off;
				if ( y_min < y_off ) y_min = y_off;
				if ( x_max > x_off + gimp_drawable_width( layers[number] ) ) x_max = x_off + gimp_drawable_width( layers[number] );
				if ( y_max > y_off + gimp_drawable_height( layers[number] ) ) y_max = y_off + gimp_drawable_height( layers[number] );
			}
		}

		gimp_image_crop( image_id, x_max - x_min, y_max - y_min, x_min, y_min );
	}

	if ( parameters.two_star )
	{
		gimp_progress_init( _("Align Layers...") );

		GtkWidget *dialog_widget;
		dialog_widget = gtk_message_dialog_new( 0,
			GTK_DIALOG_DESTROY_WITH_PARENT,
			GTK_MESSAGE_INFO,
			GTK_BUTTONS_OK_CANCEL,
			_("Please select another area for de-rotating layers.") );
		gimp_window_set_transient( GTK_WINDOW( dialog_widget ) );

		gint result = gtk_dialog_run( GTK_DIALOG( dialog_widget ) );
		gtk_widget_destroy( dialog_widget );
		if ( result == GTK_RESPONSE_OK )
		{
			active_layer = /*gimp_image_get_active_layer( image_id );
			if ( active_layer == -1 ) active_layer =*/ 0;

			if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
			{
				sel_width = sel_width - sel_pos_x;
				sel_height = sel_height - sel_pos_y;
			}
			else
			{
				sel_pos_x = 0;
				sel_pos_y = 0;
				sel_width = gimp_drawable_width( layers[active_layer] );
				sel_height = gimp_drawable_height( layers[active_layer] );
			}

			if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
			{
				if ( sel_width%2 ) sel_width--;
				if ( sel_height%2 ) sel_height--;

				reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );

				GimpPixelRgn region_source;
				gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
					0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );

				gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
			}

			gint center_x_corrected, center_y_corrected, sel_pos_x_corrected, sel_pos_y_corrected;
			gint x_off_0, y_off_0;
			gimp_drawable_offsets( layers[0], &x_off_0, &y_off_0 );

			for( number=layers_number - 1; number >= 0; number-- )
			{
				gimp_drawable_offsets( layers[number], &x_off, &y_off );
				sel_pos_x_corrected = sel_pos_x + x_off_0 - x_off;
				sel_pos_y_corrected = sel_pos_y + y_off_0 - y_off;

				if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
				{
					if ( number == active_layer && ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT ) )
					{
						move_x = 0;
						move_y = 0;
						quality = 1.;
					}
					else
					{
						get_center( layers[number], parameters.alignment_method, sel_pos_x_corrected, sel_pos_y_corrected, sel_width, sel_height,
							parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
							&move_x, &move_y, &quality );
					}

					layer_move_x[number] = move_x;
					layer_move_y[number] = move_y;
					layer_quality[number] = quality;
				}
				else
				{
					layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
					layer_move_y[number] = 0;
					layer_quality[number] = 0;
				}

				gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
			}
			free( reference_data );

			gimp_progress_init( _("Rotating layers by estimated values...") );
			gint32 selection = gimp_selection_save( image_id );
			gimp_selection_none( image_id );

			for( number=layers_number - 1; number > 0; number-- )
			{
				if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
				{
					if ( parameters.invisible_remove )
					{
						gimp_image_remove_layer( image_id, layers[number] );
					}
				}
				else
				{
					gimp_drawable_offsets( layers[number], &x_off, &y_off );
					center_x_corrected = center_x + x_off_0 - x_off;
					center_y_corrected = center_y + y_off_0 - y_off;

					gdouble distance_p0_pnumber =
						sqrt( (gdouble)((layer_move_x[0]-layer_move_x[number])*(layer_move_x[0]-layer_move_x[number]) +
							(layer_move_y[0]-layer_move_y[number])*(layer_move_y[0]-layer_move_y[number]) ) );
					gdouble distance_p0_center =
						sqrt( (gdouble)( (center_x_corrected-sel_pos_x-layer_move_x[0])*(center_x_corrected-sel_pos_x-layer_move_x[0]) +
							(center_y_corrected-sel_pos_y-layer_move_y[0])*(center_y_corrected-sel_pos_y-layer_move_y[0]) ) );

					gdouble angle = atan2( distance_p0_pnumber, distance_p0_center );

					/* (point0_x-center_x)*(pointnumber_y-point0_y) - (point0_y-center_y)*(pointnumber_x-point0_x) */
					gdouble crossproduct_3rd_component =
						(layer_move_x[0]+sel_pos_x-center_x_corrected)*(layer_move_y[number]-layer_move_y[0])-
						(layer_move_y[0]+sel_pos_y-center_y_corrected)*(layer_move_x[number]-layer_move_x[0]);

					if ( crossproduct_3rd_component > 0 ) angle = -angle;

					printf( "alignment: rotating by angle=%.5f around center x=%d, y=%d...\n", angle, center_x_corrected, center_y_corrected );

					gimp_drawable_transform_rotate( layers[number], angle, FALSE /* Auto-center */,
						center_x_corrected, center_y_corrected, GIMP_TRANSFORM_FORWARD,
						GIMP_INTERPOLATION_CUBIC, FALSE /* Supersample */, 3 /* recursion */,
						( parameters.trim_image ? GIMP_TRANSFORM_RESIZE_CROP : GIMP_TRANSFORM_RESIZE_ADJUST )/* clip results */ );
				}

				gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
			}
			gimp_selection_load( selection );

				active_layer = /*gimp_image_get_active_layer( image_id );
				if ( active_layer == -1 ) active_layer =*/ 0;

				if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
				{
					sel_width = sel_width - sel_pos_x;
					sel_height = sel_height - sel_pos_y;
				}
				else
				{
					sel_pos_x = 0;
					sel_pos_y = 0;
					sel_width = gimp_drawable_width( layers[active_layer] );
					sel_height = gimp_drawable_height( layers[active_layer] );
				}

				if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
				{
					if ( sel_width%2 ) sel_width--;
					if ( sel_height%2 ) sel_height--;

					reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );

					GimpPixelRgn region_source;
					gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
						0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );

					gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
				}

				gimp_drawable_offsets( layers[0], &x_off_0, &y_off_0 );

				for( number=layers_number - 1; number >= 0; number-- )
				{
					gimp_drawable_offsets( layers[number], &x_off, &y_off );
					sel_pos_x_corrected = sel_pos_x + x_off_0 - x_off;
					sel_pos_y_corrected = sel_pos_y + y_off_0 - y_off;

					if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
					{
						if ( number == active_layer && ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT ) )
						{
							move_x = 0;
							move_y = 0;
							quality = 1.;
						}
						else
						{
							get_center( layers[number], parameters.alignment_method, sel_pos_x_corrected, sel_pos_y_corrected, sel_width, sel_height,
								parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
								&move_x, &move_y, &quality );
						}

						layer_move_x[number] = move_x;
						layer_move_y[number] = move_y;
						layer_quality[number] = quality;
					}
					else
					{
						layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
						layer_move_y[number] = 0;
						layer_quality[number] = 0;
					}

					gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
				}
				free( reference_data );

			gimp_progress_init( _("Moving layers by estimated values...") );
			sprintf( buffer, _("Quality: %.5f"), layer_quality[0] );
			gimp_drawable_set_name( layers[0], buffer );
			for( number=layers_number - 1; number > 0; number-- )
			{
				if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
				{
					if ( parameters.invisible_remove )
					{
						gimp_image_remove_layer( image_id, layers[number] );
					}
				}
				else
				{
					/* Sanity check of calculated move values */
					if ( abs( layer_move_x[0]-layer_move_x[number] ) <= sel_width*parameters.cross_correlation/100 &&
						abs( layer_move_y[0]-layer_move_y[number] ) <= sel_height*parameters.cross_correlation/100 )
					{
						gimp_layer_translate( layers[number], layer_move_x[0]-layer_move_x[number], layer_move_y[0]-layer_move_y[number] );
						sprintf( buffer, _("Quality: %.5f"), layer_quality[number] );
						gimp_drawable_set_name( layers[number], buffer );
					}
					else
					{
						g_message( _("Requested invalid move values. Not moving layer!") );
					}
				}

				gimp_progress_update( ((gdouble)layers_number-number) / (layers_number-1) );
			}

			if ( parameters.trim_image )
			{
				x_min = 0;
				y_min = 0;
				x_max = gimp_image_width( image_id );
				y_max = gimp_image_height( image_id );
				for( number=layers_number - 1; number >= 0; number-- )
				{
					if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
					{
						gimp_drawable_offsets( layers[number], &x_off, &y_off );
						if ( x_min < x_off ) x_min = x_off;
						if ( y_min < y_off ) y_min = y_off;
						if ( x_max > x_off + gimp_drawable_width( layers[number] ) ) x_max = x_off + gimp_drawable_width( layers[number] );
						if ( y_max > y_off + gimp_drawable_height( layers[number] ) ) y_max = y_off + gimp_drawable_height( layers[number] );
					}
				}
				gimp_image_crop( image_id, x_max - x_min, y_max - y_min, x_min, y_min );
			}
		}
	}

	gimp_image_undo_group_end( image_id );

	gimp_displays_flush();
}


/*
GUI
*/

GtkObject *percent_adj;
GtkWidget *visible_layers;
GtkWidget *invisible_remove_layers;
GtkWidget *trim_image;
GtkWidget *subpixel_resolution;
GtkWidget *two_star;

static void visible_only_toggled( GtkWidget *check )
{
	parameters.visible_only = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void invisible_remove_toggled( GtkWidget *check )
{
	parameters.invisible_remove = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void trim_image_toggled( GtkWidget *check )
{
	parameters.trim_image = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void two_star_toggled( GtkWidget *check )
{
	parameters.two_star = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void subpixel_resolution_toggled( GtkWidget *check )
{
	parameters.subpixel_resolution = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

/*
main dialog
*/
static gint dialog()
{
	GtkWidget *dlg;
	GtkWidget *main_vbox;
	GtkWidget *frame;
	GtkWidget *table;
	GtkObject *adj;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( _("Align Layers"), "astro_align_layers", NULL, 0,
		gimp_standard_help_func, PLUG_IN_NAME,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		GTK_STOCK_OK, GTK_RESPONSE_OK,
		NULL);

	main_vbox = gtk_vbox_new( FALSE, 12 );
	gtk_container_set_border_width( GTK_CONTAINER( main_vbox ), 12 );
	gtk_container_add( GTK_CONTAINER( GTK_DIALOG( dlg )->vbox ), main_vbox );

	/* Generic settings */
	frame = gimp_frame_new( _("General Settings") );
	gtk_box_pack_start( GTK_BOX( main_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 5, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

		/* Scale before processing */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 0,
		_("Scale before processing (%):"), 180, 75,
		parameters.scale_percent, 100, 400, 1, 10, 0,
		TRUE, 0, 0, _("Scale layers before processing"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.scale_percent );

		/* visible layers ? */
	visible_layers = gtk_check_button_new_with_label( _("Use visible layers only") );
	gtk_table_attach_defaults( GTK_TABLE( table ), visible_layers, 0, 2, 1, 2 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( visible_layers ), parameters.visible_only );
	g_signal_connect( visible_layers, "toggled", G_CALLBACK( visible_only_toggled ), NULL );
	gtk_widget_show( visible_layers );

		/* remove invisible layers ? */
	invisible_remove_layers = gtk_check_button_new_with_label( _("Remove invisible layers") );
	gtk_table_attach_defaults( GTK_TABLE( table ), invisible_remove_layers, 0, 2, 2, 3 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( invisible_remove_layers ), parameters.invisible_remove );
	g_signal_connect( invisible_remove_layers, "toggled", G_CALLBACK( invisible_remove_toggled ), NULL );
	gtk_widget_show( invisible_remove_layers );

		/* trim image ? */
	trim_image = gtk_check_button_new_with_label( _("Trim image to overlap area") );
	gtk_table_attach_defaults( GTK_TABLE( table ), trim_image, 0, 2, 3, 4 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( trim_image ), parameters.trim_image );
	g_signal_connect( trim_image, "toggled", G_CALLBACK( trim_image_toggled ), NULL );
	gtk_widget_show( trim_image );

		/* two star alignment ? */
	two_star = gtk_check_button_new_with_label( _("Two star alignment") );
	gtk_table_attach_defaults( GTK_TABLE( table ), two_star, 0, 2, 4, 5 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( two_star ), parameters.two_star );
	g_signal_connect( two_star, "toggled", G_CALLBACK( two_star_toggled ), NULL );
	gtk_widget_show( two_star );

	/* alignment methods + settings */
	frame = gimp_frame_new( _("Alignment Method") );
	gtk_box_pack_start( GTK_BOX( main_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 2, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	GtkWidget *alignment_method = gimp_option_menu_new2( FALSE, G_CALLBACK( gimp_menu_item_update ),
		&parameters.alignment_method, (gpointer) parameters.alignment_method,
		_("Gauss fit"), (gpointer) 0, NULL,
		_("Center of brightness"), (gpointer) 1, NULL,
		_("Cross correlation"), (gpointer) 2, NULL,
		_("Cross correlation (FFT)"), (gpointer) 3, NULL,
		NULL );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Center of brightness mean value's radius (px):"), 180, 75,
		parameters.search_radius, 1, 100, 1, 10, 0,
		TRUE, 0, 0, _("Search radius, limited by the layer's border"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.search_radius );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Cross correlation search area (%):"), 180, 75,
		parameters.cross_correlation, 100, 1000, 1, 10, 0,
		TRUE, 0, 0, _("Cross correlation search area, measured in % of selected area width/height and limited by the layer's border"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.cross_correlation );

	gtk_table_attach( GTK_TABLE( table ), alignment_method, 0, 3, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( alignment_method );

	subpixel_resolution = gtk_check_button_new_with_label( _("Subpixel resolution when using cross correlation") );
	gtk_table_attach_defaults( GTK_TABLE( table ), subpixel_resolution, 0, 3, 3, 4 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( subpixel_resolution ), parameters.subpixel_resolution );
	g_signal_connect( subpixel_resolution, "toggled", G_CALLBACK( subpixel_resolution_toggled ), NULL );
	gtk_widget_show( subpixel_resolution );

	gtk_widget_show( main_vbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
