/***************************************************************************
 *   Copyright (C) 2006-2018 by Georg Hennig                               *
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
gimp background_gradient plug-in
(C) Georg Hennig <georg.hennig@web.de>
Creates a set of pivotal points containing background level,
and fitting a polynomial of order 4 to create an (artificial) flat field.
Might also useful for real flat fields.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-background-gradient"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

enum BACKGROUND_VALUE_METHOD
{
	ALL_VALUES = 0,
	DARKEST_VALUES,
	BRIGHTEST_VALUES
};

/*
parameters
*/
typedef struct
{
	gint32 box_width;
	gint32 box_height;
	gint32 sigma_iterations;
	gdouble sigma_clip;
	gint32 minimum_number_box;

	gint32 background_value_method;
	gint32 use_percentage;
} tparameter;


static tparameter parameters =
{
	100,
	100,
	1,
	1.5,
	200,
	ALL_VALUES,
	10
};

const int NUMBER_COEFFICIENTS = 25;

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
background gradient
*/
static void background_gradient( gint32 image_id );

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
		{ GIMP_PDB_INT32, "box_width", "Width of box for pivotal points" },
		{ GIMP_PDB_INT32, "box_height", "Height of box for pivotal points" },
		{ GIMP_PDB_INT32, "sigma_iterations", "Iterations of sigma clipping" },
		{ GIMP_PDB_FLOAT, "sigma_clip", "Sigma clip of values inside the box" },
		{ GIMP_PDB_INT32, "minimum_number_box", "Minimum number of not clipping points" },
		{ GIMP_PDB_INT32, "background_value_method", "Use 0=all,1=darkest,2=brightest values" },
		{ GIMP_PDB_INT32, "use_percentage", "Use % of values (only for brightest/darkest values)" }
	};

	/*  Initialize i18n support  */
	bindtextdomain( GETTEXT_PACKAGE, gimp_locale_directory() );
#ifdef HAVE_BIND_TEXTDOMAIN_CODESET
	bind_textdomain_codeset( GETTEXT_PACKAGE, "UTF-8" );
#endif
	textdomain( GETTEXT_PACKAGE );

	static GimpParamDef *return_vals  = NULL;
	static int nparams = sizeof( params )/sizeof( params[0] );
	static int nreturn_vals = 0;

	gimp_install_procedure( PLUG_IN_NAME,
		_("Create a layer containing the background for an artificial flat field"),
		_("This plug-in tries to fit a polynomial of order 4 through the background pixel values."),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Background Gradient"),
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

			gint32 drawable = param[2].data.d_int32;
			if ( !dialog( drawable ) ) return;

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 10 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				parameters.box_width = param[3].data.d_int32;
				parameters.box_height = param[4].data.d_int32;
				parameters.sigma_iterations = param[5].data.d_int32;
				parameters.sigma_clip = param[6].data.d_float;
				parameters.minimum_number_box = param[7].data.d_int32;
				parameters.background_value_method = param[8].data.d_int32;
				parameters.use_percentage = param[9].data.d_int32;
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
		background_gradient( param[1].data.d_image );
	}

	values[0].data.d_status = status;
}

/*
Default compare function, for use with qsort
*/
int compare( const int *a, const int *b )
{
	if ( *a > *b ) return 1;
	if ( *a < *b ) return -1;
	return 0;
}

/*
Returns median of numbers from val[start] to val[end]
double for median of even numbers of elements.
Useful if used for sigma median
*/

double median( guint *arr, guint start, guint end )
{
	guint number = end - start + 1;

	double ret = 0.;

	if ( number%2 != 0 ) ret = arr[start+(number-1)/2];
	else ret = (double)(arr[start+number/2-1] + arr[start+number/2])/2;

	return ret;
}

/*
Returns mean value of given data set
*/
double mean( guint *arr, guint start, guint end )
{
	guint number = end - start + 1;

	double ret = 0.;

	int i;
	for ( i=start; i<=end; i++ )
	{
		ret += arr[i];
	}

	return ret/number;
}

/*
Returns the standard deviation of given data set
*/
double sigma( guint *arr, guint start, guint end )
{
	guint number = end - start + 1;

	double sum_x = 0.;
	double sum_x_x = 0.;

	int i;
	for ( i=start; i<=end; i++ )
	{
		sum_x += arr[i];
		sum_x_x += arr[i]*arr[i];
	}

	return sqrt((double)(number*sum_x_x-sum_x*sum_x)/(number*(number-1)));
}

void fit_polynomial( const double *x, const double *y, const double *z, const double *err, int start, int n,
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
		gsl_matrix_set( X, i, 3, x[i+start]*x[i+start]*x[i+start] );
		gsl_matrix_set( X, i, 4, x[i+start]*x[i+start]*x[i+start]*x[i+start] );
		gsl_matrix_set( X, i, 5, y[i+start] );
		gsl_matrix_set( X, i, 6, y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 7, y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 8, y[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 9, x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 10, x[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 11, x[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 12, x[i+start]*y[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 13, x[i+start]*x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 14, x[i+start]*x[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 15, x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 16, x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 17, x[i+start]*x[i+start]*x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 18, x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 19, x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 20, x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 21, x[i+start]*x[i+start]*x[i+start]*x[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 22, x[i+start]*x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 23, x[i+start]*x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start] );
		gsl_matrix_set( X, i, 24, x[i+start]*x[i+start]*x[i+start]*x[i+start]*y[i+start]*y[i+start]*y[i+start]*y[i+start] );

		gsl_vector_set( Y, i, z[i+start] );
		gsl_vector_set( W, i, 1.0/(err[i+start]*err[i+start]) );
	}

	gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc( n, NUMBER_COEFFICIENTS );
	gsl_multifit_wlinear( X, W, Y, c, cov, &chisq, work );
	gsl_multifit_linear_free( work );

	printf(
		"background_gradient: f(x,y)=%g %+g*x %+g*x**2 %+g*x**3 %+g*x**4 %+g*y %+g*y**2 %+g*y**3 %+g*y**4 %+g*x*y %+g*x*y**2 %+g*x*y**3 %+g*x*y**4 %+g*x**2*y %+g*x**2*y**2 %+g*x**2*y**3 %+g*x**2*y**4 %+g*x**3*y %+g*x**3*y**2 %+g*x**3*y**3 %+g*x**3*y**4 %+g*x**4*y %+g*x**4*y**2 %+g*x**4*y**3 %+g*x**4*y**4\n",
		gsl_vector_get( c, 0 ), gsl_vector_get( c, 1 ), gsl_vector_get( c, 2 ), gsl_vector_get( c, 3 ),
		gsl_vector_get( c, 4 ), gsl_vector_get( c, 5 ), gsl_vector_get( c, 6 ), gsl_vector_get( c, 7 ),
		gsl_vector_get( c, 8 ), gsl_vector_get( c, 9 ), gsl_vector_get( c, 10 ), gsl_vector_get( c, 11 ),
		gsl_vector_get( c, 12 ), gsl_vector_get( c, 13 ), gsl_vector_get( c, 14 ), gsl_vector_get( c, 15 ),
		gsl_vector_get( c, 16 ), gsl_vector_get( c, 17 ), gsl_vector_get( c, 18 ), gsl_vector_get( c, 19 ),
		gsl_vector_get( c, 20 ), gsl_vector_get( c, 21 ), gsl_vector_get( c, 22 ), gsl_vector_get( c, 23 ),
		gsl_vector_get( c, 24 ) );

	printf ( "background_gradient: covariance matrix:\n" );
	printf ( "  [ %+.5e, %+.5e, %+.5e  \n", gsl_matrix_get( cov, 0, 0 ),
		gsl_matrix_get( cov, 0, 1 ), gsl_matrix_get( cov, 0, 2 ) );
	printf ( "    %+.5e, %+.5e, %+.5e  \n", gsl_matrix_get( cov, 1, 0 ),
		gsl_matrix_get( cov, 1, 1 ), gsl_matrix_get( cov, 1, 2 ) );
	printf ( "    %+.5e, %+.5e, %+.5e ]\n", gsl_matrix_get( cov, 2, 0 ),
		gsl_matrix_get( cov, 2, 1 ), gsl_matrix_get( cov, 2, 2 ) );
	printf ( "background_gradient: chisq = %g\n", chisq );

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

double polynomial_value( const double x, const double y, double *coeffs, int start_coeffs )
{
	return coeffs[start_coeffs+0] + coeffs[start_coeffs+1]*x + coeffs[start_coeffs+2]*x*x +
		coeffs[start_coeffs+3]*x*x*x + coeffs[start_coeffs+4]*x*x*x*x + coeffs[start_coeffs+5]*y +
		coeffs[start_coeffs+6]*y*y + coeffs[start_coeffs+7]*y*y*y + coeffs[start_coeffs+8]*y*y*y*y +
		coeffs[start_coeffs+9]*x*y + coeffs[start_coeffs+10]*x*y*y + coeffs[start_coeffs+11]*x*y*y*y +
		coeffs[start_coeffs+12]*x*y*y*y*y + coeffs[start_coeffs+13]*x*x*y + coeffs[start_coeffs+14]*x*x*y*y +
		coeffs[start_coeffs+15]*x*x*y*y*y + coeffs[start_coeffs+16]*x*x*y*y*y*y + coeffs[start_coeffs+17]*x*x*x*y +
		coeffs[start_coeffs+18]*x*x*x*y*y + coeffs[start_coeffs+19]*x*x*x*y*y*y + coeffs[start_coeffs+20]*x*x*x*y*y*y*y +
		coeffs[start_coeffs+21]*x*x*x*x*y + coeffs[start_coeffs+22]*x*x*x*x*y*y + coeffs[start_coeffs+23]*x*x*x*x*y*y*y +
		coeffs[start_coeffs+24]*x*x*x*x*y*y*y*y;
}

static void background_gradient( gint32 image_id )
{
	gimp_image_undo_group_start( image_id );

	gint32 layer_source = gimp_image_get_active_layer( image_id );
	if ( layer_source == -1 ) return;

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( layer_source ) )
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
			image_type = GIMP_RGB;
			layer_type = GIMP_RGB_IMAGE;
			bpp = 3;
			break;
	}

	gint32 layer_destination = gimp_layer_new( image_id, _("Background gradient"), gimp_drawable_width( layer_source ),
		gimp_drawable_height( layer_source ), layer_type, 100, GIMP_NORMAL_MODE );

	GimpPixelRgn region_source;
	gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layer_source ),
		0, 0, gimp_drawable_width( layer_source ), gimp_drawable_height( layer_source ), FALSE, FALSE );

	GimpPixelRgn region_destination;
	gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
		gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

	int i = 0, j = 0, k = 0, l = 0;
	guchar buffer_source[parameters.box_width*parameters.box_height*region_source.bpp];
	guint values[parameters.box_width*parameters.box_height];

	guint pivotal_points_x = (gimp_drawable_width( layer_source )/parameters.box_width
		+(gimp_drawable_width( layer_source )%parameters.box_width>0?1:0));
	guint pivotal_points_y = (gimp_drawable_height( layer_source )/parameters.box_height
		+(gimp_drawable_height( layer_source )%parameters.box_height>0?1:0));
	guint pivotal_points = pivotal_points_x*pivotal_points_y;
	double x[bpp*pivotal_points];
	double y[bpp*pivotal_points];
	double z[bpp*pivotal_points];
	double err[bpp*pivotal_points];

	guint real_pivotal_points[bpp];
	for ( i=0; i<bpp; i++ ) real_pivotal_points[i] = 0;

	guint x_start, width, y_start, height;

	gimp_progress_init( _("Calculating mean values for pivotal points...") );

	guint progress_skip = pivotal_points/50;
	if ( progress_skip == 0 ) progress_skip = 1;

	for ( i=0; i<pivotal_points; i++ )
	{
		x_start = (i%pivotal_points_x)*parameters.box_width;
		y_start = (i/pivotal_points_x)*parameters.box_height;
		width = (x_start+parameters.box_width-1<=gimp_drawable_width( layer_source ) ?
			parameters.box_width :
			gimp_drawable_width( layer_source )-x_start);
		height = (y_start+parameters.box_height-1<=gimp_drawable_height( layer_source ) ?
			parameters.box_height :
			gimp_drawable_height( layer_source )-y_start);

		gimp_pixel_rgn_get_rect( &region_source, buffer_source,
			x_start, y_start, width, height );

		for ( j=0; j<bpp; j++ )
		{
			guint offset = 0;
			for ( k=0; k<width; k++ )
			{
				for ( l=0; l<height; l++ )
				{
					if ( buffer_source[region_source.bpp*(l*width+k)+region_source.bpp-1] == 0 )
					{
						values[l*width+k] = 0;
						offset++;
					}
					else
					{
						values[l*width+k] = buffer_source[region_source.bpp*(l*width+k)+j];
					}
				}
			}

			guint values_number = width * height - offset;
			if ( parameters.background_value_method == DARKEST_VALUES || parameters.background_value_method == BRIGHTEST_VALUES )
				values_number *= (gdouble)(parameters.use_percentage) / 100;

			if ( values_number >  parameters.minimum_number_box )
			{
				qsort( values, width*height, sizeof(guint), (void *)compare );

				guint lower_limit = offset;
				guint upper_limit = width*height-1;

				if ( parameters.background_value_method == DARKEST_VALUES )
				{
					upper_limit = ( width * height - 1 - offset ) * (gdouble)(parameters.use_percentage) / 100;
				}
				else if ( parameters.background_value_method == BRIGHTEST_VALUES )
				{
					lower_limit = ( width * height - 1 - offset ) * ( 1. - (gdouble)(parameters.use_percentage) / 100 ) + offset;
				}

				for ( k=0; k<parameters.sigma_iterations; k++ )
				{
					double median_value = median( values, lower_limit, upper_limit );
					double sigma_value = sigma( values, lower_limit, upper_limit );

					guint lower_limit_value = CLAMP( ROUND( median_value - parameters.sigma_clip*sigma_value ), values[offset], 255 );
					guint upper_limit_value = CLAMP( ROUND( median_value + parameters.sigma_clip*sigma_value ), values[offset], 255 );

					for ( l=offset+lower_limit; l<=upper_limit; l++ )
					{
						if ( values[l] >= lower_limit_value )
						{
							lower_limit = l;
							break;
						}
					}
					for ( l=upper_limit; l>=offset+lower_limit; l-- )
					{
						if ( values[l] <= upper_limit_value )
						{
							upper_limit = l;
							break;
						}
					}
					if ( lower_limit > upper_limit ) upper_limit = lower_limit;
				}

				double mean_value = mean( values, lower_limit, upper_limit );

				if ( upper_limit-lower_limit+1 > parameters.minimum_number_box )
				{
					x[j*pivotal_points+real_pivotal_points[j]] = x_start + 0.5*width;
					y[j*pivotal_points+real_pivotal_points[j]] = y_start + 0.5*height;
					z[j*pivotal_points+real_pivotal_points[j]] = mean_value;
					err[j*pivotal_points+real_pivotal_points[j]] = 1.0; /* Fake error for now */
					real_pivotal_points[j] += 1;
				}
			}
		}

		if ( i%progress_skip == 0 ) gimp_progress_update( (double)(i+1)/pivotal_points );
	}

	double coeffs[bpp*NUMBER_COEFFICIENTS];
	double maximum[bpp];
	double value;

	for ( i=0; i<bpp; i++ )
	{
		if ( real_pivotal_points[i] < NUMBER_COEFFICIENTS )
		{
			g_message( _("Too few pivotal points found!") );

			gimp_drawable_delete( layer_destination );

			gimp_image_undo_group_end( image_id );

			gimp_displays_flush();

			return;
		}
	}

	gimp_progress_init( _("Fitting polynomial and calculating maximum pixel values...") );

	/* Alternative maximum determination:
	Local maxima in each line + absolute maxima on the border
	x^3 + a x^2 + b x + c = 0 */
	double max0, max1, max2;
	double a0;

	width = gimp_drawable_width( layer_destination );
	height = gimp_drawable_height( layer_destination );

	progress_skip = height/50;
	if ( progress_skip == 0 ) progress_skip = 1;
	for ( i=0; i<bpp; i++ )
	{
		fit_polynomial( x, y, z, err, i*pivotal_points, real_pivotal_points[i], coeffs, i*NUMBER_COEFFICIENTS );

		maximum[i] = 0;

		for ( k=0; k<width; k++ )
		{
			value = polynomial_value( k, 0, coeffs, i*NUMBER_COEFFICIENTS );
			if ( maximum[i] < value ) maximum[i] = value;
		}

		for ( k=1; k<height-1; k++ )
		{
			value = polynomial_value( 0, k, coeffs, i*NUMBER_COEFFICIENTS );
			if ( maximum[i] < value ) maximum[i] = value;

			max0 = 0.;
			max1 = 0.;
			max2 = 0.;

			a0 = 4*(coeffs[i*NUMBER_COEFFICIENTS+4]+coeffs[i*NUMBER_COEFFICIENTS+21]*k+coeffs[i*NUMBER_COEFFICIENTS+22]*k*k+
				coeffs[i*NUMBER_COEFFICIENTS+23]*k*k*k+coeffs[i*NUMBER_COEFFICIENTS+24]*k*k*k*k);
			gsl_poly_solve_cubic(
				3*(coeffs[i*NUMBER_COEFFICIENTS+3]+coeffs[i*NUMBER_COEFFICIENTS+17]*k+coeffs[i*NUMBER_COEFFICIENTS+18]*k*k+
					coeffs[i*NUMBER_COEFFICIENTS+19]*k*k*k+coeffs[i*NUMBER_COEFFICIENTS+20]*k*k*k*k)/a0,
				2*(coeffs[i*NUMBER_COEFFICIENTS+2]+coeffs[i*NUMBER_COEFFICIENTS+13]*k+coeffs[i*NUMBER_COEFFICIENTS+14]*k*k+
					coeffs[i*NUMBER_COEFFICIENTS+15]*k*k*k+coeffs[i*NUMBER_COEFFICIENTS+16]*k*k*k*k)/a0,
				(coeffs[i*NUMBER_COEFFICIENTS+1]+coeffs[i*NUMBER_COEFFICIENTS+9]*k+coeffs[i*NUMBER_COEFFICIENTS+10]*k*k+
					coeffs[i*NUMBER_COEFFICIENTS+11]*k*k*k+coeffs[i*NUMBER_COEFFICIENTS+12]*k*k*k*k)/a0,
				&max0, &max1, &max2 );

			if ( max0 < width-1 && max0 >= 0.5 )
			{
				value = polynomial_value( ceil( max0 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
				value = polynomial_value( floor( max0 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
			}

			if ( max1 < width-1 && max1 >= 0.5  )
			{
				value = polynomial_value( ceil( max1 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
				value = polynomial_value( floor( max1 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
			}

			if ( max2 < width-1 && max2 >= 0.5  )
			{
				value = polynomial_value( ceil( max2 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
				value = polynomial_value( floor( max2 ), k, coeffs, i*NUMBER_COEFFICIENTS );
				if ( maximum[i] < value ) maximum[i] = value;
			}

			value = polynomial_value( width-1, k, coeffs, i*NUMBER_COEFFICIENTS );
			if ( maximum[i] < value ) maximum[i] = value;

			if ( k%progress_skip == 0 ) gimp_progress_update( (double)(k+1+i*height)/(bpp*height) );
		}

		for ( k=0; k<width; k++ )
		{
			value = polynomial_value( k, height-1, coeffs, i*NUMBER_COEFFICIENTS );
			if ( maximum[i] < value ) maximum[i] = value;
		}
	}

	double max = 0.;
	for ( i=0; i<bpp; i++ ) if ( max < maximum[i] ) max = maximum[i];

	printf( "background_gradient: max: %f\n", max );

	gimp_progress_init( _("Writing normalized pixel values...") );

 	gimp_image_add_layer( image_id, layer_destination, 0 );

	guint progress = 0;
	progress_skip = 0;
	guchar *src, *s;
	gpointer pr;
	for( pr = gimp_pixel_rgns_register( 1, &region_destination ); pr != NULL; pr = gimp_pixel_rgns_process( pr ) )
	{
		src = region_destination.data;

		for ( i=0; i<region_destination.h; i++ )
		{
			s = src;

			for ( j=0; j<region_destination.w; j++ )
			{
				for ( k=0; k<bpp; k++ )
				{
					s[k] = CLAMP( max>0. ?
						ROUND( 255*polynomial_value( j+region_destination.x, i+region_destination.y, coeffs, k*NUMBER_COEFFICIENTS )/max ) :
						255, 0, 255 );
				}

				if ( region_destination.bpp > bpp ) s[region_destination.bpp-1] = 255;

				s += region_destination.bpp;
			}

			src += region_destination.rowstride;
		}

		progress += region_destination.w * region_destination.h;
		progress_skip++;
		if ( progress_skip%50 == 0 ) /* only update progress every xth tile to avoid slowdown */
			gimp_progress_update( (double)(progress) /
				( region_destination.drawable->width * region_destination.drawable->height ) );
	}

	gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
	gimp_drawable_merge_shadow( layer_destination, FALSE );
	gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );

	gimp_layer_set_mode( layer_destination, GIMP_LAYER_MODE_DIVIDE );

	gimp_image_undo_group_end( image_id );

	gimp_displays_flush();
}

/*
GUI
*/
GtkObject *minimum_adj;
GtkWidget *use_percentage_label;
GtkWidget *use_percentage;

static void spin_width( GtkWidget *spin )
{
	parameters.box_width = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );

	if ( parameters.box_width*parameters.box_height < parameters.minimum_number_box )
	{
		gtk_adjustment_set_value( GIMP_SCALE_ENTRY_SCALE_ADJ( minimum_adj ), parameters.box_width*parameters.box_height );
	}
}

static void spin_height( GtkWidget *spin )
{
	parameters.box_height = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_sigma( GtkWidget *spin )
{
	parameters.sigma_iterations = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_use_percentage( GtkWidget *spin )
{
	parameters.use_percentage = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void background_value_changed( GtkWidget *combo_box )
{
	gint active;
	gimp_int_combo_box_get_active( GIMP_INT_COMBO_BOX( combo_box ), &active );
	switch ( active )
	{
		case ALL_VALUES:
			gtk_widget_set_sensitive( use_percentage_label, FALSE );
			gtk_widget_set_sensitive( use_percentage, FALSE );
			break;
		case DARKEST_VALUES:
		case BRIGHTEST_VALUES:
		default:
			gtk_widget_set_sensitive( use_percentage_label, TRUE );
			gtk_widget_set_sensitive( use_percentage, TRUE );
			break;
	}
}

/*
main dialog
*/
static gint dialog( gint32 image_id, GimpDrawable *drawable )
{
	GtkWidget *dlg;
	GtkWidget *main_vbox;
	GtkWidget *table;
	GtkObject *adj;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( _("Background Gradient"), "astro_background_gradient", NULL, 0,
		gimp_standard_help_func, PLUG_IN_NAME,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		GTK_STOCK_OK, GTK_RESPONSE_OK,
		NULL );

	main_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( main_vbox ), 8 );
	gtk_container_add( GTK_CONTAINER( GTK_DIALOG( dlg )->vbox ), main_vbox );

	/* Box size */
	table = gtk_table_new( 2, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( main_vbox ), table );
	gtk_widget_show( table );

		/* Box width */
	GtkWidget *width_spin_label = gtk_label_new( _("Box width:") );
	GtkWidget *width_spin = gtk_spin_button_new_with_range( 10, 200, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( width_spin ), parameters.box_width );

	g_signal_connect( width_spin, "value_changed", G_CALLBACK( spin_width ), NULL );

	gtk_table_attach( GTK_TABLE( table ), width_spin_label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( width_spin_label );

	gtk_table_attach( GTK_TABLE( table ), width_spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( width_spin );

		/* Box height */
	GtkWidget *height_spin_label = gtk_label_new( _("Box height:") );
	GtkWidget *height_spin = gtk_spin_button_new_with_range( 10, 200, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( height_spin ), parameters.box_height );

	g_signal_connect( height_spin, "value_changed", G_CALLBACK( spin_height ), NULL );

	gtk_table_attach( GTK_TABLE( table ), height_spin_label, 2, 3, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( height_spin_label );

	gtk_table_attach( GTK_TABLE( table ), height_spin, 3, 4, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( height_spin );

	/* Sigma clip */
		/* Clip iterations */
	/* Box size */
	table = gtk_table_new( 1, 2, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( main_vbox ), table );
	gtk_widget_show( table );

	GtkWidget *sigma_spin_label = gtk_label_new( _("Iterations of sigma clipping:") );
	GtkWidget *sigma_spin = gtk_spin_button_new_with_range( 1, 5, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( sigma_spin ), parameters.sigma_iterations );

	g_signal_connect( sigma_spin, "value_changed", G_CALLBACK( spin_sigma ), NULL );

	gtk_table_attach( GTK_TABLE( table ), sigma_spin_label, 0, 1, 1, 2, GTK_FILL, 0, 0, 0);
	gtk_widget_show( sigma_spin_label );

	gtk_table_attach( GTK_TABLE( table ), sigma_spin, 1, 2, 1, 2, GTK_FILL, 0, 0, 0);
	gtk_widget_show( sigma_spin );

		/* x*Sigma clip */
	table = gtk_table_new( 2, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( main_vbox ), table );
	gtk_widget_show( table );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 0,
		_("Sigma clipping:"), 185, 75,
		parameters.sigma_clip, 0.5, 6.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Clip values inside the box by ... sigma"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.sigma_clip );

	/* Minimum number */
	minimum_adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Min number of non-clipping values:"), 185, 75,
		parameters.minimum_number_box, 10, 2000, 1, 5, 0,
		TRUE, 0, 0, _("Minimum number of values that do not clip inside the box"), NULL );
	g_signal_connect( minimum_adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.minimum_number_box );

	/* use all / darkest / brightest values */
	table = gtk_table_new( 1, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( main_vbox ), table );
	gtk_widget_show( table );


	GtkWidget *background_value = gimp_int_combo_box_new(
		_("All values"), ALL_VALUES,
		_("Darkest values"), DARKEST_VALUES,
		_("Brightest values"), BRIGHTEST_VALUES,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( background_value ), parameters.background_value_method );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( background_value ), parameters.background_value_method,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.background_value_method );
	g_signal_connect( background_value, "changed", G_CALLBACK( background_value_changed ), NULL );

	gtk_table_attach( GTK_TABLE( table ), background_value, 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( background_value );

	use_percentage_label = gtk_label_new( _("Percentage:") );
	use_percentage = gtk_spin_button_new_with_range( 1, 100, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( use_percentage ), parameters.use_percentage );

	g_signal_connect( use_percentage, "value_changed", G_CALLBACK( spin_use_percentage ), NULL );

	gtk_table_attach( GTK_TABLE( table ), use_percentage_label, 1, 2, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( use_percentage_label );

	gtk_table_attach( GTK_TABLE( table ), use_percentage, 2, 3, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( use_percentage );

	gint active;
	gimp_int_combo_box_get_active( GIMP_INT_COMBO_BOX( background_value ), &active );
	switch ( active )
	{
		case ALL_VALUES:
			gtk_widget_set_sensitive( use_percentage_label, FALSE );
			gtk_widget_set_sensitive( use_percentage, FALSE );
			break;
		case DARKEST_VALUES:
		case BRIGHTEST_VALUES:
		default:
			gtk_widget_set_sensitive( use_percentage_label, TRUE );
			gtk_widget_set_sensitive( use_percentage, TRUE );
			break;
	}

	gtk_widget_show( main_vbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
