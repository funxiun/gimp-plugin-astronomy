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
gimp star_rounding plug-in
(C) Georg Hennig <georg.hennig@web.de>
Performs operations on each pixel, resulting in rounder stars, if they are
stretched.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-round-stars"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

enum REPLACE_METHOD
{
	MINIMUM_2 = 0,
	MINIMUM_5,
	MINIMUM_1_MEAN_4,
	MAXIMUM_2,
	MAXIMUM_5,
	MAXIMUM_1_MEAN_4
};

const gint RADIUS_MAXIMUM = 5;

/*
parameters
*/
typedef struct
{
	gdouble radial_length;
	gdouble radial_percentage;
	gdouble radial_x3;
	gdouble radial_x2;
	gdouble radial_x1;

	gdouble linear_length;
	gdouble linear_angle;
	gdouble linear_percentage;

	gint32 replace_method;
	gboolean rgb;
	gint32 threshold;

	gboolean show_preview;
} tparameter;


static tparameter parameters =
{
	2.,
	0.,
	0.5,
	2.,
	1.,
	2.,
	90.,
	100.,
	MINIMUM_2,
	TRUE,
	10,
	TRUE
};

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
star rounding
*/
static void star_rounding( GimpPreview *preview, GimpDrawable *drawable );

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
		{ GIMP_PDB_FLOAT, "radial_length", "Maximal radial length of blurred stars" },
		{ GIMP_PDB_FLOAT, "radial_percentage", "Replace pixel value by percentage weighted mean value" },
		{ GIMP_PDB_FLOAT, "radial_x3", "Cubic coefficient for radial increase of the length of blurred stars" },
		{ GIMP_PDB_FLOAT, "radial_x2", "Quadratic coefficient for radial increase of the length of blurred stars" },
		{ GIMP_PDB_FLOAT, "radial_x1", "Linear coefficient for radial increase of the length of blurred stars" },
		{ GIMP_PDB_FLOAT, "linear_length", "Length of blurred stars" },
		{ GIMP_PDB_FLOAT, "linear_angle", "Angle of blurred stars" },
		{ GIMP_PDB_FLOAT, "linear_percentage", "Replace pixel value by percentage weighted mean value" },
		{ GIMP_PDB_INT32, "replace_method", "Replace method to apply to each pixel" },
		{ GIMP_PDB_INT32, "rgb", "Treat colors independently, instead of brightness" },
		{ GIMP_PDB_INT32, "threshold", "Upper/lower limit for pixel values to be changed" },
		{ GIMP_PDB_INT32, "show_preview", "Toggle preview" }
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
		_("Rounds longish stars"),
		_("This plug-in applies different calculations to each pixel, resulting in rounded stars."),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Round Stars"),
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

			if ( !dialog( param[1].data.d_image, gimp_drawable_get( param[2].data.d_drawable ) ) ) return;

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 15 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				parameters.radial_length = param[3].data.d_float;
				parameters.radial_percentage = param[4].data.d_float;
				parameters.radial_x3 = param[5].data.d_float;
				parameters.radial_x2 = param[6].data.d_float;
				parameters.radial_x1 = param[7].data.d_float;

				parameters.linear_length = param[8].data.d_float;
				parameters.linear_angle = param[9].data.d_float;
				parameters.linear_percentage = param[10].data.d_float;

				parameters.replace_method = param[11].data.d_int32;
				parameters.rgb = param[12].data.d_int32;
				parameters.threshold = param[13].data.d_int32;
				parameters.show_preview = param[14].data.d_int32;
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
		star_rounding( NULL, gimp_drawable_get( param[2].data.d_drawable ) );
	}

	values[0].data.d_status = status;
}

#define MIN5(a, b, c, d, e) MIN(a, MIN(b, MIN(c, MIN(d, e))))
#define MAX5(a, b, c, d, e) MAX(a, MAX(b, MAX(c, MAX(d, e))))

guint MIN5_INDEX( gdouble a, gdouble b, gdouble c, gdouble d, gdouble e )
{
	gdouble min = MIN5( a, b, c, d, e );
	if ( a == min ) return 0;
	if ( b == min ) return 1;
	if ( c == min ) return 2;
	if ( d == min ) return 3;
	return 4;
}

guint MAX5_INDEX( gdouble a, gdouble b, gdouble c, gdouble d, gdouble e )
{
	gdouble max = MAX5( a, b, c, d, e );
	if ( a == max ) return 0;
	if ( b == max ) return 1;
	if ( c == max ) return 2;
	if ( d == max ) return 3;
	return 4;
}

/*
Apply this filter to a given selection
*/
static void process_selection( GimpPixelRgn *rgn_in, GimpPixelRgn *rgn_out,
	gint32 x_start, gint32 y_start, gint32 x_size, gint32 y_size, gint bpp_effective,
	gdouble length, gdouble angle, gdouble percentage,
	gboolean radial, gdouble x_maximum, gdouble y_maximum, gdouble x3, gdouble x2, gdouble x1,
	gint32 replace_method, gboolean rgb, gint32 threshold,
	gboolean show_progress )
{
	guchar **rows_input;
	guchar *row_output;
	guchar *tmp_row = NULL;

	guchar color[bpp_effective][6];
	gdouble brightness[6];
	gint32 x_closest[4], y_closest[4];
	gdouble x_percentage, y_percentage;
	gdouble pixel_weight[4];
	gdouble x, y;
	gint32 closest_x, closest_y;

	guint index_tmp;
	gdouble radius_tmp;
	gdouble polynom0_tmp, polynom_tmp;
	gdouble brightness_tmp;

	gdouble angle_final, length_final;
	length_final = length;
	angle_final = -angle + 180.;

	gint32 i, j, k, l;

	/* gimp documentation:
		If however the plug-in accesses drawable pixel data row-by-row,
		it should set the tile cache large enough to hold the number of tiles
		per row. Double this size if your plug-in uses shadow tiles. */
	gimp_tile_cache_ntiles( 2 * ( rgn_in->w / gimp_tile_width() + 1 ) );

  rows_input = g_new( guchar *, 2 * RADIUS_MAXIMUM + 1 );
  for( i=-RADIUS_MAXIMUM; i<=RADIUS_MAXIMUM; i++ )
		rows_input[i + RADIUS_MAXIMUM] = g_new( guchar, x_size * rgn_in->bpp );
  row_output = g_new( guchar, x_size * rgn_out->bpp );

  for( i=-RADIUS_MAXIMUM; i<=RADIUS_MAXIMUM; i++ )
	{
		gimp_pixel_rgn_get_row( rgn_in, rows_input[RADIUS_MAXIMUM+i],
			x_start, y_start + CLAMP( i, 0, y_size - 1 ), x_size );
	}

	for( i=1; i<=y_size; i++ )
	{
		/* Process current array */
		for ( j=0; j<x_size; j++ )
		{
			if ( radial ) /* otherwise it's already set previously */
			{
				radius_tmp = sqrt( (j+x_start-x_maximum/2)*(j+x_start-x_maximum/2) +
					(i+y_start-y_maximum/2)*(i+y_start-y_maximum/2) );
				polynom_tmp = x3 * pow( radius_tmp, 3. ) +
					x2 * pow( radius_tmp, 2. ) +
					x1 * radius_tmp;
				polynom0_tmp = x3 * pow( sqrt( pow( x_maximum/2, 2. ) + pow( y_maximum/2, 2. ) ), 3. ) +
					x2 * pow( sqrt( pow( x_maximum/2, 2. ) + pow( y_maximum/2, 2. ) ), 2. ) +
					x1 * sqrt( pow( x_maximum/2, 2. ) + pow( y_maximum/2, 2. ) );

				length_final = fabs( length ) * polynom_tmp / polynom0_tmp;

				angle_final = 180. * atan2( j+x_start-x_maximum/2, i+y_start-y_maximum/2 ) / M_PI;
				if ( length < 0. ) angle_final += 180.;
			}

			/* Value of current pixel */
			for ( k=0; k<bpp_effective; k++ ) color[k][0] = rows_input[RADIUS_MAXIMUM][rgn_in->bpp*j+k];
			if ( rgn_in->bpp == 4 || rgn_in->bpp == 3 )
			{
				brightness[0] = 0.30 * color[0][0] + 0.59 * color[1][0] + 0.11 * color[2][0];
			}
			else
			{
				brightness[0] = color[0][0];
			}

			/* Find the 4 pixels closest to the vector <current pixel> + <length_final, angle_final> */
			x = j + length_final * sin( angle_final * M_PI / 180. );
			y = RADIUS_MAXIMUM + length_final * cos(  angle_final * M_PI / 180. );

			x_closest[0] = floor( x );
			y_closest[0] = floor( y );
			closest_x = x_closest[0];
			closest_y = y_closest[0];
			x_closest[1] = ceil( x );
			y_closest[1] = floor( y );
			if ( sqrt( (x-x_closest[1])*(x-x_closest[1]) + (y-y_closest[1])*(y-y_closest[1]) ) <
				sqrt( (x-closest_x)*(x-closest_x) + (y-closest_y)*(y-closest_y) ) )
			{
				closest_x = x_closest[1];
				closest_y = y_closest[1];
			}
			x_closest[2] = floor( x );
			y_closest[2] = ceil( y );
			if ( sqrt( (x-x_closest[2])*(x-x_closest[2]) + (y-y_closest[2])*(y-y_closest[2]) ) <
				sqrt( (x-closest_x)*(x-closest_x) + (y-closest_y)*(y-closest_y) ) )
			{
				closest_x = x_closest[2];
				closest_y = y_closest[2];
			}
			x_closest[3] = ceil( x );
			y_closest[3] = ceil( y );
			if ( sqrt( (x-x_closest[3])*(x-x_closest[3]) + (y-y_closest[3])*(y-y_closest[3]) ) <
				sqrt( (x-closest_x)*(x-closest_x) + (y-closest_y)*(y-closest_y) ) )
			{
				closest_x = x_closest[3];
				closest_y = y_closest[3];
			}

			x_percentage = fabs( x - floor( x ) ); /* x_percentage for 0, 2; 1-x_percentage for 1, 3 */
			y_percentage = fabs( y - floor( y ) ); /* y_percentage for 0, 1; 1-y_percentage for 2, 3 */
			pixel_weight[0] = 1. - x_percentage * y_percentage;
			pixel_weight[1] = 1. - ( 1. - x_percentage ) * y_percentage;
			pixel_weight[2] = 1. - x_percentage * ( 1. - y_percentage );
			pixel_weight[3] = 1. - ( 1. - x_percentage ) * ( 1. - y_percentage );

			if ( x_closest[0] >= 0 && x_closest[1] >= 0 &&
				x_closest[2] >= 0 && x_closest[3] >= 0 &&
				x_closest[0] < x_size && x_closest[1] < x_size &&
				x_closest[2] < x_size && x_closest[3] < x_size &&
				y_closest[0] >= 0 && y_closest[1] >= 0 &&
				y_closest[2] >= 0 && y_closest[3] >= 0 &&
				y_closest[0] < y_size && y_closest[1] < y_size &&
				y_closest[2] < y_size && y_closest[3] < y_size &&
				( ( brightness[0] > threshold &&
					( replace_method == MINIMUM_2 || replace_method == MINIMUM_5 ||
					replace_method == MINIMUM_1_MEAN_4 ) ) ||
					( brightness[0] < threshold &&
					( replace_method == MAXIMUM_2 || replace_method == MAXIMUM_5 ||
					replace_method == MAXIMUM_1_MEAN_4 ) ) ) )
			{
				for ( k=0; k<4; k++ )
				{
					for ( l=0; l<bpp_effective; l++ )
					{
						color[l][k+1] = rows_input[y_closest[k]][rgn_in->bpp*x_closest[k]+l];
					}

					if ( rgn_in->bpp == 4 || rgn_in->bpp == 3 )
					{
						brightness[k+1] = 0.30 * color[0][k+1] + 0.59 * color[1][k+1] + 0.11 * color[2][k+1];
					}
					else
					{
						brightness[k+1] = color[0][k+1];
					}
				}

				for ( l=0; l<bpp_effective; l++ )
				{
					color[l][5] = rows_input[closest_y][rgn_in->bpp*closest_x+l];
				}
				if ( rgn_in->bpp == 4 || rgn_in->bpp == 3 )
				{
					brightness[5] = 0.30 * color[0][5] + 0.59 * color[1][5] + 0.11 * color[2][5];
				}
				else
				{
					brightness[5] = color[0][5];
				}

				if ( rgb )
				{
					for ( k=0; k<bpp_effective; k++ )
					{
						switch ( replace_method )
						{
							case MINIMUM_2:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MIN( color[k][0], color[k][5] ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							case MINIMUM_5:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MIN5( color[k][0], color[k][1],
									color[k][2], color[k][3], color[k][4] ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							case MINIMUM_1_MEAN_4:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MIN( color[k][0], ROUND(
									pixel_weight[0] * color[k][1] + pixel_weight[1] * color[k][2] +
									pixel_weight[2] * color[k][3] + pixel_weight[3] * color[k][4] ) ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							case MAXIMUM_2:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MAX( color[k][0], color[k][5] ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							case MAXIMUM_5:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MAX5( color[k][0], color[k][1],
									color[k][2], color[k][3], color[k][4] ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							case MAXIMUM_1_MEAN_4:
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * MAX( color[k][0], ROUND(
									pixel_weight[0] * color[k][1] + pixel_weight[1] * color[k][2] +
									pixel_weight[2] * color[k][3] + pixel_weight[3] * color[k][4] ) ) ) / 100. +
									(gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
								break;
							default:
								row_output[rgn_out->bpp*j+k] = 0;
								break;
						}
					}
				}
				else /* brightness instead of rgb */
				{
					switch ( replace_method )
					{
						case MINIMUM_2:
							if ( brightness[0] > 0. && brightness[0] > brightness[5] )
								for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness[5] ) /
										( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else
								for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						case MINIMUM_5:
							index_tmp = MIN5_INDEX( brightness[0], brightness[1], brightness[2], brightness[3], brightness[4] );
							if ( brightness[0] > 0. )
								for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness[index_tmp] ) /
										( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						case MINIMUM_1_MEAN_4:
							brightness_tmp = pixel_weight[0] * brightness[1] + pixel_weight[1] * brightness[2] +
								pixel_weight[2] * brightness[3] + pixel_weight[3] * brightness[4];
							brightness_tmp = MIN( brightness_tmp, brightness[0] );
							if ( brightness[0] > 0. ) for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness_tmp ) /
									( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						case MAXIMUM_2:
							if ( brightness[0] > 0. && brightness[0] < brightness[5] )
								for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness[5] ) /
										( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else
								for ( k=0; k<bpp_effective; k++ )
									row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						case MAXIMUM_5:
							index_tmp = MAX5_INDEX( brightness[0], brightness[1], brightness[2], brightness[3], brightness[4] );
							if ( brightness[0] > 0. ) for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness[index_tmp] ) /
									( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						case MAXIMUM_1_MEAN_4:
							brightness_tmp = pixel_weight[0] * brightness[1] + pixel_weight[1] * brightness[2] +
								pixel_weight[2] * brightness[3] + pixel_weight[3] * brightness[4];
							brightness_tmp = MAX( brightness_tmp, brightness[0] );
							if ( brightness[0] > 0. ) for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = (gdouble)( percentage * color[k][0] * brightness_tmp ) /
									( brightness[0] * 100. ) + (gdouble)( ( 100. - percentage ) * color[k][0] ) / 100. ;
							else for ( k=0; k<bpp_effective; k++ )
								row_output[rgn_out->bpp*j+k] = color[k][0];
							break;
						default:
							row_output[rgn_out->bpp*j+k] = 0;
							break;
					}
				}
			}
			else
			{
				for ( k=0; k<bpp_effective; k++ )
					row_output[rgn_out->bpp*j+k] = rows_input[RADIUS_MAXIMUM][rgn_in->bpp*j+k];
			}

			if ( rgn_in->bpp > bpp_effective )
				row_output[rgn_out->bpp*j+rgn_out->bpp-1] = rows_input[RADIUS_MAXIMUM][rgn_in->bpp*j+rgn_in->bpp-1];
		}

		/* Write result to pixel region */
		gimp_pixel_rgn_set_row( rgn_out, row_output, x_start, i + y_start - 1, x_size );

		/* Next line is the new last line of the array */
		gimp_pixel_rgn_get_row( rgn_in, rows_input[0],
			x_start, MIN( i + RADIUS_MAXIMUM + y_start, y_start + y_size - 1 ), x_size );

		tmp_row = rows_input[0];

		for( j = 1; j < 2 * RADIUS_MAXIMUM + 1; j++ ) rows_input[j-1] = rows_input[j];
		rows_input[2*RADIUS_MAXIMUM] = tmp_row;

		if ( show_progress && i % 10 == 0 ) gimp_progress_update( (gdouble)i / (gdouble)y_size );
	}

	for( i=0; i<2*RADIUS_MAXIMUM+1; i++ ) g_free( rows_input[i] );
  g_free( rows_input );
  g_free( row_output );
}

/*
Apply this filter to the drawable (alpha channel is ignored)
*/
static void star_rounding( GimpPreview *preview, GimpDrawable *drawable )
{
	if ( !preview ) gimp_progress_init( _("Round stars...") );

	gint x_start, y_start, x_size, y_size;
	if ( preview )
	{
		gimp_preview_get_position( preview, &x_start, &y_start );
		gimp_preview_get_size( preview, &x_size, &y_size );
	}
	else
	{
		if ( gimp_drawable_mask_bounds( drawable->drawable_id, &x_start, &y_start, &x_size, &y_size ) )
		{
			x_size = x_size - x_start;
			y_size = y_size - y_start;
		}
		else
		{
			x_start = 0;
			y_start = 0;
			x_size = drawable->width;
			y_size = drawable->height;
		}
	}

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( drawable->drawable_id ) )
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

	GimpPixelRgn region_source;
	gimp_pixel_rgn_init( &region_source, drawable, x_start, y_start,
		x_size, y_size, FALSE, FALSE );

	GimpPixelRgn region_destination;
	gimp_pixel_rgn_init( &region_destination, drawable, x_start, y_start,
		x_size, y_size, preview == NULL, TRUE );

	process_selection( &region_source, &region_destination, x_start, y_start, x_size, y_size, bpp,
		parameters.radial_length, 0., parameters.radial_percentage,
		TRUE, drawable->width, drawable->height, parameters.radial_x3, parameters.radial_x2, parameters.radial_x1,
		parameters.replace_method, parameters.rgb, parameters.threshold, preview == NULL );

	gimp_pixel_rgn_init( &region_source, drawable, x_start, y_start,
		x_size, y_size, FALSE, TRUE ); /* read from shadow tiles now */

	process_selection( &region_source, &region_destination, x_start, y_start, x_size, y_size, bpp,
		parameters.linear_length, parameters.linear_angle, parameters.linear_percentage,
		FALSE, 0., 0., 0., 0., 0.,
		parameters.replace_method, parameters.rgb, parameters.threshold, preview == NULL );

	if ( preview )
	{
		gimp_drawable_preview_draw_region( GIMP_DRAWABLE_PREVIEW( preview ), &region_destination );
	}
	else
	{
		gimp_drawable_flush( drawable );
		gimp_drawable_merge_shadow( drawable->drawable_id, TRUE );
		gimp_drawable_update( drawable->drawable_id, x_start, y_start, x_size, y_size );

		gimp_displays_flush();
	}
}

/*
GUI
*/

static void spin_x3( GtkWidget *spin )
{
	parameters.radial_x3 = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_x2( GtkWidget *spin )
{
	parameters.radial_x2 = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_x1( GtkWidget *spin )
{
	parameters.radial_x1 = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void rgb_toggled( GtkWidget *check )
{
	parameters.rgb = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

/*
main dialog
*/
static gint dialog( gint32 image_id, GimpDrawable *drawable )
{
	GtkWidget *dlg;
	GtkWidget *main_vbox;
	GtkWidget *main_hbox;
	GtkWidget *sub_vbox;
	GtkWidget *frame;
	GtkWidget *table;
	GtkWidget *label;
	GtkWidget *spin;
	GtkObject *adj;
	GtkWidget *preview;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( _("Rounds Stars"), "astro_round_stars", NULL, 0,
		gimp_standard_help_func, PLUG_IN_NAME,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		GTK_STOCK_OK, GTK_RESPONSE_OK,
		NULL );

	main_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( main_vbox ), 8 );
	gtk_container_add( GTK_CONTAINER( GTK_DIALOG( dlg )->vbox ), main_vbox );

	main_hbox = gtk_hbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( main_hbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_vbox ), main_hbox, FALSE, FALSE, 0 );

	preview = gimp_drawable_preview_new( drawable, &parameters.show_preview );
	gtk_box_pack_start( GTK_BOX( main_hbox ), preview, FALSE, FALSE, 0 );
	gtk_widget_show( preview );

	g_signal_connect( preview, "invalidated", G_CALLBACK( star_rounding ), drawable );

	sub_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( sub_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), sub_vbox, FALSE, FALSE, 0 );


	/* Radial correction */
	frame = gimp_frame_new( _("Radial Correction") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 1 );
	gtk_widget_show( frame );

	table = gtk_table_new( 1, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

		/* Length */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 0,
		_("Maximum length:"), 185, 75,
		parameters.radial_length, -5., 5., 0.1, 1., 1,
		TRUE, 0, 0, _("Length of the blurred stars"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.radial_length );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

		/* Percentage */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Percentage:"), 185, 75,
		parameters.radial_percentage, 0., 100., 1., 10., 1,
		TRUE, 0, 0, _("Replace pixel value by <percentage> of new pixel value"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.radial_percentage );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );


	table = gtk_table_new( 7, 1, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), table, FALSE, FALSE, 0 );
	gtk_widget_show( table );

	label = gtk_label_new( "   " );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );

	spin = gtk_spin_button_new_with_range( 0., 5., 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.radial_x3 );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_x3 ), NULL );
	gtk_table_attach( GTK_TABLE( table ), spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("* x^3 + ") );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );

	spin = gtk_spin_button_new_with_range( 0., 5., 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.radial_x2 );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_x2 ), NULL );
	gtk_table_attach( GTK_TABLE( table ), spin, 3, 4, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("* x^2 + ") );
	gtk_table_attach( GTK_TABLE( table ), label, 4, 5, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );

	spin = gtk_spin_button_new_with_range( 0., 5., 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.radial_x1 );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_x1 ), NULL );
	gtk_table_attach( GTK_TABLE( table ), spin, 5, 6, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("* x") );
	gtk_table_attach( GTK_TABLE( table ), label, 6, 7, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );


	/* Linear correction */
	frame = gimp_frame_new( _("Linear Correction") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 1 );
	gtk_widget_show( frame );

	table = gtk_table_new( 1, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

		/* Length */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 0,
		_("Length:"), 185, 75,
		parameters.linear_length, 0., 5., 0.1, 1., 1,
		TRUE, 0, 0, _("Length of the blurred stars"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.linear_length );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

		/* Angle */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Angle:"), 185, 75,
		parameters.linear_angle, 0., 360., 1., 10., 1,
		TRUE, 0, 0, _("Angle of the blurred stars"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.linear_angle );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

		/* Percentage */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Percentage:"), 185, 75,
		parameters.linear_percentage, 0., 100., 1., 10., 1,
		TRUE, 0, 0, _("Replace pixel value by <percentage> of new pixel value"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.linear_percentage );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );


	/* Replace method */
	frame = gimp_frame_new( _("Replacement Method") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 2 );
	gtk_widget_show( frame );

	table = gtk_table_new( 2, 1, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	GtkWidget *replace_method = gimp_int_combo_box_new(
		_("Minimum of pixels 0 and 1"), MINIMUM_2,
		_("Minimum of pixels 0-4"), MINIMUM_5,
		_("Minimum of pixel 0 and mean value of pixels 1-4"), MINIMUM_1_MEAN_4,
		_("Maximum of pixels 0 and 1"), MAXIMUM_2,
		_("Maximum of pixels 0-4"), MAXIMUM_5,
		_("Maximum of pixel 0 and mean value of pixels 1-4"), MAXIMUM_1_MEAN_4,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( replace_method ), parameters.replace_method );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( replace_method ), parameters.replace_method,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.replace_method );
	g_signal_connect_swapped( replace_method, "changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	gtk_table_attach( GTK_TABLE( table ), replace_method, 0, 2, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( replace_method );

		/* RGB or brightness ? */
	GtkWidget *rgb_button = gtk_check_button_new_with_label( _("Treat colors seperately (instead of brightness)") );
	gtk_table_attach_defaults( GTK_TABLE( table ), rgb_button, 0, 2, 1, 2 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( rgb_button ), parameters.rgb );
	g_signal_connect( rgb_button, "toggled", G_CALLBACK( rgb_toggled ), NULL );
	g_signal_connect_swapped( rgb_button, "toggled", G_CALLBACK( gimp_preview_invalidate ), preview );
	gtk_widget_show( rgb_button );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Threshold:"), 185, 75,
		parameters.threshold, 0, 255, 1, 10, 0,
		TRUE, 0, 0, _("Replace pixel value only if above threshold (for minimum) or below threshowld (for maximum)"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.threshold );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	gimp_preview_invalidate( GIMP_PREVIEW( preview ) );

	gtk_widget_show( sub_vbox );
	gtk_widget_show( main_hbox );
	gtk_widget_show( main_vbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
