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
gimp merge plug-in
(C) Georg Hennig <georg.hennig@web.de>
Merges all (or all visible) layers using arithmethic, geometric, median, arithmetic
of median +- sigma or sigma median (1 or 2 pass) algorithm.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if !defined(WIN32)
/* on GNU systems, nothing necessary */
#else
char *strndup(const char *s, size_t n);
#endif

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-merge-layers"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

enum MERGE_METHOD
{
	ARITHMETIC = 0,
	GEOMETRIC,
	MEDIAN,
	SIGMA_MEDIAN_ARITHMETIC,
	SIGMA_MEDIAN_1,
	SIGMA_MEDIAN_2
};

/*
parameters
*/
typedef struct
{
	gint32 merge_method;
	gdouble sigma_merge_method;

	gboolean visible_only;

	gdouble multiply_constant;
	gint32 add_constant;

	gboolean show_preview;
} tparameter;


static tparameter parameters =
{
	SIGMA_MEDIAN_ARITHMETIC,
	1.0,
	TRUE,
	1.0,
	0,
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
merge layers
*/
static void merge_layers();

static void merge_selection( gint32 *layers_visible, gint32 layers_number_visible, GimpPixelRgn *region_destination,
	gboolean visible_only, gint32 x_start, gint32 y_start, gint32 x_size, gint32 y_size,
	gint32 method, gint bpp_effective, gdouble sigma_merge_method, gdouble multiply_constant, gint32 add_constant );

int compare( const int *a, const int *b );
double median( gint *arr, guint start, guint end );

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

gint32 image_id;

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
		{ GIMP_PDB_INT32, "merge_method", "Layer Merge Method (0,1,2,3,4)" },
		{ GIMP_PDB_FLOAT, "sigma", "Clip to x sigma" },
		{ GIMP_PDB_INT32, "visible_only", "Process visible layers only" },
		{ GIMP_PDB_FLOAT, "multiply_constant", "Multiply values before merging" },
		{ GIMP_PDB_INT32, "add_constant", "Add value to values before merging" },
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
		_("Merge layers using different mean value algorithms"),
		_("This plug-in merges layers using different mean value algorithms. "),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Merge Layers"),
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

	image_id = param[1].data.d_image;

	gint32 layers_number, number;
	gint32 *layers;

	layers = gimp_image_get_layers( image_id, &layers_number );

	gint layer_width[layers_number];
	gint layer_height[layers_number];
	gint layer_x_offset[layers_number];
	gint layer_y_offset[layers_number];

	gimp_image_undo_disable( image_id );
	for( number = 0; number < layers_number; number++ )
	{
		layer_width[number] = gimp_drawable_width( layers[number] );
		layer_height[number] = gimp_drawable_height( layers[number] );
		gimp_drawable_offsets( layers[number], &layer_x_offset[number], &layer_y_offset[number] );
		gimp_layer_resize_to_image_size( layers[number] );
	}

	switch(run_mode)
	{
		case GIMP_RUN_INTERACTIVE:
			gimp_get_data( PLUG_IN_NAME, &parameters );

			if ( !dialog( param[1].data.d_image, gimp_drawable_get( param[2].data.d_drawable ) ) )
			{
				for( number = layers_number-1; number >= 0; number-- )
				{
					gimp_layer_resize( layers[number], layer_width[number], layer_height[number],
						-layer_x_offset[number], -layer_y_offset[number] );
				}
				gimp_image_undo_enable( image_id );
				return;
			}

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 9 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				image_id = param[1].data.d_image;

				parameters.merge_method = param[3].data.d_int32;
				parameters.sigma_merge_method = param[4].data.d_float;
				parameters.visible_only = param[5].data.d_int32;
				parameters.multiply_constant = param[6].data.d_float;
				parameters.add_constant = param[7].data.d_int32;
				parameters.show_preview = param[8].data.d_int32;
			}
			break;
		case GIMP_RUN_WITH_LAST_VALS:
			gimp_get_data( PLUG_IN_NAME, &parameters );
			image_id = param[1].data.d_image;

			break;
		default:
			status = GIMP_PDB_CALLING_ERROR;
			break;
	}

	if( status == GIMP_PDB_SUCCESS )
	{
		merge_layers();
	}

	for( number = layers_number-1; number >= 0; number-- )
	{
		gimp_layer_resize( layers[number], layer_width[number], layer_height[number],
			-layer_x_offset[number], -layer_y_offset[number] );
	}
	gimp_image_undo_enable( image_id );

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

double median( gint *arr, guint start, guint end )
{
	guint number = end - start + 1;

	double ret = 0;

	if ( number%2 != 0 ) ret = arr[start+(number-1)/2];
	else ret = (double)(arr[start+number/2-1] + arr[start+number/2])/2;

	return ret;
}

/*
Merges layer (visible or all) from image image_id to destination_layer, 
Source image to merge layers, destination layer (must exist and belong to an image),
*/
static void merge_selection( gint32 *layers_visible, gint32 layers_number_visible, GimpPixelRgn *region_destination,
	gboolean visible_only, gint32 x_start, gint32 y_start, gint32 x_size, gint32 y_size,
	gint32 method, gint bpp_effective, gdouble sigma_merge_method, gdouble multiply_constant, gint32 add_constant )
{
	gint32 number;

	/* initialize regions of layers to merge (source layers) */
	/* Last rgn will become the destination layer! */
	GimpPixelRgn **rgn;
	rgn = g_malloc( (layers_number_visible+1)*sizeof( GimpPixelRgn* ) );
	gint32 number_region = 0;

	for( number = layers_number_visible - 1; number >= 0; number-- )
	{
		rgn[number] = ( GimpPixelRgn * ) g_malloc( sizeof( GimpPixelRgn ) );
		gimp_pixel_rgn_init( rgn[number], gimp_drawable_get( layers_visible[number] ),
			x_start, y_start, x_size, y_size, FALSE, FALSE );
		number_region++;
	}

	rgn[layers_number_visible] = region_destination;

	guchar *src[layers_number_visible+1], *s[layers_number_visible+1];
	gint32 x, y;
	gpointer pr;
	gint i, j;
	gint visible_pixels;
	gint visible_number;
	double sum, sum_square;
	double product;
	double sigma, median_tmp;
	gint start, end;
	double alpha_sum;

	guint progress = 0;
	guint progress_counter = 0;
	/* loop through all source layers */
	for( pr = gimp_pixel_rgns_register2( layers_number_visible+1, rgn ); pr != NULL; pr = gimp_pixel_rgns_process( pr ) )
	{
		for ( number = 0; number <= layers_number_visible; number++ ) src[number] = rgn[number]->data;

		for ( y = 0; y < rgn[layers_number_visible]->h; y++ )
		{
			for ( number = 0; number <= layers_number_visible; number++ ) s[number] = src[number];

			for ( x = 0; x < rgn[layers_number_visible]->w; x++ )
			{
				/* The actual pixels are s[number=0,...,layers_number_visible-1][i=0,...,bpp-1] */
				/* The destination pixels s[layers_number_visible][i=0,...,bpp-1] */
				for ( i = 0; i < bpp_effective; i++ )
				{
					switch( parameters.merge_method )
					{
						case ARITHMETIC: /* Arithmetic mean value */
							/* value = (Sum of all values)/(Number of all values) */
							sum = 0.;
							alpha_sum = 0.;
							for ( number = 0; number < layers_number_visible; number++ )
							{
								sum += ( (double)(s[number][i])*multiply_constant + add_constant );
								/* If current layer has alpha channel, only include a weighted pixel value */
								if ( rgn[number]->bpp > bpp_effective ) alpha_sum += (double)(s[number][rgn[number]->bpp-1])/255;
								else alpha_sum += 1.;
							}
							if ( alpha_sum < 1e-10 ) alpha_sum = 1e-10;
							/* CLAMP(x, low, high)  (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x))) */
							s[layers_number_visible][i] =
								(guchar) CLAMP( ROUND( sum/alpha_sum ), 0, 255 );
							break;
						case GEOMETRIC: /* Geometric mean value */
							/* value = (Number of all values)th root of (Product of all values) */
							product = 1.;
							alpha_sum = 0.;

							for ( number = 0; number < layers_number_visible; number++ )
							{
								if ( rgn[number]->bpp > bpp_effective )
								{
									alpha_sum += (double)(s[number][rgn[number]->bpp-1])/255;
								}
								else alpha_sum += 1.;
							}
							if ( alpha_sum < 1e-10 ) alpha_sum = 1e-10;

							for ( number = 0; number < layers_number_visible; number++ )
							{
							/* Calculate root now to avoid INF (for example for 200 layers this might happen soon) */
								if ( (double)(s[number][i])*multiply_constant + add_constant != 0 )
									product *= pow( (double)(s[number][i])*multiply_constant + add_constant, 1./alpha_sum );
							}
							s[layers_number_visible][i] =
								(guchar) CLAMP( ROUND( product ), 0, 255 );
							break;
						case MEDIAN: /* Median */
							/* Note: Alpha channel weighted pixel values don't make any sense on median */
							/* However, we can skip all invisible layers */
							{
								visible_pixels = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 ) visible_pixels += 1;
									}
									else visible_pixels += 1;
								}

								gint vals[visible_pixels]; /* Negative values are allowed! */

								/* Sort all values and take the one in the middle (or the arithmetic mean value of the two values in the middle) */
								visible_number = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 )
										{
											vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
											visible_number++;
										}
									}
									else
									{
										vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
										visible_number++;
									}
								}
								qsort( vals, visible_pixels, sizeof(guint), (void *)compare );
								s[layers_number_visible][i] =
									(guchar) CLAMP( ROUND( median( vals, 0, visible_pixels-1 ) ), 0, 255 );
							}
							break;
						case SIGMA_MEDIAN_ARITHMETIC: /* Arithmetic mean of x sigma around median */
							/* Median of all values, that are inside the range of (Median of all values)-x*sigma and (Median of all values)+x*sigma */
							{
								visible_pixels = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 ) visible_pixels += 1;
									}
									else visible_pixels += 1;
								}

								gint vals[visible_pixels]; /* Negative values are allowed! */

								sum = 0.;
								sum_square = 0.;
								visible_number = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 )
										{
											sum += ( (double)(s[number][i])*multiply_constant + add_constant );
											sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
											vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
											visible_number++;
										}
									}
									else
									{
										sum += ( (double)(s[number][i])*multiply_constant + add_constant );
										sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
										vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
										visible_number++;
									}
								}
								qsort( vals, visible_pixels, sizeof(guint), (void *)compare );
								sigma = sqrt( ( visible_pixels*sum_square - sum*sum )/( visible_pixels * ( visible_pixels - 1 ) ) );
								median_tmp = median( vals, 0, visible_pixels-1 );

								start = 0;
								end = visible_pixels-1;
								for ( j=0; j < visible_pixels; j++ )
								{
									if ( vals[j] >= ROUND( median_tmp - sigma_merge_method*sigma ) )
									{
										start = j;
										break;
									}
								}
								for ( j = visible_pixels-1; j >= 0; j-- )
								{
									if ( vals[j] <= ROUND( median_tmp + sigma_merge_method*sigma ) )
									{
										end = j;
										break;
									}
								}
								if ( start > end ) end = start;

								sum = 0.;
								for ( j=start; j<=end; j++ ) sum += vals[j];
								sum /= end-start+1;

								s[layers_number_visible][i] =
									(guchar) CLAMP( ROUND( sum ), 0, 255 );
							}
							break;
						case SIGMA_MEDIAN_1: /* Sigma median (1 pass) */
							/* Median of all values, that are inside the range of (Median of all values)-sigma and (Median of all values)+sigma */
							{
								visible_pixels = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 ) visible_pixels += 1;
									}
									else visible_pixels += 1;
								}

								gint vals[visible_pixels]; /* Negative values are allowed! */

								sum = 0.;
								sum_square = 0.;
								visible_number = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 )
										{
											sum += ( (double)(s[number][i])*multiply_constant + add_constant );
											sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
											vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
											visible_number++;
										}
									}
									else
									{
										sum += ( (double)(s[number][i])*multiply_constant + add_constant );
										sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
										vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
										visible_number++;
									}
								}
								qsort( vals, visible_pixels, sizeof(guint), (void *)compare );
								sigma = sqrt( ( visible_pixels*sum_square - sum*sum )/( visible_pixels * ( visible_pixels - 1 ) ) );
								median_tmp = median( vals, 0, visible_pixels-1 );

								start = 0;
								end = visible_pixels-1;
								for ( j=0; j < visible_pixels; j++ )
								{
									if ( vals[j] >= ROUND( median_tmp - sigma_merge_method*sigma ) )
									{
										start = j;
										break;
									}
								}
								for ( j = visible_pixels-1; j >= 0; j-- )
								{
									if ( vals[j] <= ROUND( median_tmp + sigma_merge_method*sigma ) )
									{
										end = j;
										break;
									}
								}
								if ( start > end ) end = start;

								s[layers_number_visible][i] =
									(guchar) CLAMP( ROUND( median( vals, start, end ) ), 0, 255 );
							}
							break;
						case SIGMA_MEDIAN_2: /* Sigma median (2 passes) */
							/* Median of all values, that are inside the range of (Sigma Median)-sigma and (Sigma Median)+sigma */
							{
								visible_pixels = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 ) visible_pixels += 1;
									}
									else visible_pixels += 1;
								}

								gint vals[visible_pixels]; /* Negative values are allowed! */

								sum = 0.;
								sum_square = 0.;
								visible_number = 0;
								for ( number = 0; number < layers_number_visible; number++ )
								{
									if ( rgn[number]->bpp > bpp_effective )
									{
										if ( s[number][rgn[number]->bpp-1] > 0 )
										{
											sum += ( (double)(s[number][i])*multiply_constant + add_constant );
											sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
											vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
											visible_number++;
										}
									}
									else
									{
										sum += ( (double)(s[number][i])*multiply_constant + add_constant );
										sum_square += SQR( (double)(s[number][i])*multiply_constant + add_constant );
										vals[visible_number] = ( (double)(s[number][i])*multiply_constant + add_constant );
										visible_number++;
									}
								}
								qsort( vals, visible_pixels, sizeof(guint), (void *)compare );
								median_tmp = median( vals, 0, visible_pixels-1 );

								sigma = sqrt( ( visible_pixels*sum_square - sum*sum )/( visible_pixels * ( visible_pixels - 1 ) ) );

								start = 0;
								end = visible_pixels-1;
								for ( j=0; j < visible_pixels; j++ )
								{
									if ( vals[j] >= ROUND( median_tmp - sigma_merge_method*sigma ) )
									{
										start = j;
										break;
									}
								}
								for ( j = visible_pixels-1; j >= 0; j-- )
								{
									if ( vals[j] <= ROUND( median_tmp + sigma_merge_method*sigma ) )
									{
										end = j;
										break;
									}
								}
								if ( start > end ) end = start;

								if ( start < end )
								{
									median_tmp = median( vals, start, end );

									sum = 0.;
									sum_square = 0.;
									for ( j = start; j <= end; j++ )
									{
										sum += vals[j];
										sum_square += vals[j]*vals[j];
									}
									sigma = sqrt( ( ( end - start + 1 )*sum_square - sum*sum )/( ( end - start + 1 ) * ( end - start ) ) );

									start = 0;
									end = visible_pixels-1;
									for ( j=0; j < visible_pixels; j++ )
									{
										if ( vals[j] >= ROUND( median_tmp - sigma_merge_method*sigma ) )
										{
											start = j;
											break;
										}
									}
									for ( j = visible_pixels-1; j >= 0; j-- )
									{
										if ( vals[j] <= ROUND( median_tmp + sigma_merge_method*sigma ) )
										{
											end = j;
											break;
										}
									}
									if ( start > end ) end = start;
								}

								s[layers_number_visible][i] =
									(guchar) CLAMP( ROUND( median( vals, start, end ) ), 0, 255 );
							}
							break;
						default:
							break;
					}
				}
				/* If destination region has alpha channel, set it to visible */
				if ( region_destination->bpp > bpp_effective ) s[layers_number_visible][region_destination->bpp-1] = 255;

				for ( number = 0; number <= layers_number_visible; number++ ) s[number] += rgn[number]->bpp;
			}

			for ( number = 0; number <= layers_number_visible; number++ ) src[number] += rgn[number]->rowstride;
		}

		progress += rgn[layers_number_visible]->w * rgn[layers_number_visible]->h;
		progress_counter++;
		if ( progress_counter%25 == 0 ) /* only update progress every xth tile to avoid slowdown */
			gimp_progress_update( (double)(progress) /
				( rgn[layers_number_visible]->drawable->width * rgn[layers_number_visible]->drawable->height ) );
	}

	for ( number = 0; number < layers_number_visible; number++ )
	{
		g_free( rgn[number] );
	}
	g_free( rgn );
}

/*
merge layers
*/
static void merge_layers()
{
	gint32 layers_number, number;
	gint32 *layers;

	layers = gimp_image_get_layers( image_id, &layers_number );

	gimp_progress_init( _("Merge Layers...") );

	/* Number of visible layers */
	gint32 layers_number_visible = layers_number;
	if ( parameters.visible_only )
	{
		for( number = 0; number < layers_number; number++ )
		{
			if ( !gimp_drawable_get_visible( layers[number] ) )
			{
				layers_number_visible--;
			}
		}
	}

	if ( layers_number_visible < 1 ) return;

	gint32 layers_visible[layers_number_visible];
	gint32 number_visible = 0;

	for( number = layers_number - 1; number >= 0; number-- )
	{
		if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
		{
			layers_visible[number_visible] = layers[number];
			number_visible++;
		}
	}

	gint32 x_start, y_start, x_size, y_size;
	if ( gimp_drawable_mask_bounds( layers_visible[0], &x_start, &y_start, &x_size, &y_size ) )
	{
		x_size = x_size - x_start;
		y_size = y_size - y_start;
	}
	else
	{
		x_start = 0;
		y_start = 0;
		x_size = gimp_drawable_width( layers_visible[0] );
		y_size = gimp_drawable_height( layers_visible[0] );
	}

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( layers[0] ) )
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

	/* destination layer */
	gint32 image_destination = gimp_image_new( x_size, y_size, image_type );

	gimp_display_new( image_destination );
	gchar *name = "";
	switch( parameters.merge_method )
	{
		case ARITHMETIC:
			name = _("Arithmetic mean");
			break;
		case GEOMETRIC:
			name = _("Geometric mean");
			break;
		case MEDIAN:
			name = _("Median");
			break;
		case SIGMA_MEDIAN_ARITHMETIC:
			name = _("Sigma median / arithmetic");
			break;
		case SIGMA_MEDIAN_1:
			name = _("Sigma median (1 pass)");
			break;
		case SIGMA_MEDIAN_2:
			name = _("Sigma median (2 passes)");
			break;
		default:
			break;
	}
	gint32 layer_destination = gimp_layer_new( image_destination, name, gimp_image_width( image_destination ),
		gimp_image_height( image_destination ), layer_type, 100, GIMP_NORMAL_MODE );

	gimp_image_add_layer( image_destination, layer_destination, 0 );

	gimp_drawable_fill( layer_destination, GIMP_WHITE_FILL );

	gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
	gimp_drawable_merge_shadow( layer_destination, FALSE );
	gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );

	GimpPixelRgn region_destination;
	gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
		gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

	merge_selection( layers_visible, layers_number_visible, &region_destination, parameters.visible_only,
		x_start, y_start, x_size, y_size, parameters.merge_method, bpp, parameters.sigma_merge_method,
		parameters.multiply_constant, parameters.add_constant );

	gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
	gimp_drawable_merge_shadow( layer_destination, TRUE );
	gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );

	gimp_displays_flush();
}


/*
GUI
*/

GtkWidget *visible_layers;

/* Preview callback */
static void preview_callback( GimpPreview *preview, gpointer *data )
{
	GimpDrawable *drawable = gimp_drawable_preview_get_drawable( GIMP_DRAWABLE_PREVIEW( preview ) );

	gint32 x_pos, y_pos, x_size, y_size;
	gimp_preview_get_position( preview, &x_pos, &y_pos );
	gimp_preview_get_size( preview, &x_size, &y_size );

	gint32 layers_number, number;
	gint32 *layers;

	layers = gimp_image_get_layers( image_id, &layers_number );

	/* Number of visible layers */
	gint32 layers_number_visible = layers_number;
	if ( parameters.visible_only )
	{
		for( number = 0; number < layers_number; number++ )
		{
			if ( !gimp_drawable_get_visible( layers[number] ) )
			{
				layers_number_visible--;
			}
		}
	}

	if ( layers_number_visible < 1 ) return;

	gint32 layers_visible[layers_number_visible];
	gint32 number_visible = 0;

	for( number = layers_number - 1; number >= 0; number-- )
	{
		if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
		{
			layers_visible[number_visible] = layers[number];
			number_visible++;
		}
	}

	gint bpp = 0;
	GimpImageBaseType image_type;
	GimpImageType layer_type;
	switch( gimp_drawable_type( layers[0] ) )
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

	GimpPixelRgn region_destination;
	gimp_pixel_rgn_init( &region_destination, drawable, x_pos, y_pos, x_size, y_size, TRUE, TRUE );

	merge_selection( layers_visible, layers_number_visible, &region_destination, parameters.visible_only,
		x_pos, y_pos, x_size, y_size, parameters.merge_method, bpp, parameters.sigma_merge_method,
		parameters.multiply_constant, parameters.add_constant );

	gimp_drawable_flush( drawable );
	gimp_drawable_update( drawable->drawable_id, 0, 0, drawable->width, drawable->height );

	/* Read from shadow tiles ("FALSE, TRUE") */
	gimp_pixel_rgn_init( &region_destination, drawable, 0, 0, drawable->width, drawable->height, FALSE, TRUE );
	gimp_drawable_preview_draw_region( GIMP_DRAWABLE_PREVIEW( preview ), &region_destination );
}

static void visible_only_toggled( GtkWidget *check )
{
	parameters.visible_only = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
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

	dlg = gimp_dialog_new( _("Merge Layers"), "astro_merge_layers", NULL, 0,
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

	g_signal_connect( preview, "invalidated", G_CALLBACK( preview_callback ), NULL );

	sub_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( sub_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), sub_vbox, FALSE, FALSE, 0 );

	/* Generic settings */
	frame = gimp_frame_new( _("General Settings") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 1 );
	gtk_widget_show( frame );

	table = gtk_table_new( 1, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

		/* visible layers ? */
	visible_layers = gtk_check_button_new_with_label( _("Use visible layers only") );
	gtk_table_attach_defaults( GTK_TABLE( table ), visible_layers, 0, 2, 0, 1 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( visible_layers ), parameters.visible_only );
	g_signal_connect( visible_layers, "toggled", G_CALLBACK( visible_only_toggled ), NULL );
	g_signal_connect_swapped( visible_layers, "toggled", G_CALLBACK( gimp_preview_invalidate ), preview );
	gtk_widget_show( visible_layers );

	/* Merge methods + settings */
	frame = gimp_frame_new( _("Merge Method (Mean Value)") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 2 );
	gtk_widget_show( frame );

	table = gtk_table_new( 2, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	GtkWidget *merge_method = gimp_int_combo_box_new(
		_("Arithmetic mean"), ARITHMETIC,
		_("Geometric mean"), GEOMETRIC,
		_("Median"), MEDIAN,
		_("Arithmetic mean of X sigma around median"), SIGMA_MEDIAN_ARITHMETIC,
		_("Median of X sigma around median (1 pass)"), SIGMA_MEDIAN_1,
		_("Median of X sigma around median (2 passes)"), SIGMA_MEDIAN_2,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( merge_method ), parameters.merge_method );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( merge_method ), parameters.merge_method,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.merge_method );
	g_signal_connect_swapped( merge_method, "changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	gtk_table_attach( GTK_TABLE( table ), merge_method, 0, 2, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( merge_method );

		/* Sigma */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Sigma:"), 185, 75,
		parameters.sigma_merge_method, 0.5, 3., 0.1, 5, 3,
		TRUE, 0, 0, _("Clip X sigma around median (unused for arithmetic/geometric mean or median)"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.sigma_merge_method );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	/* Preprocess values */
	frame = gimp_frame_new( _("Preprocess Pixel Values Before Merging") );
	gtk_box_pack_start( GTK_BOX( sub_vbox ), frame, FALSE, FALSE, 2 );
	gtk_widget_show( frame );

	table = gtk_table_new( 2, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

		/* Multiply with constant */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 0,
		_("Multiply:"), 185, 75,
		parameters.multiply_constant, -10, 10, 0.1, 5, 3,
		TRUE, 0, 0, _("Multiply each pixel value with this constant before processing"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.multiply_constant );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

		/* Add a constant */
	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Add a constant:"), 185, 75,
		parameters.add_constant, -255, 255, 1, 10, 0,
		TRUE, 0, 0, _("Add this constant to each pixel value before processing"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.add_constant );
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
