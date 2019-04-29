/***************************************************************************
 *   Copyright (C) 2008 by Georg Hennig                                    *
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

#include "plugin-intl.h"

#define PLUG_IN_NAME "astro-alignment-autopano"
#define PLUG_IN_VERSION "astro-alignment-autopano 0.7"
#define PLUG_IN_DATE "08.2008"

/*
parameters
*/
typedef struct
{
	/* generatekeys executable */
	gchar *generatekeys;
	/* autopano executable */
	gchar *autopano;
	/* directory to store tmp files */
	gchar *tmp_directory;

	/* generatekeys options */
	gint32 generatekeys_maximum_size;

	/* autopano options */
	/* Options */
	/* --ransac <on|off|1|0>   Switch RANSAC filtration on or off (default: on) */
	gboolean autopano_ransac;
	/* --maxmatches <matches>  Use no more than the given number of matches
                             (default: 16, use zero for unlimited) */
	gint32 autopano_maxmatches;
	/* --disable-areafilter    Do not use max-area filtration, which is default.
                             See manpage for details. */
	gboolean autopano_disable_areafilter;
	/* --integer-coordinates   Truncate match coordinates to integer numbers. */
	/* --absolute-pathnames <on|off|1|0>   Use the absolute pathname of the image
                             file in the PTO output file. Disabled by default. */

	/* Alignment options */
	/* --align                 Automatically pre-align images in PTO file. */
	/* --bottom-is-left */
	/* --bottom-is-right       Use in case the automatic algorithm fails. */
	/* --generate-horizon <c>  Generate up to 'c' horizon lines. */

	/* Refinement options */
	/* --refine                Refine the found control points using the
                             original images. */
	gboolean autopano_refine;
	/* --refine-by-middle      Use the best middle point to refine (default). */
	/* --refine-by-mean        Use the mean of the patches control points. */
	gint32 autopano_refinement_method;
	/* --keep-unrefinable <on|off|1|0>
                             Keep unrefinable matches (default: on). */
	gboolean autopano_keep_unrefinable;

	gboolean crop_overlap;
} tparameter;

enum REFINEMENT_METHOD
{
	REFINEMENT_MIDDLE = 0,
	REFINEMENT_MEAN
};

static tparameter parameters =
{
	"generatekeys", /* gchar *generatekeys; */
	"autopano", /* gchar *autopano; */
	"/tmp", /* gchar *tmp_directory; */
	800, /* gint32 generatekeys_maximum_size; */
	TRUE, /* gboolean autopano_ransac; */
	16, /* gint32 autopano_maxmatches; */
	TRUE, /* gboolean autopano_disable_areafilter; */
	TRUE, /* gboolean autopano_refine; */
	REFINEMENT_MIDDLE, /* gint32 autopano_refinement_method; */
	TRUE, /* gboolean autopano_keep_unrefinable; */

	FALSE /* gboolean crop_overlap; */
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
Open aligned layers
*/
static void open_aligned_layers();

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
		{ GIMP_PDB_STRING, "generatekeys", "generatekeys executable" },
		{ GIMP_PDB_STRING, "autopano", "autopano executable" },
		{ GIMP_PDB_STRING, "tmp_directory", "" },
		{ GIMP_PDB_INT32, "generatekeys_maximum_size", "" },
		{ GIMP_PDB_INT32, "autopano_ransac", "" },
		{ GIMP_PDB_INT32, "autopano_maxmatches", "" },
		{ GIMP_PDB_INT32, "autopano_disable_areafilter", "" },
		{ GIMP_PDB_INT32, "autopano_refine", "" },
		{ GIMP_PDB_INT32, "autopano_refinement_method", "" },
		{ GIMP_PDB_INT32, "autopano_keep_unrefinable", "" },
		{ GIMP_PDB_INT32, "crop_overlap", "" }
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

	gimp_install_procedure( "plug_in_"PLUG_IN_NAME,
		"autopano layer alignment",
		_("This plug-in opens images as layers and translates/rotates them so they fit best."),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Align images with autopano..."),
		"RGB*,GRAY*",
		GIMP_PLUGIN,
		nparams,
		nreturn_vals,
		params,
		return_vals );

	gimp_plugin_menu_register( "plug_in_"PLUG_IN_NAME, _("<Toolbox>/Xtns") );
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
			gimp_get_data( "plug_in_"PLUG_IN_NAME, &parameters );

			if ( !dialog() ) return;

			gimp_set_data( "plug_in_"PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 12 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				parameters.generatekeys = param[1].data.d_string;
				parameters.autopano = param[2].data.d_string;
				parameters.tmp_directory = param[3].data.d_string;

				parameters.generatekeys_maximum_size = param[4].data.d_int32;

				parameters.autopano_ransac = param[5].data.d_int32;
				parameters.autopano_maxmatches = param[6].data.d_int32;
				parameters.autopano_disable_areafilter = param[7].data.d_int32;
				parameters.autopano_refine = param[8].data.d_int32;
				parameters.autopano_refinement_method = param[9].data.d_int32;
				parameters.autopano_keep_unrefinable = param[10].data.d_int32;
				parameters.crop_overlap = param[11].data.d_int32;
			}
			break;
		case GIMP_RUN_WITH_LAST_VALS:
			gimp_get_data( "plug_in_"PLUG_IN_NAME, &parameters );
			break;
		default:
			status = GIMP_PDB_CALLING_ERROR;
			break;
	}

	if( status == GIMP_PDB_SUCCESS )
	{
		open_aligned_layers();
	}

	values[0].data.d_status = status;
}

/*
open aligned layers
*/
static void open_aligned_layers()
{
	GtkWidget *file_open_dialog;

	file_open_dialog = gtk_file_chooser_dialog_new( _("Open image files"),
		NULL,
		GTK_FILE_CHOOSER_ACTION_OPEN,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		GTK_STOCK_OPEN, GTK_RESPONSE_ACCEPT,
		NULL );
	gtk_file_chooser_set_select_multiple( file_open_dialog, TRUE );

	GSList *filename_list = NULL;

	if( gtk_dialog_run( GTK_DIALOG( file_open_dialog ) ) == GTK_RESPONSE_ACCEPT )
	{
// 		char *filename;
// 		filename = gtk_file_chooser_get_filename( GTK_FILE_CHOOSER( dialog ) );
// 		open_file (filename);
// 		g_free (filename);

		filename_list = gtk_file_chooser_get_filenames( file_open_dialog );
	}
	else
	{
		printf( "No files selected. Aborting.\n" );
		return;
	}

	gtk_widget_destroy( file_open_dialog );




g_slist_free(), and the filenames with g_free().





// 	gimp_image_undo_group_start( image_id );

// 	gint32 layers_number, number;
// 	gint32 *layers;
// 
// 	layers = gimp_image_get_layers( image_id, &layers_number );
// 
// 	if ( layers_number < 2 )
// 	{
// 		gimp_image_undo_group_end( image_id );
// 
// 		gimp_displays_flush();
// 		return;
// 	}
// 
// 	gint layer_move_x[layers_number];
// 	gint layer_move_y[layers_number];
// 	gdouble layer_quality[layers_number];
// 	gint move_x = 0, move_y = 0;
// 	gdouble quality = 0.;
// 
// 	gimp_progress_init( _("Align Layers...") );
// 	gint sel_pos_x, sel_pos_y, sel_width, sel_height;
// 
// 	gimp_image_resize_to_layers( image_id );
// 	for( number = 0; number < layers_number; number++ )
// 	{
// 		gimp_layer_resize_to_image_size( layers[number] );
// 	}
// 
// /* FIXME: first visible layer or selected layer? */
// 	gint32 active_layer = /*gimp_image_get_active_layer( image_id );
// 	if ( active_layer == -1 ) active_layer =*/ 0;
// 
// 	if ( parameters.scale_percent > 100 )
// 	{
// 		gint new_size_x = (parameters.scale_percent*gimp_drawable_width( layers[active_layer] ))/100;
// 		gint new_size_y = (parameters.scale_percent*gimp_drawable_height( layers[active_layer] ))/100;
// 
// 		gimp_image_scale( image_id, new_size_x, new_size_y ); /* This also scales the selection! */
// 	}
// 
// 	if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
// 	{
// 		sel_width = sel_width - sel_pos_x;
// 		sel_height = sel_height - sel_pos_y;
// 	}
// 	else
// 	{
// 		sel_pos_x = 0;
// 		sel_pos_y = 0;
// 		sel_width = gimp_drawable_width( layers[active_layer] );
// 		sel_height = gimp_drawable_height( layers[active_layer] );
// 	}
// 
// 	guchar *reference_data = NULL;
// 	if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
// 	{
// 		if ( sel_width%2 ) sel_width--;
// 		if ( sel_height%2 ) sel_height--;
// 
// 		reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );
// 
// 		GimpPixelRgn region_source;
// 		gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
// 			0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );
// 
// 		gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
// 	}
// 
// 	for( number=layers_number - 1; number >= 0; number-- )
// 	{
// 		if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
// 		{
// 			if ( number == active_layer && parameters.alignment_method == CROSS_CORRELATION ) /* Don't skip this for CROSS_CORRELATION_FFT ! */
// 			{
// 				move_x = 0;
// 				move_y = 0;
// 				quality = 1.;
// 			}
// 			else
// 			{
// 				get_center( layers[number], parameters.alignment_method, sel_pos_x, sel_pos_y, sel_width, sel_height,
// 					parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
// 					&move_x, &move_y, &quality );
// 			}
// 
// 			layer_move_x[number] = move_x;
// 			layer_move_y[number] = move_y;
// 			layer_quality[number] = quality;
// 		}
// 		else
// 		{
// 			layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
// 			layer_move_y[number] = 0;
// 			layer_quality[number] = 0;
// 		}
// 
// 		gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
// 	}
// 	free( reference_data );
// 
// 	gimp_progress_init( _("Moving layers by estimated values...") );
// 	char buffer[1024];
// /* Calculate the "quality" always in the same way. Contrast? */
// 	sprintf( buffer, _("Quality: %.5f"), layer_quality[0] );
// 	gimp_drawable_set_name( layers[0], buffer );
// 	for( number=layers_number - 1; number > 0; number-- )
// 	{
// 		if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
// 		{
// 			if ( parameters.invisible_remove )
// 			{
// 				gimp_image_remove_layer( image_id, layers[number] );
// 			}
// 		}
// 		else
// 		{
// 			/* Sanity check of calculated move values */
// 			if ( abs( layer_move_x[0]-layer_move_x[number] ) <= sel_width*parameters.cross_correlation/100 &&
// 				abs( layer_move_y[0]-layer_move_y[number] ) <= sel_height*parameters.cross_correlation/100 )
// 			{
// 				gimp_layer_translate( layers[number], layer_move_x[0]-layer_move_x[number], layer_move_y[0]-layer_move_y[number] );
// 				sprintf( buffer, _("Quality: %.5f"), layer_quality[number] );
// 				gimp_drawable_set_name( layers[number], buffer );
// 			}
// 			else
// 			{
// 				g_message( _("Requested invalid move values. Not moving layer!") );
// 			}
// 		}
// 
// 		gimp_progress_update( ((gdouble)layers_number-number) / (layers_number-1) );
// 	}
// 
// 	gimp_image_resize_to_layers( image_id ); /* Resized and moved layers require a bigger image */
// 
// 	/* Position of selection might have changed during resize */
// 	gint center_x = sel_pos_x + sel_width/2;
// 	gint center_y = sel_pos_y + sel_height/2;
// 
// 	gint x_min = 0;
// 	gint y_min = 0;
// 	gint x_max = gimp_image_width( image_id );
// 	gint y_max = gimp_image_height( image_id );
// 	gint x_off, y_off;
// 
// 	if ( parameters.trim_image && !parameters.two_star )
// 	{
// 		for( number=layers_number - 1; number >= 0; number-- )
// 		{
// 			if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
// 			{
// 				gimp_drawable_offsets( layers[number], &x_off, &y_off );
// 				if ( x_min < x_off ) x_min = x_off;
// 				if ( y_min < y_off ) y_min = y_off;
// 				if ( x_max > x_off + gimp_drawable_width( layers[number] ) ) x_max = x_off + gimp_drawable_width( layers[number] );
// 				if ( y_max > y_off + gimp_drawable_height( layers[number] ) ) y_max = y_off + gimp_drawable_height( layers[number] );
// 			}
// 		}
// 
// 		gimp_image_crop( image_id, x_max - x_min, y_max - y_min, x_min, y_min );
// 	}
// 
// 	if ( parameters.two_star )
// 	{
// 		gimp_progress_init( _("Align Layers...") );
// 
// 		GtkWidget *dialog_widget;
// 		dialog_widget = gtk_message_dialog_new( 0,
// 			GTK_DIALOG_DESTROY_WITH_PARENT,
// 			GTK_MESSAGE_INFO,
// 			GTK_BUTTONS_OK_CANCEL,
// 			_("Please select another area for de-rotating layers.") );
// 		gimp_window_set_transient( GTK_WINDOW( dialog_widget ) );
// 
// 		gint result = gtk_dialog_run( GTK_DIALOG( dialog_widget ) );
// 		gtk_widget_destroy( dialog_widget );
// 		if ( result == GTK_RESPONSE_OK )
// 		{
// 			active_layer = /*gimp_image_get_active_layer( image_id );
// 			if ( active_layer == -1 ) active_layer =*/ 0;
// 
// 			if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
// 			{
// 				sel_width = sel_width - sel_pos_x;
// 				sel_height = sel_height - sel_pos_y;
// 			}
// 			else
// 			{
// 				sel_pos_x = 0;
// 				sel_pos_y = 0;
// 				sel_width = gimp_drawable_width( layers[active_layer] );
// 				sel_height = gimp_drawable_height( layers[active_layer] );
// 			}
// 
// 			if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
// 			{
// 				if ( sel_width%2 ) sel_width--;
// 				if ( sel_height%2 ) sel_height--;
// 
// 				reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );
// 
// 				GimpPixelRgn region_source;
// 				gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
// 					0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );
// 
// 				gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
// 			}
// 
// 			gint center_x_corrected, center_y_corrected, sel_pos_x_corrected, sel_pos_y_corrected;
// 			gint x_off_0, y_off_0;
// 			gimp_drawable_offsets( layers[0], &x_off_0, &y_off_0 );
// 
// 			for( number=layers_number - 1; number >= 0; number-- )
// 			{
// 				gimp_drawable_offsets( layers[number], &x_off, &y_off );
// 				sel_pos_x_corrected = sel_pos_x + x_off_0 - x_off;
// 				sel_pos_y_corrected = sel_pos_y + y_off_0 - y_off;
// 
// 				if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
// 				{
// 					if ( number == active_layer && ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT ) )
// 					{
// 						move_x = 0;
// 						move_y = 0;
// 						quality = 1.;
// 					}
// 					else
// 					{
// 						get_center( layers[number], parameters.alignment_method, sel_pos_x_corrected, sel_pos_y_corrected, sel_width, sel_height,
// 							parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
// 							&move_x, &move_y, &quality );
// 					}
// 
// 					layer_move_x[number] = move_x;
// 					layer_move_y[number] = move_y;
// 					layer_quality[number] = quality;
// 				}
// 				else
// 				{
// 					layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
// 					layer_move_y[number] = 0;
// 					layer_quality[number] = 0;
// 				}
// 
// 				gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
// 			}
// 			free( reference_data );
// 
// 			gimp_progress_init( _("Rotating layers by estimated values...") );
// 			gint32 selection = gimp_selection_save( image_id );
// 			gimp_selection_none( image_id );
// 
// 			for( number=layers_number - 1; number > 0; number-- )
// 			{
// 				if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
// 				{
// 					if ( parameters.invisible_remove )
// 					{
// 						gimp_image_remove_layer( image_id, layers[number] );
// 					}
// 				}
// 				else
// 				{
// 					gimp_drawable_offsets( layers[number], &x_off, &y_off );
// 					center_x_corrected = center_x + x_off_0 - x_off;
// 					center_y_corrected = center_y + y_off_0 - y_off;
// 
// 					gdouble distance_p0_pnumber =
// 						sqrt( (gdouble)((layer_move_x[0]-layer_move_x[number])*(layer_move_x[0]-layer_move_x[number]) +
// 							(layer_move_y[0]-layer_move_y[number])*(layer_move_y[0]-layer_move_y[number]) ) );
// 					gdouble distance_p0_center =
// 						sqrt( (gdouble)( (center_x_corrected-sel_pos_x-layer_move_x[0])*(center_x_corrected-sel_pos_x-layer_move_x[0]) +
// 							(center_y_corrected-sel_pos_y-layer_move_y[0])*(center_y_corrected-sel_pos_y-layer_move_y[0]) ) );
// 
// 					gdouble angle = atan2( distance_p0_pnumber, distance_p0_center );
// 
// 					/* (point0_x-center_x)*(pointnumber_y-point0_y) - (point0_y-center_y)*(pointnumber_x-point0_x) */
// 					gdouble crossproduct_3rd_component =
// 						(layer_move_x[0]+sel_pos_x-center_x_corrected)*(layer_move_y[number]-layer_move_y[0])-
// 						(layer_move_y[0]+sel_pos_y-center_y_corrected)*(layer_move_x[number]-layer_move_x[0]);
// 
// 					if ( crossproduct_3rd_component > 0 ) angle = -angle;
// 
// 					printf( "alignment: rotating by angle=%.5f around center x=%d, y=%d...\n", angle, center_x_corrected, center_y_corrected );
// 
// 					gimp_drawable_transform_rotate( layers[number], angle, FALSE /* Auto-center */,
// 						center_x_corrected, center_y_corrected, GIMP_TRANSFORM_FORWARD,
// 						GIMP_INTERPOLATION_CUBIC, FALSE /* Supersample */, 3 /* recursion */,
// 						( parameters.trim_image ? GIMP_TRANSFORM_RESIZE_CROP : GIMP_TRANSFORM_RESIZE_ADJUST )/* clip results */ );
// 				}
// 
// 				gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
// 			}
// 			gimp_selection_load( selection );
// 
// 				active_layer = /*gimp_image_get_active_layer( image_id );
// 				if ( active_layer == -1 ) active_layer =*/ 0;
// 
// 				if ( gimp_drawable_mask_bounds( layers[active_layer], &sel_pos_x, &sel_pos_y, &sel_width, &sel_height ) )
// 				{
// 					sel_width = sel_width - sel_pos_x;
// 					sel_height = sel_height - sel_pos_y;
// 				}
// 				else
// 				{
// 					sel_pos_x = 0;
// 					sel_pos_y = 0;
// 					sel_width = gimp_drawable_width( layers[active_layer] );
// 					sel_height = gimp_drawable_height( layers[active_layer] );
// 				}
// 
// 				if ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT )
// 				{
// 					if ( sel_width%2 ) sel_width--;
// 					if ( sel_height%2 ) sel_height--;
// 
// 					reference_data = malloc( gimp_drawable_bpp( layers[active_layer] )*sel_width*sel_height*sizeof(guchar) );
// 
// 					GimpPixelRgn region_source;
// 					gimp_pixel_rgn_init( &region_source, gimp_drawable_get( layers[active_layer] ),
// 						0, 0, gimp_drawable_width( layers[active_layer] ), gimp_drawable_height( layers[active_layer] ), FALSE, FALSE );
// 
// 					gimp_pixel_rgn_get_rect( &region_source, reference_data, sel_pos_x, sel_pos_y, sel_width, sel_height );
// 				}
// 
// 				gimp_drawable_offsets( layers[0], &x_off_0, &y_off_0 );
// 
// 				for( number=layers_number - 1; number >= 0; number-- )
// 				{
// 					gimp_drawable_offsets( layers[number], &x_off, &y_off );
// 					sel_pos_x_corrected = sel_pos_x + x_off_0 - x_off;
// 					sel_pos_y_corrected = sel_pos_y + y_off_0 - y_off;
// 
// 					if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
// 					{
// 						if ( number == active_layer && ( parameters.alignment_method == CROSS_CORRELATION || parameters.alignment_method == CROSS_CORRELATION_FFT ) )
// 						{
// 							move_x = 0;
// 							move_y = 0;
// 							quality = 1.;
// 						}
// 						else
// 						{
// 							get_center( layers[number], parameters.alignment_method, sel_pos_x_corrected, sel_pos_y_corrected, sel_width, sel_height,
// 								parameters.search_radius, reference_data, gimp_drawable_bpp( layers[active_layer] ), parameters.cross_correlation,
// 								&move_x, &move_y, &quality );
// 						}
// 
// 						layer_move_x[number] = move_x;
// 						layer_move_y[number] = move_y;
// 						layer_quality[number] = quality;
// 					}
// 					else
// 					{
// 						layer_move_x[number] = 0; /* Skip invisible layers if user wants to */
// 						layer_move_y[number] = 0;
// 						layer_quality[number] = 0;
// 					}
// 
// 					gimp_progress_update( ((gdouble)layers_number-number) / layers_number );
// 				}
// 				free( reference_data );
// 
// 			gimp_progress_init( _("Moving layers by estimated values...") );
// 			sprintf( buffer, _("Quality: %.5f"), layer_quality[0] );
// 			gimp_drawable_set_name( layers[0], buffer );
// 			for( number=layers_number - 1; number > 0; number-- )
// 			{
// 				if ( !gimp_drawable_get_visible( layers[number] ) && parameters.visible_only )
// 				{
// 					if ( parameters.invisible_remove )
// 					{
// 						gimp_image_remove_layer( image_id, layers[number] );
// 					}
// 				}
// 				else
// 				{
// 					/* Sanity check of calculated move values */
// 					if ( abs( layer_move_x[0]-layer_move_x[number] ) <= sel_width*parameters.cross_correlation/100 &&
// 						abs( layer_move_y[0]-layer_move_y[number] ) <= sel_height*parameters.cross_correlation/100 )
// 					{
// 						gimp_layer_translate( layers[number], layer_move_x[0]-layer_move_x[number], layer_move_y[0]-layer_move_y[number] );
// 						sprintf( buffer, _("Quality: %.5f"), layer_quality[number] );
// 						gimp_drawable_set_name( layers[number], buffer );
// 					}
// 					else
// 					{
// 						g_message( _("Requested invalid move values. Not moving layer!") );
// 					}
// 				}
// 
// 				gimp_progress_update( ((gdouble)layers_number-number) / (layers_number-1) );
// 			}
// 
// 			if ( parameters.trim_image )
// 			{
// 				x_min = 0;
// 				y_min = 0;
// 				x_max = gimp_image_width( image_id );
// 				y_max = gimp_image_height( image_id );
// 				for( number=layers_number - 1; number >= 0; number-- )
// 				{
// 					if ( !parameters.visible_only || gimp_drawable_get_visible( layers[number] ) )
// 					{
// 						gimp_drawable_offsets( layers[number], &x_off, &y_off );
// 						if ( x_min < x_off ) x_min = x_off;
// 						if ( y_min < y_off ) y_min = y_off;
// 						if ( x_max > x_off + gimp_drawable_width( layers[number] ) ) x_max = x_off + gimp_drawable_width( layers[number] );
// 						if ( y_max > y_off + gimp_drawable_height( layers[number] ) ) y_max = y_off + gimp_drawable_height( layers[number] );
// 					}
// 				}
// 				gimp_image_crop( image_id, x_max - x_min, y_max - y_min, x_min, y_min );
// 			}
// 		}
// 	}

// 	gimp_image_undo_group_end( image_id );

	gimp_displays_flush();
}


/*
GUI
*/

GtkWidget *toggle_crop_overlap;
GtkWidget *toggle_autopano_ransac;
GtkWidget *toggle_autopano_disable_areafilter;
GtkWidget *toggle_autopano_refine;
GtkWidget *toggle_autopano_keep_unrefinable;

GtkWidget *spin_generatekeys_maximum_size;
GtkWidget *spin_autopano_maxmatches;

GtkWidget *refinement_method;

static void crop_overlap_toggled( GtkWidget *check )
{
	parameters.crop_overlap = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void autopano_ransac_toggled( GtkWidget *check )
{
	parameters.autopano_ransac = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void autopano_disable_areafilter_toggled( GtkWidget *check )
{
	parameters.autopano_disable_areafilter = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void autopano_refine_toggled( GtkWidget *check )
{
	parameters.autopano_refine = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );

	gtk_widget_set_sensitive( toggle_autopano_keep_unrefinable, parameters.autopano_refine );
	gtk_widget_set_sensitive( refinement_method, parameters.autopano_refine );
}

static void autopano_keep_unrefinable_toggled( GtkWidget *check )
{
	parameters.autopano_keep_unrefinable = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void generatekeys_maximum_size_spin( GtkWidget *spin )
{
	parameters.generatekeys_maximum_size = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void autopano_maxmatches_spin( GtkWidget *spin )
{
	parameters.autopano_maxmatches = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}


/*
main dialog
*/
static gint dialog()
{
	GtkWidget *dlg;
	GtkWidget *main_hbox;
	GtkWidget *left_vbox;
	GtkWidget *right_vbox;
	GtkWidget *frame;
	GtkWidget *table;
	GtkWidget *button;
	GtkWidget *label;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( PLUG_IN_VERSION, PLUG_IN_NAME, NULL, 0,
		gimp_standard_help_func, "plug-in-"PLUG_IN_NAME,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		_("Open & Align"), GTK_RESPONSE_OK,
		NULL);

	/* General layout */
	main_hbox = gtk_hbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( main_hbox ), 8 );
	gtk_container_add( GTK_CONTAINER( GTK_DIALOG( dlg )->vbox ), main_hbox );

	left_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( left_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), left_vbox, FALSE, FALSE, 0 );

	GtkWidget *vertical_line = gtk_vseparator_new();
	gtk_box_pack_start( GTK_BOX( main_hbox ), vertical_line, FALSE, FALSE, 0 );
	gtk_widget_show( vertical_line );

	right_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( right_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), right_vbox, FALSE, FALSE, 0 );

	/* File selection */
	GtkWidget *file_list = gtk_tree_view_new();
	gtk_box_pack_start( GTK_BOX( left_vbox ), file_list, FALSE, FALSE, 0 );
	gtk_widget_show( file_list );

	/* Generic settings */
	frame = gimp_frame_new( _("Generic Settings") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 4, 2, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("generatekeys executable:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( label );

	button = gtk_file_chooser_button_new( _("Select a file"), GTK_FILE_CHOOSER_ACTION_OPEN );
	gtk_file_chooser_set_filename( GTK_FILE_CHOOSER( button ), "/usr/bin/generatekeys" );
	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );

	label = gtk_label_new( _("autopano executable:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 1, 2, GTK_FILL, 0, 0, 0);
	gtk_widget_show( label );

	button = gtk_file_chooser_button_new( _("Select a file"), GTK_FILE_CHOOSER_ACTION_OPEN );
	gtk_file_chooser_set_filename( GTK_FILE_CHOOSER( button ), "/usr/bin/autopano" );
	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 1, 2, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );

	label = gtk_label_new( _("Temporary folder:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0);
	gtk_widget_show( label );

	button = gtk_file_chooser_button_new( _("Select a folder"), GTK_FILE_CHOOSER_ACTION_SELECT_FOLDER );
	gtk_file_chooser_set_filename( GTK_FILE_CHOOSER( button ), "/tmp" );
	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 2, 3, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );

	toggle_crop_overlap = gtk_check_button_new_with_label( _("Crop images to overlap") );
	gtk_table_attach( GTK_TABLE( table ), toggle_crop_overlap, 0, 2, 3, 4, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( toggle_crop_overlap ), parameters.crop_overlap );
	g_signal_connect( toggle_crop_overlap, "toggled", G_CALLBACK( crop_overlap_toggled ), NULL );
	gtk_widget_show( toggle_crop_overlap );

	/* generatekeys options */
	frame = gimp_frame_new( _("generatekeys Settings") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 1, 2, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("Reduce image size below:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	spin_generatekeys_maximum_size = gtk_spin_button_new_with_range( 100, 10000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin_generatekeys_maximum_size ), parameters.generatekeys_maximum_size );
	g_signal_connect( spin_generatekeys_maximum_size, "value_changed", G_CALLBACK( generatekeys_maximum_size_spin ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0);
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), spin_generatekeys_maximum_size, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin_generatekeys_maximum_size );

	/* autopano options */
	frame = gimp_frame_new( _("autopano Settings") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 7, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	toggle_autopano_ransac = gtk_check_button_new_with_label( _("Ransac filtration") );
	gtk_table_attach( GTK_TABLE( table ), toggle_autopano_ransac, 0, 3, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( toggle_autopano_ransac ), parameters.autopano_ransac );
	g_signal_connect( toggle_autopano_ransac, "toggled", G_CALLBACK( autopano_ransac_toggled ), NULL );
	gtk_widget_show( toggle_autopano_ransac );

	label = gtk_label_new( _("Maximum number of matching points:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	spin_autopano_maxmatches = gtk_spin_button_new_with_range( 0, 200, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin_autopano_maxmatches ), parameters.autopano_maxmatches );
	g_signal_connect( spin_autopano_maxmatches, "value_changed", G_CALLBACK( autopano_maxmatches_spin ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 2, 1, 2, GTK_FILL, 0, 0, 0);
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), spin_autopano_maxmatches, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin_autopano_maxmatches );

	toggle_autopano_disable_areafilter = gtk_check_button_new_with_label( _("Disable areafilter") );
	gtk_table_attach( GTK_TABLE( table ), toggle_autopano_disable_areafilter, 0, 3, 2, 3, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( toggle_autopano_disable_areafilter ), parameters.autopano_disable_areafilter );
	g_signal_connect( toggle_autopano_disable_areafilter, "toggled", G_CALLBACK( autopano_disable_areafilter_toggled ), NULL );
	gtk_widget_show( toggle_autopano_disable_areafilter );

	toggle_autopano_refine = gtk_check_button_new_with_label( _("Refine") );
	gtk_table_attach( GTK_TABLE( table ), toggle_autopano_refine, 0, 3, 3, 4, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( toggle_autopano_refine ), parameters.autopano_refine );
	g_signal_connect( toggle_autopano_refine, "toggled", G_CALLBACK( autopano_refine_toggled ), NULL );
	gtk_widget_show( toggle_autopano_refine );

	refinement_method = gimp_option_menu_new2( FALSE, G_CALLBACK( gimp_menu_item_update ),
		&parameters.autopano_refinement_method, (gpointer) parameters.autopano_refinement_method,
		_("Refine by middle"), (gpointer) REFINEMENT_MIDDLE, NULL,
		_("Refine by mean"), (gpointer) REFINEMENT_MEAN, NULL,
		NULL );
	gtk_table_attach( GTK_TABLE( table ), refinement_method, 0, 3, 4, 5, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( refinement_method );

	toggle_autopano_keep_unrefinable = gtk_check_button_new_with_label( _("Keep unrefinable") );
	gtk_table_attach( GTK_TABLE( table ), toggle_autopano_keep_unrefinable, 1, 2, 5, 6, GTK_EXPAND | GTK_FILL, 0, 0, 0 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( toggle_autopano_keep_unrefinable ), parameters.autopano_keep_unrefinable );
	g_signal_connect( toggle_autopano_keep_unrefinable, "toggled", G_CALLBACK( autopano_keep_unrefinable_toggled ), NULL );
	gtk_widget_show( toggle_autopano_keep_unrefinable );

	gtk_widget_show( left_vbox );
	gtk_widget_show( right_vbox );
	gtk_widget_show( main_hbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
