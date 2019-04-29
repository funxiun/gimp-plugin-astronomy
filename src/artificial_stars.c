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
gimp artificial stars plug-in
(C) Georg Hennig <georg.hennig@web.de>
Creates an artificial star distribution (background stars, object stars -
distributed randomly, with a radial dependency, or as a globular
cluster, and foreground stars.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "artificial_stars_temperature.h"

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-artificial-stars"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

enum PSF
{
	PSF_DELTA_PEAK = 0,
	PSF_GAUSS,
	PSF_GAUSS_4_DIFFRACTIONS,
	PSF_GAUSS_6_DIFFRACTIONS,
	PSF_GAUSS_12_DIFFRACTIONS
};

enum DENSITY
{
	DENSITY_RANDOM = 0,
	DENSITY_GAUSS,
	DENSITY_PLUMMER
};

enum SAMPLE_DISTRIBUTIONS
{
	STANDARD_HD = 0,
	STANDARD_LD,
	GLOBULAR_HD,
	GLOBULAR_LD,
	OPEN_HD,
	OPEN_LD
};

/*
parameters
*/
typedef struct
{
	gint32 psf;
	gdouble sigma_psf;
	gdouble diffraction_percentage;
	gdouble diffraction_angle;
	gdouble diffraction_length;
	gdouble diffraction_color;
	gint32 noise;
	gint32 background;
	gdouble burnout;
	gdouble shininess;

	gboolean split_layers;
	gboolean object_mask;

	guint32 random_seed;
	gboolean random_seed_bool;

	gint32 background_stars;
	gint32 object_stars;
	gint32 foreground_stars;

	gdouble background_brightness_mean;
	gdouble object_brightness_mean;
	gdouble foreground_brightness_mean;

	gdouble background_brightness_sigma;
	gdouble object_brightness_sigma;
	gdouble foreground_brightness_sigma;

	gint32 background_color_mean;
	gint32 object_color_mean;
	gint32 foreground_color_mean;

	gint32 background_color_sigma;
	gint32 object_color_sigma;
	gint32 foreground_color_sigma;

	gint32 object_x;
	gint32 object_y;

	gint32 object_density;
	gint32 object_radius;

	gboolean show_preview;
} tparameter;


static tparameter parameters =
{
	PSF_GAUSS,
	1.6,
	0.5,
	30.0,
	1.0,
	1.0,
	20,
	10,
	0.5,
	1.5,

/* split layers, object mask */
	FALSE,
	FALSE,

/* random seed */
	0,
	FALSE,

/* number */
	1000,
	500,
	5,

/* brightness */
	40,
	90,
	150,

/* brightness sigma */
	30,
	50,
	100,

/* color */
	7500,
	7500,
	7500,

/* color sigma */
	1000,
	1000,
	1000,

/* x, y */
	0,
	0,

/* density, radius */
	DENSITY_RANDOM,
	20,

	TRUE
};

/*
Prototypes
*/
gint32 background_stars_number = 0;
gdouble *background_stars = NULL;
gint32 object_stars_number = 0;
gdouble *object_stars = NULL;
gint32 foreground_stars_number = 0;
gdouble *foreground_stars = NULL;
gint32 all_stars_number = 0;
gdouble *all_stars = NULL;

gint32 sample_distribution = 0;

gdouble brightness_star_norm = 0.;
gdouble brightness_star_diffraction = 0.;

void star_new( gdouble **list, gint32 *list_number, const gint32 list_number_new );
void star_free( gdouble **list, gint32 *list_number );
void star_add( gdouble **list, const gint32 list_number, const gint32 list_position,
	const gdouble x, const gdouble y, const gdouble brightness, const gdouble color );
void star_get( gdouble **list, const gint32 list_position,
	gdouble *x, gdouble *y, gdouble *brightness, gdouble *color );

/*
communication to the GIMP
*/
static void query( void );

static void run( const gchar *name, gint nparams, const GimpParam *param,
	gint *nreturn_vals, GimpParam **return_vals );

/*
create stars
*/
void create_star_distribution();
void create_background_star_distribution();
void create_object_star_distribution();
void create_foreground_star_distribution();

void find_star_norm();
void find_star_norm_sort();

void create_stars_psf( gdouble *star_distribution, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, gdouble **list, const gint32 list_number, const gboolean show_progress );
void create_stars_noise( gdouble *star_distribution, gint32 x_size, gint32 y_size, const gboolean show_progress );

static void create_stars();

static void create_stars_selection( GimpPixelRgn *region_destination,
	const gboolean draw_background_stars, const gboolean draw_object_stars, const gboolean draw_foreground_stars,
	const gboolean show_progress );

static void recalculation_necessary();
static void recalculation_done();

/*
user interface
*/
GtkWidget *preview;
GtkWidget *progress_bar;

GtkWidget *recalculate_button;
gboolean recalculating = FALSE;
gboolean cancel_recalculation = FALSE;

gboolean distribution_dirty = TRUE;
gboolean dialog_closed = FALSE;

GtkWidget *background_number_spin;
GtkWidget *background_brightness_spin, *background_brightness_sigma_spin;
GtkWidget *background_color_spin, *background_color_sigma_spin;

GtkWidget *object_number_spin;
GtkWidget *object_brightness_spin, *object_brightness_sigma_spin;
GtkWidget *object_color_spin, *object_color_sigma_spin;
GtkWidget *object_center_x_spin, *object_center_y_spin;
GtkWidget *object_density_combo;
GtkWidget *sample_distribution_combo;
GtkObject *object_radius_adj;

GtkWidget *foreground_number_spin;
GtkWidget *foreground_brightness_spin, *foreground_brightness_sigma_spin;
GtkWidget *foreground_color_spin, *foreground_color_sigma_spin;

void reporter( const gchar *string, const gdouble fraction );

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
		{ GIMP_PDB_INT32, "psf", "Point spread function (0,1,2,3,4)" },
		{ GIMP_PDB_FLOAT, "sigma", "Sigma of psf gauss function [px]" },
		{ GIMP_PDB_FLOAT, "diffraction_percentage", "Diffraction lines of % stars" },
		{ GIMP_PDB_FLOAT, "diffraction_angle", "Angle of diffraction lines [deg]" },
		{ GIMP_PDB_FLOAT, "diffraction_length", "Length of diffraction lines [a.u.]" },
		{ GIMP_PDB_FLOAT, "diffraction_color", "Color interference of diffraction lines [a.u.]" },
		{ GIMP_PDB_INT32, "noise", "Noise in % of sqrt(N)" },
		{ GIMP_PDB_INT32, "background", "Background pixel value [absolute value]" },
		{ GIMP_PDB_FLOAT, "burnout", "Burn out % stars" },
		{ GIMP_PDB_FLOAT, "shininess", "Enlarge burnt out stars [a.u.]" },
		{ GIMP_PDB_INT32, "split_layers", "Split background, object and foreground stars to layers" },
		{ GIMP_PDB_INT32, "object_mask", "Create a mask using object stars" },
		{ GIMP_PDB_INT32, "random_seed", "Seed number" },
		{ GIMP_PDB_INT32, "random_seed_bool", "Use a random seed number" },
		{ GIMP_PDB_INT32, "background_stars", "Number of background stars" },
		{ GIMP_PDB_INT32, "object_stars", "Number of object stars" },
		{ GIMP_PDB_INT32, "foreground_stars", "Number of foreground stars" },
		{ GIMP_PDB_FLOAT, "background_brightness_mean", "Mean value of brightness of background stars [relative units]" },
		{ GIMP_PDB_FLOAT, "object_brightness_mean", "Mean value of brightness of object stars [relative units]" },
		{ GIMP_PDB_FLOAT, "foreground_brightness_mean", "Mean value of brightness of foreground stars [relative units]" },
		{ GIMP_PDB_FLOAT, "background_brightness_sigma", "Sigma of brightness of background stars" },
		{ GIMP_PDB_FLOAT, "object_brightness_sigma", "Sigma of brightness of object stars" },
		{ GIMP_PDB_FLOAT, "foreground_brightness_sigma", "Sigma of brightness of foreground stars" },
		{ GIMP_PDB_INT32, "background_color_mean", "Mean value of color of background stars [K]" },
		{ GIMP_PDB_INT32, "object_color_mean", "Mean value of color of object stars [K]" },
		{ GIMP_PDB_INT32, "foreground_color_mean", "Mean value of color of foreground stars [K]" },
		{ GIMP_PDB_INT32, "background_color_sigma", "Sigma of color of background stars" },
		{ GIMP_PDB_INT32, "object_color_sigma", "Sigma of color of object stars" },
		{ GIMP_PDB_INT32, "foreground_color_sigma", "Sigma of color of foreground stars" },
		{ GIMP_PDB_INT32, "object_x", "Center of object star distribution - x coordinate [px]" },
		{ GIMP_PDB_INT32, "object_y", "Center of object star distribution - y coordinate [px]" },
		{ GIMP_PDB_INT32, "object_density", "Object distribution (0,1,2)" },
		{ GIMP_PDB_INT32, "object_radius", "Radius of object distribution" },
		{ GIMP_PDB_INT32, "show_preview", "Show preview" }
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

	gimp_install_procedure(PLUG_IN_NAME,
		_("Create an artificial star distribution"),
		_("This plug-in creates an artificial star distribution. "),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		"2008, 2012",
		_("Artificial Stars"),
		"RGB*",
		GIMP_PLUGIN,
		nparams,
		nreturn_vals,
		params,
		return_vals );

	gimp_plugin_menu_register(PLUG_IN_NAME, _("<Image>/Filters/Astronomy") );
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

	switch( run_mode )
	{
		case GIMP_RUN_INTERACTIVE:
			gimp_get_data( PLUG_IN_NAME, &parameters );

			if ( !dialog( param[1].data.d_image, gimp_drawable_get( param[2].data.d_drawable ) ) )
			{
				return;
			}

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 37 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				image_id = param[1].data.d_image;

				parameters.psf = param[3].data.d_int32;
				parameters.sigma_psf = param[4].data.d_float;
				parameters.diffraction_percentage = param[5].data.d_float;
				parameters.diffraction_angle = param[6].data.d_float;
				parameters.diffraction_length = param[7].data.d_float;
				parameters.diffraction_color = param[8].data.d_float;
				parameters.noise = param[9].data.d_int32;
				parameters.background = param[10].data.d_int32;
				parameters.burnout = param[11].data.d_float;
				parameters.shininess = param[12].data.d_float;
				parameters.split_layers = param[13].data.d_int32;
				parameters.object_mask = param[14].data.d_int32;
				parameters.random_seed = param[15].data.d_int32;
				parameters.random_seed_bool = param[16].data.d_int32;
				parameters.background_stars = param[17].data.d_int32;
				parameters.object_stars = param[18].data.d_int32;
				parameters.foreground_stars = param[19].data.d_int32;
				parameters.background_brightness_mean = param[20].data.d_float;
				parameters.object_brightness_mean = param[21].data.d_float;
				parameters.foreground_brightness_mean = param[22].data.d_float;
				parameters.background_brightness_sigma = param[23].data.d_float;
				parameters.object_brightness_sigma = param[24].data.d_float;
				parameters.foreground_brightness_sigma = param[25].data.d_float;
				parameters.background_color_mean = param[26].data.d_int32;
				parameters.object_color_mean = param[27].data.d_int32;
				parameters.foreground_color_mean = param[28].data.d_int32;
				parameters.background_color_sigma = param[29].data.d_int32;
				parameters.object_color_sigma = param[30].data.d_int32;
				parameters.foreground_color_sigma = param[31].data.d_int32;
				parameters.object_x = param[32].data.d_int32;
				parameters.object_y = param[33].data.d_int32;
				parameters.object_density = param[34].data.d_int32;
				parameters.object_radius = param[35].data.d_int32;
				parameters.show_preview = param[36].data.d_int32;
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
		dialog_closed = TRUE;
		if ( distribution_dirty )
		{

			GtkWidget *message_dialog = gtk_message_dialog_new( NULL,
				GTK_DIALOG_DESTROY_WITH_PARENT, GTK_MESSAGE_QUESTION, GTK_BUTTONS_YES_NO,
				_( "Do you want to recalculate the star distribution before actually drawing the stars?" ) );
			gint response = gtk_dialog_run( GTK_DIALOG( message_dialog ) );
			gtk_widget_destroy( message_dialog );

			if ( response == GTK_RESPONSE_YES ) create_star_distribution();
		}
		create_stars();
	}

	values[0].data.d_status = status;
}

void star_new( gdouble **list, gint32 *list_number, const gint32 list_number_new )
{
	*list = malloc( list_number_new * sizeof(gdouble) * 4 );

	*list_number = list_number_new;
}

void star_free( gdouble **list, gint32 *list_number )
{
	if ( !*list ) return;

	if ( *list_number == 0 ) return;

	free( *list );
	*list = NULL;
	*list_number = 0;
}

void star_add( gdouble **list, const gint32 list_number, const gint32 list_position,
	const gdouble x, const gdouble y, const gdouble brightness, const gdouble color )
{
	if ( list_position >= list_number ) return;

	(*list)[list_position*4+0] = x;
	(*list)[list_position*4+1] = y;
	(*list)[list_position*4+2] = brightness;
	(*list)[list_position*4+3] = color;
}

void star_get( gdouble **list, const gint32 list_position,
	gdouble *x, gdouble *y, gdouble *brightness, gdouble *color )
{
	if ( x ) *x = (*list)[list_position*4+0];
	if ( y ) *y = (*list)[list_position*4+1];
	if ( brightness ) *brightness = (*list)[list_position*4+2];
	if ( color ) *color = (*list)[list_position*4+3];
}

inline gdouble gauss( const gdouble x, const gdouble y, const gdouble A, const gdouble x0, const gdouble y0,
	const gdouble sigma_x, const gdouble sigma_y )
{
	/* z = A * exp(-0.5*(x-x0)^2/sigma_x^2) * exp(-0.5*(y-y0)^2/sigma_y^2) + b */
	return A * exp( -0.5*(x-x0)*(x-x0)/(sigma_x*sigma_x) ) * exp( -0.5*(y-y0)*(y-y0)/(sigma_y*sigma_y) );
}

void create_star_distribution()
{
	star_free( &background_stars, &background_stars_number );
	star_free( &object_stars, &object_stars_number );
	star_free( &foreground_stars, &foreground_stars_number );

	star_new( &background_stars, &background_stars_number, parameters.background_stars );
	star_new( &object_stars, &object_stars_number, parameters.object_stars );
	star_new( &foreground_stars, &foreground_stars_number, parameters.foreground_stars );

	if ( !cancel_recalculation )
	{
		reporter( _("Creating background star distribution"), 0. );
		create_background_star_distribution();
	}

	if ( !cancel_recalculation )
	{
		reporter( _("Creating object star distribution"), 0. );
		create_object_star_distribution();
	}

	if ( !cancel_recalculation )
	{
		reporter( _("Creating foreground star distribution"), 0. );
		create_foreground_star_distribution();
	}

	if ( !cancel_recalculation )
	{
		find_star_norm_sort();

		recalculation_done();

		recalculating = FALSE;
		cancel_recalculation = TRUE;
		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Recalculate distribution") );
		gimp_preview_invalidate( GIMP_PREVIEW( preview ) );
	}

	reporter( "", 0. );
}

int comp( const void *ptr1, const void *ptr2 )
{
	const gdouble *p_1 = (gdouble*)ptr1;
	const gdouble *p_2 = (gdouble*)ptr2;

	if ( ( p_2[2] - p_1[2] ) > 0. ) return 1;
	if ( ( p_2[2] - p_1[2] ) < 0. ) return -1;
	return 0;
}

void find_star_norm_sort()
{
	reporter( _("Sorting stars by brightness"), 0. );

	star_free( &all_stars, &all_stars_number );

	/* Determine the normalize value */
	star_new( &all_stars, &all_stars_number, background_stars_number + object_stars_number + foreground_stars_number );
	if ( all_stars_number == 0 ) return;

	gdouble x, y, brightness, color;

	gint32 i;
	for ( i=0; i<background_stars_number; i++ )
	{
		star_get( &background_stars, i, &x, &y, &brightness, &color );

		star_add( &all_stars, all_stars_number, i, x, y, brightness, color );
	}
	for ( i=0; i<object_stars_number; i++ )
	{
		star_get( &object_stars, i, &x, &y, &brightness, &color );

		star_add( &all_stars, all_stars_number, i+background_stars_number, x, y, brightness, color );
	}
	for ( i=0; i<foreground_stars_number; i++ )
	{
		star_get( &foreground_stars, i, &x, &y, &brightness, &color );

		star_add( &all_stars, all_stars_number, i+background_stars_number+object_stars_number, x, y, brightness, color );
	}

	qsort( all_stars, all_stars_number, 4 * sizeof( gdouble ), comp );

	find_star_norm();

	reporter( "", 0. );
}

void find_star_norm()
{
	if ( all_stars_number == 0 ) return;

	gdouble brightness;

	star_get( &all_stars, (gint32)((parameters.burnout*all_stars_number)/100),
		NULL, NULL, &brightness, NULL );

	brightness_star_norm = brightness;

	star_get( &all_stars, (gint32)((parameters.diffraction_percentage*all_stars_number)/100),
		NULL, NULL, &brightness, NULL );

	brightness_star_diffraction = brightness;
}

void create_background_star_distribution()
{
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed/*time( NULL )*/ );

	gdouble x, y, brightness, color;

	gint i;
	for ( i=0; i<parameters.background_stars; i++ )
	{
		if ( cancel_recalculation ) break;

		x = gimp_image_width( image_id )*gsl_rng_uniform( r );
		y = gimp_image_height( image_id )*gsl_rng_uniform( r );

		brightness = parameters.background_brightness_mean + gsl_ran_gaussian( r, parameters.background_brightness_sigma );
		color = parameters.background_color_mean + gsl_ran_gaussian( r, parameters.background_color_sigma );

		star_add( &background_stars, background_stars_number, i, x, y, brightness, color );

		if ( i%(1+parameters.background_stars/100) ) reporter( "", (gdouble)(i+1)/parameters.background_stars );
	}

	gsl_rng_free( r );
}

void create_object_star_distribution()
{
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed+10/*time( NULL )*/ );

	gdouble x, y, brightness, color;
	gdouble radius, angle1, angle2;

	gint i;

	switch( parameters.object_density )
	{
		case DENSITY_RANDOM:
		{
			for ( i=0; i<parameters.object_stars; i++ )
			{
				if ( cancel_recalculation ) break;

				x = gimp_image_width( image_id )*gsl_rng_uniform( r );
				y = gimp_image_height( image_id )*gsl_rng_uniform( r );

				brightness = -1.;
				while ( brightness < 0. ) brightness = parameters.object_brightness_mean + gsl_ran_gaussian( r, parameters.object_brightness_sigma );
				color = parameters.object_color_mean + gsl_ran_gaussian( r, parameters.object_color_sigma );

				star_add( &object_stars, object_stars_number, i, x, y, brightness, color );

				if ( i%(1+parameters.object_stars/100) ) reporter( "", (gdouble)(i+1)/parameters.object_stars );
			}

			break;
		}
		case DENSITY_GAUSS:
		{
			for ( i=0; i<parameters.object_stars; i++ )
			{
				if ( cancel_recalculation ) break;

				radius = abs( gsl_ran_gaussian( r, parameters.object_radius * GSL_MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) / 100 ) ) );
				angle1 = gsl_rng_uniform( r ) * 2 * M_PI;

				x = parameters.object_x + radius * cos( angle1 );
				y = parameters.object_y + radius * sin( angle1 );

				brightness = -1.;
				while ( brightness < 0. ) brightness = parameters.object_brightness_mean + gsl_ran_gaussian( r, parameters.object_brightness_sigma );
				color = parameters.object_color_mean + gsl_ran_gaussian( r, parameters.object_color_sigma );

				star_add( &object_stars, object_stars_number, i, x, y, brightness, color );

				if ( i%(1+parameters.object_stars/100) ) reporter( "", (gdouble)(i+1)/parameters.object_stars );
			}

			break;
		}
		case DENSITY_PLUMMER:
		{
			for ( i=0; i<parameters.object_stars; i++ )
			{
				if ( cancel_recalculation ) break;

				do
				{
					/* 20: arbitrary radius, where P(random*radius)->0, here P(random*radius)<0.1% */
					radius = gsl_rng_uniform( r ) * 20;
				}
				while ( 5.38 * radius * radius / pow( 1. + radius * radius, 2.5 ) < gsl_rng_uniform( r ) );

				angle1 = gsl_rng_uniform( r ) * 2 * M_PI;

				do
				{
					angle2 = gsl_rng_uniform( r ) * M_PI;
				}
				while ( sin( angle2 ) <= gsl_rng_uniform( r ) );

				x = parameters.object_x + radius * ( (gdouble)( parameters.object_radius * GSL_MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ) / ( 100 * 20 ) ) * sin( angle2 ) * cos( angle1 ) ;
				y = parameters.object_y + radius * ( (gdouble)( parameters.object_radius * GSL_MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ) / ( 100 * 20 ) ) * sin( angle2 ) * sin( angle1 );

				brightness = -1.;
				while ( brightness < 0. ) brightness = parameters.object_brightness_mean + gsl_ran_gaussian( r, parameters.object_brightness_sigma );
				color = parameters.object_color_mean + gsl_ran_gaussian( r, parameters.object_color_sigma );

				star_add( &object_stars, object_stars_number, i, x, y, brightness, color );

				if ( i%(1+parameters.object_stars/100) ) reporter( "", (gdouble)(i+1)/parameters.object_stars );
			}

			break;
		}
		default:
			break;
	}

	gsl_rng_free( r );
}

void create_foreground_star_distribution()
{
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed+20/*time( NULL )*/ );

	gdouble x, y, brightness, color;

	gint i;
	for ( i=0; i<parameters.foreground_stars; i++ )
	{
		if ( cancel_recalculation ) break;

		x = gimp_image_width( image_id )*gsl_rng_uniform( r );
		y = gimp_image_height( image_id )*gsl_rng_uniform( r );

		brightness = parameters.foreground_brightness_mean + gsl_ran_gaussian( r, parameters.foreground_brightness_sigma );
		color = parameters.foreground_color_mean + gsl_ran_gaussian( r, parameters.foreground_color_sigma );

		star_add( &foreground_stars, foreground_stars_number, i, x, y, brightness, color );

		if ( i%(1+parameters.foreground_stars/100) ) reporter( "", (gdouble)(i+1)/parameters.foreground_stars );
	}

	gsl_rng_free( r );
}

void create_stars_psf( gdouble *star_distribution, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, gdouble **list, const gint32 list_number, const gboolean show_progress )
{
	gdouble x, y;
	gdouble brightness, color;
	gdouble r, g, b;

	gint32 i, x_it, y_it;

	if ( parameters.psf == PSF_DELTA_PEAK )
	{
		for ( i=0; i<list_number; i++ )
		{
			star_get( list, i, &x, &y, &brightness, &color );

			if ( (gint32)x >= x_start &&
				(gint32)y >= y_start &&
				(gint32)x < x_start+x_size &&
				(gint32)y < y_start+y_size )
			{
				temperature_to_rgb_relative( color, &r, &g, &b );

				star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+0] =
					GSL_MAX(
						brightness * r,
						star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+0]
					);

				star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+1] =
					GSL_MAX(
						brightness * g,
						star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+1]
					);

				star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+2] =
					GSL_MAX(
						brightness * b,
						star_distribution[3*((gint32)(y-y_start)*x_size+(gint32)(x-x_start))+2]
					);
			}

			if ( show_progress && i % (1+list_number/100) ) gimp_progress_update( (gdouble)(i) / list_number );
		}
	}
	else if ( parameters.psf == PSF_GAUSS || parameters.psf == PSF_GAUSS_4_DIFFRACTIONS ||
		parameters.psf == PSF_GAUSS_6_DIFFRACTIONS || parameters.psf == PSF_GAUSS_12_DIFFRACTIONS )
	{
		gdouble draw_width;
		gdouble sigma_new;
		gdouble brightness_new;

		for ( i=0; i<list_number; i++ )
		{
			star_get( list, i, &x, &y, &brightness, &color );
			temperature_to_rgb_relative( color, &r, &g, &b );

			sigma_new = parameters.sigma_psf;
			brightness_new = brightness;
			if ( brightness - brightness_star_norm > 0. )
			{
				sigma_new += 2 * sigma_new * parameters.shininess * ( brightness - brightness_star_norm ) / brightness_star_norm;
				brightness_new -= brightness_new * 16 * pow( parameters.shininess, 0.2 ) * ( brightness - brightness_star_norm ) / ( brightness_star_norm * brightness_star_norm );
			}
			draw_width = 4*sigma_new;

			if ( (gint32)(x+draw_width) >= x_start &&
				(gint32)(y+draw_width) >= y_start &&
				(gint32)(x-draw_width) < x_start+x_size &&
				(gint32)(y-draw_width) < y_start+y_size )
			{
				for ( x_it=(gint32)(x-draw_width); x_it<=(gint32)(x+draw_width); x_it++ )
				{
					for ( y_it=(gint32)(y-draw_width); y_it<=(gint32)(y+draw_width); y_it++ )
					{
						if ( x_it >= x_start &&
							y_it >= y_start &&
							x_it < x_start+x_size &&
							y_it < y_start+y_size )
						{
							star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+0] =
								GSL_MAX(
									r * gauss( x_it, y_it, brightness, x, y, sigma_new, sigma_new ),
									star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+0]
								);

							star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+1] =
								GSL_MAX(
									g * gauss( x_it, y_it, brightness, x, y, sigma_new, sigma_new ),
									star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+1]
								);

							star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+2] =
								GSL_MAX(
									b * gauss( x_it, y_it, brightness, x, y, sigma_new, sigma_new ),
									star_distribution[3*((y_it-y_start)*x_size+(x_it-x_start))+2]
								);
						}
					}
				}
			}

			if ( parameters.psf == PSF_GAUSS_4_DIFFRACTIONS ||
				parameters.psf == PSF_GAUSS_6_DIFFRACTIONS ||
				parameters.psf == PSF_GAUSS_12_DIFFRACTIONS )
			{

				gint32 radius_it, angle_it, offset;
				gdouble x_pos, y_pos, radius, x_offset, y_offset;
				gdouble brightness_new_new;
				gdouble radius_half;

				const gint32 oversample = 5;

				gdouble angle_start = M_PI*parameters.diffraction_angle/180;

				const gint32 number_of_diffractions = ( parameters.psf == PSF_GAUSS_4_DIFFRACTIONS ) ? 4 :
					( parameters.psf == PSF_GAUSS_6_DIFFRACTIONS ) ? 6 : 12;
				const gdouble multiplier = ( parameters.psf == PSF_GAUSS_4_DIFFRACTIONS ) ? 1. :
					( parameters.psf == PSF_GAUSS_6_DIFFRACTIONS ) ? 0.6 : 0.3;


				if ( brightness - brightness_star_diffraction > 0. )
				{
					radius_half = -1.;

					for ( radius_it=1; radius_it<multiplier*parameters.diffraction_length*oversample*24*sigma_new; radius_it++ )
					{
						for ( angle_it=0; angle_it<number_of_diffractions; angle_it++ )
						{
							for ( offset=floor(-sigma_new*oversample); offset<=ceil(sigma_new*oversample); offset++ )
							{
								x_pos = radius_it*cos(angle_it*2.*M_PI/number_of_diffractions+angle_start)/oversample;
								y_pos = radius_it*sin(angle_it*2.*M_PI/number_of_diffractions+angle_start)/oversample;

								x_offset = offset*cos(angle_it*2.*M_PI/number_of_diffractions+angle_start-M_PI/2)/oversample;
								y_offset = offset*sin(angle_it*2.*M_PI/number_of_diffractions+angle_start-M_PI/2)/oversample;

								if ( (gint32)((x_pos+x_offset)+x) >= x_start &&
									(gint32)((y_pos+y_offset)+y) >= y_start &&
									(gint32)((x_pos+x_offset)+x) < x_start+x_size &&
									(gint32)((y_pos+y_offset)+y) < y_start+y_size )
								{

/*
gnuplot> param=1.0
gnuplot> half=130
gnuplot> plot [1:600][0:255] "artificial_star_spike" using 1:2 smooth csplines, "artificial_star_spike" using 1:3 smooth csplines, "artificial_star_spike" using 1:4 smooth csplines, (450./sqrt(param))*exp(-x/(3*27./sqrt(param)))+28, (450./sqrt(param))*exp(-x/(3*27./sqrt(param)))+0.5*(1./sqrt(param))*27*half*(450./sqrt(param))*sin(1.3*sqrt(param)*(x-half)/27+pi)/((x+1)**2)+28
*/

									radius = sqrt( (x_pos+x_offset)*(x_pos+x_offset) + (y_pos+y_offset)*(y_pos+y_offset) );

									brightness_new_new = CLAMP(r*(brightness/sqrt(1./r))*exp(-radius/(multiplier*3*parameters.diffraction_length*sigma_new/sqrt(1./r))),0.,r*brightness_star_norm);

									if ( brightness_new_new < 0.5*brightness_star_norm )
									{
										if ( radius_half < 0. ) radius_half = radius;

										brightness_new_new +=
											0.5*(1.-exp(-radius/(8*sigma_new)))*parameters.diffraction_color*brightness_new_new*sin(0.25*(radius-radius_half)/parameters.sigma_psf+1.0*M_PI);
									}

									brightness_new_new *=
										0.6*exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(0.1*parameters.sigma_psf*parameters.sigma_psf)) + 
										0.4*exp((fabs(offset)/(10*oversample))/(4*parameters.sigma_psf*parameters.sigma_psf))*
										exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(parameters.sigma_psf*parameters.sigma_psf));

									star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+0] =
										GSL_MAX(
											brightness_new_new,
											star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+0]
										);

									brightness_new_new = CLAMP(g*(brightness/sqrt(1./g))*exp(-radius/(multiplier*3*parameters.diffraction_length*sigma_new/sqrt(1./g))),0.,g*brightness_star_norm);

									if ( brightness_new_new < 0.5*brightness_star_norm )
									{
										if ( radius_half < 0. ) radius_half = radius;

										brightness_new_new += 
											0.5*(1.-exp(-radius/(8*sigma_new)))*parameters.diffraction_color*brightness_new_new*sin(0.25*(radius-radius_half)/parameters.sigma_psf+1.5*M_PI);
									}

									brightness_new_new *=
										0.6*exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(0.1*parameters.sigma_psf*parameters.sigma_psf)) + 
										0.4*exp((fabs(offset)/(10*oversample))/(4*parameters.sigma_psf*parameters.sigma_psf))*
										exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(parameters.sigma_psf*parameters.sigma_psf));

									star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+1] =
										GSL_MAX(
											brightness_new_new,
											star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+1]
										);

									brightness_new_new = CLAMP(b*(brightness/sqrt(1./b))*exp(-radius/(multiplier*3*parameters.diffraction_length*sigma_new/sqrt(1./b))),0.,b*brightness_star_norm);

									if ( brightness_new_new < 0.5*brightness_star_norm )
									{
										if ( radius_half < 0. ) radius_half = radius;

											brightness_new_new +=
												0.5*(1.-exp(-radius/(8*sigma_new)))*parameters.diffraction_color*brightness_new_new*sin(0.25*(radius-radius_half)/parameters.sigma_psf+2.0*M_PI);
									}

									brightness_new_new *=
										0.6*exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(0.1*parameters.sigma_psf*parameters.sigma_psf)) + 
										0.4*exp((fabs(offset)/(10*oversample))/(4*parameters.sigma_psf*parameters.sigma_psf))*
										exp(-0.5*((gdouble)(offset*offset)/(oversample*oversample))/(parameters.sigma_psf*parameters.sigma_psf));

									star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+2] =
										GSL_MAX(
											brightness_new_new,
											star_distribution[3*((gint32)((y_pos+y_offset)+y-y_start)*x_size+(gint32)((x_pos+x_offset)+x-x_start))+2]
										);
								}
							}
						}
					}
				}
			}

			if ( show_progress && i % (1+list_number/100) ) gimp_progress_update( (gdouble)(i) / list_number );
		}
	}
}

void create_stars_noise( gdouble *star_distribution, gint32 x_size, gint32 y_size, const gboolean show_progress )
{
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed+30/*time( NULL )*/ );

	gint32 i;
	for ( i=0; i<3*x_size*y_size; i++ )
	{
		if ( 255*star_distribution[i]/brightness_star_norm < 255 )
			star_distribution[i] += gsl_ran_gaussian( r, parameters.noise*sqrt(star_distribution[i])/100 );

		if ( show_progress && i % (1+3*x_size*y_size/100) )
			gimp_progress_update( (gdouble)(i) / 3*x_size*y_size );
	}

	gsl_rng_free( r );
}

/*
Create stars into destination
*/
static void create_stars_selection( GimpPixelRgn *region_destination,
	const gboolean draw_background_stars, const gboolean draw_object_stars, const gboolean draw_foreground_stars,
	const gboolean show_progress )
{
	gint32 x, y;
	gint32 i;

	/* create the stars before actually drawing them */
	gdouble *star_distribution;
	star_distribution = g_malloc( sizeof(gdouble)*3*region_destination->w*region_destination->h );

	reporter( "artificial_stars: Creating background...", -1.0 );
	if ( show_progress ) gimp_progress_init( _("Creating background...") );
	/* Background first */
	if ( draw_background_stars )
	{
		for ( i=0; i<3*region_destination->w*region_destination->h; i++ )
		{
			star_distribution[i] = (gdouble)parameters.background * brightness_star_norm / 255;

			if ( show_progress && i % (1+3*region_destination->w*region_destination->h/100) )
				gimp_progress_update( (gdouble)(i) / 3*region_destination->w*region_destination->h );
		}
	}
	else
	{
		for( i=0; i<3*region_destination->w*region_destination->h; i++ )
		{
			star_distribution[i] = 0.;

			if ( show_progress && i % (1+3*region_destination->w*region_destination->h/100) )
				gimp_progress_update( (gdouble)(i) / 3*region_destination->w*region_destination->h );
		}
	}

	/* Background stars */
	if ( draw_background_stars )
	{
		reporter( "artificial_stars: Creating background stars...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating background stars...") );

		create_stars_psf( star_distribution, region_destination->x, region_destination->y,
			region_destination->w, region_destination->h, &background_stars, background_stars_number, show_progress );
	}

	/* Object stars */
	if ( draw_object_stars )
	{
		reporter( "artificial_stars: Creating object stars...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating object stars...") );

		create_stars_psf( star_distribution, region_destination->x, region_destination->y,
			region_destination->w, region_destination->h, &object_stars, object_stars_number, show_progress );
	}

	/* Foreground stars */
	if ( draw_foreground_stars )
	{
		reporter( "artificial_stars: Creating foreground stars...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating foreground stars...") );

		create_stars_psf( star_distribution, region_destination->x, region_destination->y,
			region_destination->w, region_destination->h, &foreground_stars, foreground_stars_number, show_progress );
	}

	/* Noise */
	if ( parameters.noise > 0 )
	{
		reporter( "artificial_stars: Creating noise...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating noise...") );

		create_stars_noise( star_distribution, region_destination->w, region_destination->h, show_progress );
	}

	reporter( "artificial_stars: Writing normalized values to buffer...", -1.0 );
	if ( show_progress ) gimp_progress_init( _("Writing normalized values to buffer...") );

	guchar *buffer;
	buffer = g_malloc( sizeof(guchar)*region_destination->bpp*region_destination->w*region_destination->h );
	for ( y=0; y<region_destination->h; y++ )
	{
		for ( x=0; x<region_destination->w; x++ )
		{
			for ( i=0; i<3; i++ )
			{
				buffer[region_destination->bpp*(y*region_destination->w+x)+i] =
					CLAMP( 255 * star_distribution[3*(y*region_destination->w+x)+i] / brightness_star_norm, 0, 255 );
			}
			if ( region_destination->bpp > 3 )
				buffer[region_destination->bpp*(y*region_destination->w+x)+region_destination->bpp-1] = 255;
		}

		if ( show_progress && y % (1+region_destination->h/100) )
			gimp_progress_update( (gdouble)(y) / region_destination->h );
	}

	reporter( "artificial_stars: Setting buffer to image...", -1.0 );
	gimp_pixel_rgn_set_rect( region_destination, buffer, region_destination->x, region_destination->y,
		region_destination->w, region_destination->h );

	if ( show_progress ) gimp_progress_end();

	g_free( star_distribution );
	g_free( buffer );
}

/*
create stars
*/
static void create_stars()
{
	gimp_image_undo_group_start( image_id );

	gint32 layer_destination;
	GimpPixelRgn region_destination;
	if ( parameters.split_layers )
	{
		layer_destination = gimp_layer_new( image_id, _("Artificial background stars"), gimp_image_width( image_id ),
			gimp_image_height( image_id ), GIMP_RGB_IMAGE, 100, GIMP_NORMAL_MODE );

		gimp_image_add_layer( image_id, layer_destination, 0 );

		gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
			gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

		create_stars_selection( &region_destination, TRUE, FALSE, FALSE, TRUE );

		gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
		gimp_drawable_merge_shadow( layer_destination, TRUE );
		gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );

		layer_destination = gimp_layer_new( image_id, _("Artificial object stars"), gimp_image_width( image_id ),
			gimp_image_height( image_id ), GIMP_RGB_IMAGE, 100, GIMP_LIGHTEN_ONLY_MODE );

		gimp_image_add_layer( image_id, layer_destination, 0 );

		gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
			gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

		create_stars_selection( &region_destination, FALSE, TRUE, FALSE, TRUE );

		gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
		gimp_drawable_merge_shadow( layer_destination, TRUE );
		gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );

		layer_destination = gimp_layer_new( image_id, _("Artificial foreground stars"), gimp_image_width( image_id ),
			gimp_image_height( image_id ), GIMP_RGB_IMAGE, 100, GIMP_LIGHTEN_ONLY_MODE );

		gimp_image_add_layer( image_id, layer_destination, 0 );

		gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
			gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

		create_stars_selection( &region_destination, FALSE, FALSE, TRUE, TRUE );

		gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
		gimp_drawable_merge_shadow( layer_destination, TRUE );
		gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );
	}
	else
	{
		layer_destination = gimp_layer_new( image_id, _("Artificial stars"), gimp_image_width( image_id ),
			gimp_image_height( image_id ), GIMP_RGB_IMAGE, 100, GIMP_NORMAL_MODE );

		gimp_image_add_layer( image_id, layer_destination, 0 );

		gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
			gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

		create_stars_selection( &region_destination, TRUE, TRUE, TRUE, TRUE );

		gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
		gimp_drawable_merge_shadow( layer_destination, TRUE );
		gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
		gimp_drawable_height( layer_destination ) );
	}

	star_free( &background_stars, &background_stars_number );
	star_free( &object_stars, &object_stars_number );
	star_free( &foreground_stars, &foreground_stars_number );

	gimp_image_undo_group_end( image_id );

	gimp_displays_flush();
}


/*
GUI
*/

/* Preview callback */
static void preview_callback( GimpPreview *preview, gpointer *data )
{
	GimpDrawable *drawable = gimp_drawable_preview_get_drawable( GIMP_DRAWABLE_PREVIEW( preview ) );

	gint32 x_pos, y_pos, x_size, y_size;
	gimp_preview_get_position( preview, &x_pos, &y_pos );
	gimp_preview_get_size( preview, &x_size, &y_size );

	GimpPixelRgn region_destination;
	gimp_pixel_rgn_init( &region_destination, drawable, x_pos, y_pos, x_size, y_size, TRUE, TRUE );

	create_stars_selection( &region_destination, TRUE, TRUE, TRUE, FALSE );

	gimp_drawable_flush( drawable );
	gimp_drawable_update( drawable->drawable_id, 0, 0, drawable->width, drawable->height );

	/* Read from shadow tiles ("FALSE, TRUE") */
	gimp_pixel_rgn_init( &region_destination, drawable, 0, 0, drawable->width, drawable->height, FALSE, TRUE );
	gimp_drawable_preview_draw_region( GIMP_DRAWABLE_PREVIEW( preview ), &region_destination );
}

static void spin_background_stars( GtkWidget *spin )
{
	parameters.background_stars = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_background_brightness_mean( GtkWidget *spin )
{
	parameters.background_brightness_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_background_brightness_sigma( GtkWidget *spin )
{
	parameters.background_brightness_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_background_color_mean( GtkWidget *spin )
{
	parameters.background_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_background_color_sigma( GtkWidget *spin )
{
	parameters.background_color_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_stars( GtkWidget *spin )
{
	parameters.object_stars = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_brightness_mean( GtkWidget *spin )
{
	parameters.object_brightness_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_brightness_sigma( GtkWidget *spin )
{
	parameters.object_brightness_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_color_mean( GtkWidget *spin )
{
	parameters.object_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_color_sigma( GtkWidget *spin )
{
	parameters.object_color_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_x( GtkWidget *spin )
{
	parameters.object_x = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_y( GtkWidget *spin )
{
	parameters.object_y = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_foreground_stars( GtkWidget *spin )
{
	parameters.foreground_stars = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_foreground_brightness_mean( GtkWidget *spin )
{
	parameters.foreground_brightness_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_foreground_brightness_sigma( GtkWidget *spin )
{
	parameters.foreground_brightness_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_foreground_color_mean( GtkWidget *spin )
{
	parameters.foreground_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_foreground_color_sigma( GtkWidget *spin )
{
	parameters.foreground_color_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_random_seed( GtkWidget *spin )
{
	parameters.random_seed = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void split_layers_toggled( GtkWidget *check )
{
	parameters.split_layers = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void object_mask_toggled( GtkWidget *check )
{
	parameters.object_mask = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void recalculate_distribution_clicked( GtkWidget *button )
{
	if ( !recalculating )
	{
		recalculating = TRUE;
		cancel_recalculation = FALSE;

		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Cancel recalculation") );

		create_star_distribution();
	}
	else
	{
		recalculating = FALSE;
		cancel_recalculation = TRUE;

		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Recalculate distribution") );
	}
}

static void recalculation_necessary()
{
	distribution_dirty = TRUE;

	GdkColor color;
	color.red = 55000;
	color.green = 12000;
	color.blue = 15000;
	gtk_widget_modify_bg( recalculate_button, GTK_STATE_NORMAL, &color );
	gtk_widget_modify_bg( recalculate_button, GTK_STATE_PRELIGHT, &color );
	gtk_widget_modify_bg( recalculate_button, GTK_STATE_ACTIVE, &color );
}

static void recalculation_done()
{
	distribution_dirty = FALSE;

	gtk_widget_modify_bg( recalculate_button, GTK_STATE_NORMAL, NULL );
	gtk_widget_modify_bg( recalculate_button, GTK_STATE_PRELIGHT, NULL );
	gtk_widget_modify_bg( recalculate_button, GTK_STATE_ACTIVE, NULL );
}

static void sample_distribution_clicked( GtkWidget *button )
{
	gdouble number_percentage;
	switch ( sample_distribution )
	{
		case STANDARD_HD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)40000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 0. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 15. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)5000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 45. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 30. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 7000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 1000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_RANDOM );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)2000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 80. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 35. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		case STANDARD_LD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)2500*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 15. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 10. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)500*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 45. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 30. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 7000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 1000 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_RANDOM );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)15*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 90. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 40. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		case GLOBULAR_HD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)5000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 0. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 45. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)50000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 100. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 40. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 6500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 800 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_PLUMMER );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( GIMP_SCALE_ENTRY_SPINBUTTON( object_radius_adj ) ), 50 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)1000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 90. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 60. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		case GLOBULAR_LD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)5000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 0. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 45. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)10000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 100. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 40. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 6500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 800 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_PLUMMER );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( GIMP_SCALE_ENTRY_SPINBUTTON( object_radius_adj ) ), 15 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)1000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 90. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 60. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		case OPEN_HD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)40000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 0. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 25. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7200 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)1200*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 80. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 35. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 7100 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_GAUSS );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( GIMP_SCALE_ENTRY_SPINBUTTON( object_radius_adj ) ), 6 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)2000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 100. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 40. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7200 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		case OPEN_LD:
			number_percentage = (gdouble)(gimp_image_width( image_id )*gimp_image_height( image_id ))/3.3e6;

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), (gint32)40000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 0. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 25. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7200 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), (gint32)300*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 110. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 30. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 7100 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 500 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
			gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_PLUMMER );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( GIMP_SCALE_ENTRY_SPINBUTTON( object_radius_adj ) ), 180 );

			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), (gint32)2000*number_percentage );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 100. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 40. );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7200 );
			gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
			break;
		default:
			break;
	}
  

}

void reporter( const gchar *string, const gdouble fraction )
{
	if ( strcmp( string, "" ) == 0 )
	{
		if ( fraction == 0. )
		{
			if ( dialog_closed )
			{
				gimp_progress_init( string );
			}
			else
			{
				gtk_progress_bar_set_text( GTK_PROGRESS_BAR( progress_bar ), string );
				gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
				gtk_main_iteration_do( FALSE );
			}
		}
		else if ( fraction > 0. && fraction <= 1. )
		{
			if ( dialog_closed )
			{
				gimp_progress_update( fraction );
			}
			else
			{
				gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
				gtk_main_iteration_do( FALSE );
			}
		}
	}
	else
	{
		if ( fraction >= 0. && fraction <= 1. )
		{
			if ( dialog_closed )
			{
				gimp_progress_init( string );
				gimp_progress_update( fraction );
			}
			else
			{
				gtk_progress_bar_set_text( GTK_PROGRESS_BAR( progress_bar ), string );
				gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
				gtk_main_iteration_do( FALSE );
			}
		}
		else
		{
			printf( "artificial_stars: %s\n", string );
		}
	}
}

/*
main dialog
*/
static gint dialog( gint32 image_id, GimpDrawable *drawable )
{
	GtkWidget *dlg;
	GtkWidget *notebook;
	GtkWidget *main_hbox;
	GtkWidget *left_vbox;
	GtkWidget *right_vbox;
	GtkWidget *label;
	GtkWidget *button;
	GtkWidget *object_box;
	GtkWidget *check_box;
	GtkWidget *random_seed;
	GtkWidget *frame;
	GtkWidget *table;
	GtkObject *adj;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( _("Artificial Stars"), "astro_artificial_stars", NULL, 0,
		gimp_standard_help_func, PLUG_IN_NAME,
		GTK_STOCK_CANCEL, GTK_RESPONSE_CANCEL,
		GTK_STOCK_OK, GTK_RESPONSE_OK,
		NULL );

/* General layout */
	main_hbox = gtk_hbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( main_hbox ), 8 );
	gtk_container_add( GTK_CONTAINER( GTK_DIALOG( dlg )->vbox ), main_hbox );

	left_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( left_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), left_vbox, FALSE, FALSE, 0 );

/*	GtkWidget *vertical_line = gtk_vseparator_new();
	gtk_box_pack_start( GTK_BOX( main_hbox ), vertical_line, FALSE, FALSE, 0 );
	gtk_widget_show( vertical_line ); */

	right_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( right_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), right_vbox, FALSE, FALSE, 0 );

/* Preview (left side) */

	preview = gimp_drawable_preview_new( drawable, &parameters.show_preview );
	gtk_box_pack_start( GTK_BOX( left_vbox ), preview, FALSE, FALSE, 0 );
	gtk_widget_show( preview );

	g_signal_connect( preview, "invalidated", G_CALLBACK( preview_callback ), NULL );

/* Star distribution (left side) */

	frame = gimp_frame_new( _("Star Distribution") );
	gtk_box_pack_start( GTK_BOX( left_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 3, 1, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	recalculate_button = gtk_button_new_with_label( _("Recalculate distribution") );
	gtk_table_attach( GTK_TABLE( table ), recalculate_button, 0, 1, 0, 1, GTK_FILL | GTK_EXPAND, 0, 0, 0 );

	recalculation_necessary();

	gtk_widget_show( recalculate_button );
	g_signal_connect( recalculate_button, "clicked", G_CALLBACK( recalculate_distribution_clicked ), NULL );

	progress_bar = gtk_progress_bar_new();
	gtk_widget_set_size_request( progress_bar, -1, 24 );
	gtk_table_attach( GTK_TABLE( table ), progress_bar, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( progress_bar );
	gtk_progress_bar_set_orientation( GTK_PROGRESS_BAR( progress_bar ), GTK_PROGRESS_LEFT_TO_RIGHT );

	random_seed = gimp_random_seed_new( &parameters.random_seed, &parameters.random_seed_bool );
	gtk_table_attach( GTK_TABLE( table ), random_seed, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( random_seed );
	g_signal_connect( GIMP_RANDOM_SEED_SPINBUTTON( random_seed ), "value_changed", G_CALLBACK( spin_random_seed ), NULL );
	g_signal_connect_swapped( GIMP_RANDOM_SEED_SPINBUTTON( random_seed ), "value_changed", G_CALLBACK( recalculation_necessary ), NULL );

/* Sample Distributions */

	frame = gimp_frame_new( _("Sample Distributions") );
	gtk_box_pack_start( GTK_BOX( left_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 1, 2, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	sample_distribution_combo = gimp_int_combo_box_new(
		_("Standard (high density)"), STANDARD_HD,
		_("Standard (low density)"), STANDARD_LD,
		_("Globular cluster (high density)"), GLOBULAR_HD,
		_("Globular cluster (low density)"), GLOBULAR_LD,
		_("Open cluster (high density)"), OPEN_HD,
		_("Open cluster (low density)"), OPEN_LD,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( sample_distribution_combo ), sample_distribution );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( sample_distribution_combo ), sample_distribution,
		G_CALLBACK( gimp_int_combo_box_get_active ), &sample_distribution );
	gtk_table_attach( GTK_TABLE( table ), sample_distribution_combo, 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( sample_distribution_combo );

	button = gtk_button_new_with_label( _("Apply") );
	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	g_signal_connect( button, "clicked", G_CALLBACK( sample_distribution_clicked ), NULL );
	g_signal_connect_swapped( button, "clicked", G_CALLBACK( recalculation_necessary ), NULL );


/*  Create notebook  */
	notebook = gtk_notebook_new ();
	gtk_box_pack_start (GTK_BOX (main_hbox), notebook, FALSE, FALSE, 0);
	gtk_widget_show (notebook);

/*  Distribution options page  */
	right_vbox = gtk_box_new (GTK_ORIENTATION_VERTICAL, 12);
	gtk_container_set_border_width (GTK_CONTAINER (right_vbox), 12);
	gtk_notebook_append_page (GTK_NOTEBOOK (notebook), right_vbox,
                            gtk_label_new_with_mnemonic (_("_Distribution Options")));
	gtk_widget_show (right_vbox);

	frame = gimp_frame_new( _("Background Stars") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 3, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("Amount:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	background_number_spin = gtk_spin_button_new_with_range( 0, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), parameters.background_stars );
	g_signal_connect( background_number_spin, "value_changed", G_CALLBACK( spin_background_stars ), NULL );
	g_signal_connect_swapped( background_number_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), background_number_spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( background_number_spin );

	label = gtk_label_new( _("Brightness mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	background_brightness_spin = gtk_spin_button_new_with_range( 0, 1000, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), parameters.background_brightness_mean );
	g_signal_connect( background_brightness_spin, "value_changed", G_CALLBACK( spin_background_brightness_mean ), NULL );
	g_signal_connect_swapped( background_brightness_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), background_brightness_spin, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( background_brightness_spin );

	label = gtk_label_new( _("Brightness sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	background_brightness_sigma_spin = gtk_spin_button_new_with_range( 0, 500, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), parameters.background_brightness_sigma );
	g_signal_connect( background_brightness_sigma_spin, "value_changed", G_CALLBACK( spin_background_brightness_sigma ), NULL );
	g_signal_connect_swapped( background_brightness_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), background_brightness_sigma_spin, 3, 4, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( background_brightness_sigma_spin );

	label = gtk_label_new( _("Color mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	background_color_spin = gtk_spin_button_new_with_range( 1000, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), parameters.background_color_mean );
	g_signal_connect( background_color_spin, "value_changed", G_CALLBACK( spin_background_color_mean ), NULL );
	g_signal_connect_swapped( background_color_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), background_color_spin, 1, 2, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( background_color_spin );

	label = gtk_label_new( _("Color sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	background_color_sigma_spin = gtk_spin_button_new_with_range( 1, 10000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), parameters.background_color_sigma );
	g_signal_connect( background_color_sigma_spin, "value_changed", G_CALLBACK( spin_background_color_sigma ), NULL );
	g_signal_connect_swapped( background_color_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), background_color_sigma_spin, 3, 4, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( background_color_sigma_spin );


	frame = gimp_frame_new( _("Object Stars") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	object_box = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( object_box ), 0 );
	gtk_container_add( GTK_CONTAINER( frame ), object_box );

	table = gtk_table_new( 3, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_box_pack_start( GTK_BOX( object_box ), table, FALSE, FALSE, 0 );
	gtk_widget_show( table );

	label = gtk_label_new( _("Amount:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_number_spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), parameters.object_stars );
	g_signal_connect( object_number_spin, "value_changed", G_CALLBACK( spin_object_stars ), NULL );
	g_signal_connect_swapped( object_number_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_number_spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_number_spin );

	label = gtk_label_new( _("Brightness mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_brightness_spin = gtk_spin_button_new_with_range( 0, 1000, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), parameters.object_brightness_mean );
	g_signal_connect( object_brightness_spin, "value_changed", G_CALLBACK( spin_object_brightness_mean ), NULL );
	g_signal_connect_swapped( object_brightness_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_brightness_spin, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_brightness_spin );

	label = gtk_label_new( _("Brightness sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_brightness_sigma_spin = gtk_spin_button_new_with_range( 0, 500, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), parameters.object_brightness_sigma );
	g_signal_connect( object_brightness_sigma_spin, "value_changed", G_CALLBACK( spin_object_brightness_sigma ), NULL );
	g_signal_connect_swapped( object_brightness_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_brightness_sigma_spin, 3, 4, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_brightness_sigma_spin );

	label = gtk_label_new( _("Color mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_color_spin = gtk_spin_button_new_with_range( 1000, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), parameters.object_color_mean );
	g_signal_connect( object_color_spin, "value_changed", G_CALLBACK( spin_object_color_mean ), NULL );
	g_signal_connect_swapped( object_color_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_color_spin, 1, 2, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_color_spin );

	label = gtk_label_new( _("Color sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_color_sigma_spin = gtk_spin_button_new_with_range( 1, 10000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), parameters.object_color_sigma );
	g_signal_connect( object_color_sigma_spin, "value_changed", G_CALLBACK( spin_object_color_sigma ), NULL );
	g_signal_connect_swapped( object_color_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_color_sigma_spin, 3, 4, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_color_sigma_spin );

	if ( parameters.object_x == 0 ) parameters.object_x = gimp_image_width( image_id )/2;
	label = gtk_label_new( _("Object center: X:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_center_x_spin = gtk_spin_button_new_with_range( 0, gimp_image_width( image_id ), 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), parameters.object_x );
	g_signal_connect( object_center_x_spin, "value_changed", G_CALLBACK( spin_object_x ), NULL );
	g_signal_connect_swapped( object_center_x_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_center_x_spin, 1, 2, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_center_x_spin );

	if ( parameters.object_y == 0 ) parameters.object_y = gimp_image_height( image_id )/2;
	label = gtk_label_new( _("Y:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	object_center_y_spin = gtk_spin_button_new_with_range( 0, gimp_image_height( image_id ), 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), parameters.object_y );
	g_signal_connect( object_center_y_spin, "value_changed", G_CALLBACK( spin_object_y ), NULL );
	g_signal_connect_swapped( object_center_y_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), object_center_y_spin, 3, 4, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( object_center_y_spin );

	table = gtk_table_new( 2, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_box_pack_start( GTK_BOX( object_box ), table, FALSE, FALSE, 0 );
	gtk_widget_show( table );

	label = gtk_label_new( _("Density profile:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_label_set_justify( GTK_LABEL( label ), GTK_JUSTIFY_LEFT );
	gtk_widget_show( label );

	object_density_combo = gimp_int_combo_box_new(
		_("Random"), DENSITY_RANDOM,
		_("Gauss"), DENSITY_GAUSS,
		_("Plummer"), DENSITY_PLUMMER,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), parameters.object_density );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( object_density_combo ), parameters.object_density,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.object_density );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( object_density_combo ), parameters.object_density,
		G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), object_density_combo, 1, 2, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( object_density_combo );

	object_radius_adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Radius parameter:"), 125, 75,
		parameters.object_radius, 0, 200, 1, 5, 0,
		TRUE, 0, 0, _("Radius of object in % of the image (unused for random distribution)"), NULL );
	g_signal_connect( object_radius_adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.object_radius );
	g_signal_connect_swapped( object_radius_adj, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );

	gtk_widget_show( object_box );


	frame = gimp_frame_new( _("Foreground Stars") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 2, 6, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("Amount:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	foreground_number_spin = gtk_spin_button_new_with_range( 0, 10000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), parameters.foreground_stars );
	g_signal_connect( foreground_number_spin, "value_changed", G_CALLBACK( spin_foreground_stars ), NULL );
	g_signal_connect_swapped( foreground_number_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), foreground_number_spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( foreground_number_spin );

	label = gtk_label_new( _("Brightness mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	foreground_brightness_spin = gtk_spin_button_new_with_range( 0, 1000, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), parameters.foreground_brightness_mean );
	g_signal_connect( foreground_brightness_spin, "value_changed", G_CALLBACK( spin_foreground_brightness_mean ), NULL );
	g_signal_connect_swapped( foreground_brightness_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), foreground_brightness_spin, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( foreground_brightness_spin );

	label = gtk_label_new( _("Brightness sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	foreground_brightness_sigma_spin = gtk_spin_button_new_with_range( 0, 500, 0.1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), parameters.foreground_brightness_sigma );
	g_signal_connect( foreground_brightness_sigma_spin, "value_changed", G_CALLBACK( spin_foreground_brightness_sigma ), NULL );
	g_signal_connect_swapped( foreground_brightness_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), foreground_brightness_sigma_spin, 3, 4, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( foreground_brightness_sigma_spin );

	label = gtk_label_new( _("Color mean:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	foreground_color_spin = gtk_spin_button_new_with_range( 1000, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), parameters.foreground_color_mean );
	g_signal_connect( foreground_color_spin, "value_changed", G_CALLBACK( spin_foreground_color_mean ), NULL );
	g_signal_connect_swapped( foreground_color_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), foreground_color_spin, 1, 2, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( foreground_color_spin );

	label = gtk_label_new( _("Color sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	foreground_color_sigma_spin = gtk_spin_button_new_with_range( 1, 10000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), parameters.foreground_color_sigma );
	g_signal_connect( foreground_color_sigma_spin, "value_changed", G_CALLBACK( spin_foreground_color_sigma ), NULL );
	g_signal_connect_swapped( foreground_color_sigma_spin, "value_changed", G_CALLBACK( recalculation_necessary ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 2, 3, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), foreground_color_sigma_spin, 3, 4, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( foreground_color_sigma_spin );

/* Rendering options */

	right_vbox = gtk_box_new (GTK_ORIENTATION_VERTICAL, 12);
	gtk_container_set_border_width (GTK_CONTAINER (right_vbox), 12);
	gtk_notebook_append_page (GTK_NOTEBOOK (notebook), right_vbox,
                            gtk_label_new_with_mnemonic (_("R_endering")));
	gtk_widget_show (right_vbox);

	frame = gimp_frame_new( _("Options") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 12, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("PSF:") );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_widget_show( label );

	GtkWidget *psf = gimp_int_combo_box_new(
		_("Delta peak"), PSF_DELTA_PEAK,
		_("Gauss"), PSF_GAUSS,
		_("Gauss with 4 diffraction lines"), PSF_GAUSS_4_DIFFRACTIONS,
		_("Gauss with 6 diffraction lines"), PSF_GAUSS_6_DIFFRACTIONS,
		_("Gauss with 12 diffraction lines"), PSF_GAUSS_12_DIFFRACTIONS,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( psf ), parameters.psf );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( psf ), parameters.psf,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.psf );
	g_signal_connect_swapped( psf, "changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	gtk_table_attach( GTK_TABLE( table ), psf, 1, 3, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( psf );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 1,
		_("Sigma:"), 185, 75,
		parameters.sigma_psf, 0.1, 10., 0.1, 5, 3,
		TRUE, 0, 0, _("Sigma of gauss function in px"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.sigma_psf );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Stars with diffraction lines:"), 185, 75,
		parameters.diffraction_percentage, 0.0, 100.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Percentage of stars that have diffraction lines"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.diffraction_percentage );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( find_star_norm ), NULL );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 3,
		_("Diffraction angle:"), 185, 75,
		parameters.diffraction_angle, 0.0, 90., 0.1, 5, 1,
		TRUE, 0, 0, _("Angle of diffraction lines in degrees"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.diffraction_angle );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 4,
		_("Diffraction length:"), 185, 75,
		parameters.diffraction_length, 0.0, 10.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Length of diffraction lines"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.diffraction_length );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 5,
		_("Diffraction color:"), 185, 75,
		parameters.diffraction_color, 0.0, 2.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Diffraction lines usually show a color variation"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.diffraction_color );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 6,
		_("Noise:"), 185, 75,
		parameters.noise, 0, 200, 1, 5, 0,
		TRUE, 0, 0, _("Noise in % of sqrt(N) = photon noise"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.noise );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 7,
		_("Background:"), 185, 75,
		parameters.background, 0, 100, 1, 5, 0,
		TRUE, 0, 0, _("Background pixel value offset"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.background );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 8,
		_("Burn out:"), 185, 75,
		parameters.burnout, 0.00, 100.00, 0.01, 5, 2,
		TRUE, 0, 0, _("% of stars that burn out (0=normalize to brightest star, 100=to darkest one)"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.burnout );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( find_star_norm ), NULL );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 9,
		_("Shininess:"), 185, 75,
		parameters.shininess, 0.0, 2.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Burnt out stars have larger sigma"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.shininess );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	check_box = gtk_check_button_new_with_label( _("Split foreground, object and background to layers") );
	gtk_table_attach_defaults( GTK_TABLE( table ), check_box, 0, 2, 10, 11 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( check_box ), parameters.split_layers );
	g_signal_connect( check_box, "toggled", G_CALLBACK( split_layers_toggled ), NULL );
	g_signal_connect_swapped( check_box, "toggled", G_CALLBACK( gimp_preview_invalidate ), preview );
	gtk_widget_show( check_box );

	check_box = gtk_check_button_new_with_label( _("Create an illumination mask") );
	gtk_table_attach_defaults( GTK_TABLE( table ), check_box, 0, 2, 11, 12 );
	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( check_box ), parameters.object_mask );
	g_signal_connect( check_box, "toggled", G_CALLBACK( object_mask_toggled ), NULL );
	g_signal_connect_swapped( check_box, "toggled", G_CALLBACK( gimp_preview_invalidate ), preview );
/*	gtk_widget_show( check_box ); */


	gtk_widget_show( left_vbox );
	gtk_widget_show( right_vbox );
	gtk_widget_show( main_hbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
