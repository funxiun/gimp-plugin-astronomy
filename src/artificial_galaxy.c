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
gimp artificial galaxy plug-in
(C) Georg Hennig <georg.hennig@web.de>
Creates an artificial galaxy (elliptical - boxy or disky, simple spiral or
spiral barred with bulge.
*/

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include <libgimp/gimp.h>
#include <libgimp/gimpui.h>
#include <gtk/gtk.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "artificial_stars_temperature.h"

#include "plugin-intl.h"

#define PLUG_IN_NAME "gimp-plugin-astro-artificial-galaxy"
#define PLUG_IN_VERSION "0.10"
#define PLUG_IN_DATE "09.2018"

double CLAMP_DOUBLE( const double a, const double b, const double c )
{
	return ( a < b ) ? b : ( ( a > c ) ? c : a );
}

enum GALAXY_TYPE
{
	GALAXY_SPIRAL = 0,
	GALAXY_SPIRALBARRED,
	GALAXY_ELLIPTICAL
};

enum SPIRAL_TYPE
{
	SPIRAL_LOGARITHMIC,
	SPIRAL_ARCHIMEDIC,
	SPIRAL_CIRCLES
};

/*
parameters
*/
typedef struct
{
	gdouble phi;
	gdouble psi;
	gdouble theta;

	gdouble transparency_bright;
	gdouble transparency_dark;

	gint32 noise;
	gdouble multiplier;

	gboolean split_layers;

	guint32 random_seed;
	gboolean random_seed_bool;

	guint32 galaxy_type;

	gint32 object_x;
	gint32 object_y;


	gint32 spiral_bright1_objects;
	gint32 spiral_bright2_objects;
	gint32 spiral_bright3_objects;
	gint32 spiral_dark_objects;

	gint32 spiral_spiral_type;

	gdouble spiral_bulge_radius;
	gdouble spiral_spiral_radius;
	gdouble spiral_spiral_thickness;
	gdouble spiral_spiral_turns;

	gdouble spiral_spiral_randomness;
	gdouble spiral_angle_randomness;
	gdouble spiral_radius_randomness;

	gdouble spiral_spiral_lumpiness;
	gdouble spiral_disc_lumpiness;

	gint32 spiral_bulge_color_mean;

	GimpRGB spiral_color_bright1;
	gint32 spiral_color_bright1_sigma;
	GimpRGB spiral_color_bright2;
	gint32 spiral_color_bright2_sigma;
	GimpRGB spiral_color_bright3;
	gint32 spiral_color_bright3_sigma;
	GimpRGB spiral_color_dark;
	gint32 spiral_color_dark_sigma;


	gint32 spiralbarred_bright1_objects;
	gint32 spiralbarred_bright2_objects;
	gint32 spiralbarred_bright3_objects;
	gint32 spiralbarred_dark_objects;

	gint32 spiralbarred_spiral_type;

	gdouble spiralbarred_bulge_radius;
	gdouble spiralbarred_bar_length;
	gdouble spiralbarred_spiral_radius;
	gdouble spiralbarred_spiral_thickness;
	gdouble spiralbarred_spiral_turns;

	gdouble spiralbarred_spiral_randomness;
	gdouble spiralbarred_angle_randomness;
	gdouble spiralbarred_radius_randomness;

	gdouble spiralbarred_spiral_lumpiness;
	gdouble spiralbarred_disc_lumpiness;

	gint32 spiralbarred_bulge_color_mean;

	GimpRGB spiralbarred_color_bright1;
	gint32 spiralbarred_color_bright1_sigma;
	GimpRGB spiralbarred_color_bright2;
	gint32 spiralbarred_color_bright2_sigma;
	GimpRGB spiralbarred_color_bright3;
	gint32 spiralbarred_color_bright3_sigma;
	GimpRGB spiralbarred_color_dark;
	gint32 spiralbarred_color_dark_sigma;


	gdouble elliptical_radius;
	gdouble elliptical_excentricity;
	gdouble elliptical_boxiness;

	gint32 elliptical_color_mean;


	gboolean show_preview;
} tparameter;


static tparameter parameters =
{
/* phi, psi, theta */
	0.,
	0.,
	0.,

/* transparency bright, dark */
	50.,
	0.,

/* noise, multiplier */
	20,
	1.,

/* split layers */
	FALSE,

/* random seed */
	0,
	FALSE,

/* galaxy type */
	GALAXY_ELLIPTICAL,

/* x, y */
	0,
	0,


/* spiral bright1, bright2, bright3, dark objects */
	1000,
	1000,
	1000,
	1000,

/* spiral spiral type */
	SPIRAL_LOGARITHMIC,

/* spiral bulge, spiral radius, spiral thickness, spiral turns */
	10.,
	50.,
	5.,
	1.,

/* spiral spiral, angle, radius randomness */
	0.2,
	0.4,
	0.1,

/* spiral spiral, disc lumpiness */
	0.4,
	0.7,

/* spiral bulge color mean */
	6000,

/* spiral color bright1, bright2, bright3, dark */
	{ 0.67, 0.94, 1.00, 1. },
	2,
	{ 0.71, 0.61, 0.50, 1. },
	4,
	{ 1.00, 0.40, 0.46, 1. },
	2,
	{ 0.25, 0.23, 0.22, 1. },
	5,


/* spiral barred bright1, bright2, bright3, dark objects */
	1000,
	1000,
	1000,
	1000,

/* spiral barred spiral type */
	SPIRAL_LOGARITHMIC,

/* spiral barred bulge, bar, spiral radius, spiral thickness, spiral turns */
	10.,
	15.,
	50.,
	5.,
	1.,

/* spiral barred spiral, angle, radius randomness */
	0.2,
	0.4,
	0.1,

/* spiral barred spiral, disc lumpiness */
	0.4,
	0.7,

/* spiral barred bulge color mean */
	6000,

/* spiral barred color bright1, bright2, bright3, dark with sigma */
	{ 0.67, 0.94, 1.00, 1. },
	2,
	{ 0.71, 0.61, 0.50, 1. },
	4,
	{ 1.00, 0.40, 0.46, 1. },
	2,
	{ 0.25, 0.23, 0.22, 1. },
	5,


/* elliptical radius, excentricity, boxiness */
	20,
	0.8,
	0.02,

/* elliptical color mean */
	4600,


/* preview */
	TRUE
};

/*
Prototypes
*/
gint32 objects_bright1_number = 0;
gdouble *objects_bright1 = NULL;
gint32 objects_bright2_number = 0;
gdouble *objects_bright2 = NULL;
gint32 objects_bright3_number = 0;
gdouble *objects_bright3 = NULL;
gint32 objects_dark_number = 0;
gdouble *objects_dark = NULL;
gint32 objects_all_number = 0;
gdouble *objects_all = NULL;

gdouble brightness_object_norm = 0.;

void object_new( gdouble **list, gint32 *list_number, const gint32 list_number_new );
void object_free( gdouble **list, gint32 *list_number );
void object_add( gdouble **list, const gint32 list_number, const gint32 list_position,
	const gdouble x, const gdouble y, const gdouble z, const gdouble brightness,
	const gdouble r, const gdouble g, const gdouble b,
	const gdouble sigma_x, const gdouble sigma_y, const gdouble rotation );
void object_get( gdouble **list, const gint32 list_position,
	gdouble *x, gdouble *y, gdouble *z, gdouble *brightness,
	gdouble *r, gdouble *g, gdouble *b,
	gdouble *sigma_x, gdouble *sigma_y, gdouble *rotation );

/*
communication to the GIMP
*/
static void query( void );

static void run( const gchar *name, gint nparams, const GimpParam *param,
	gint *nreturn_vals, GimpParam **return_vals );

/*
create galaxy
*/
void create_object_distribution();

void draw_bulge_or_ellipse( const gint galaxy_type, gdouble *image_buffer, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, const gboolean show_progress );

void draw_objects( gdouble *image_buffer, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, gdouble **list, const gint32 list_number, const gboolean show_progress );

void normalize( gdouble *image_buffer, gint32 x_size, gint32 y_size, const gboolean show_progress );

void create_noise( gdouble *image_buffer, gint32 x_size, gint32 y_size, const gboolean show_progress );

static void create_galaxy();

static void create_galaxy_selection( GimpPixelRgn *region_destination, const gboolean show_progress );

/*
user interface
*/
GtkWidget *preview;
GtkWidget *progress_bar;

GtkWidget *recalculate_button;
gboolean recalculating = FALSE;
gboolean cancel_recalculation = FALSE;

GtkWidget *galaxy_type;

GtkWidget *galaxy_specific_frame;
GtkWidget *galaxy_current_table;

GtkWidget *galaxy_spiral_table;
GtkWidget *galaxy_spiralbarred_table;
GtkWidget *galaxy_elliptical_table;

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
		{ GIMP_PDB_FLOAT, "phi", "Phi angle" },
		{ GIMP_PDB_FLOAT, "psi", "Psi angle" },
		{ GIMP_PDB_FLOAT, "theta", "Theta angle" },
		{ GIMP_PDB_FLOAT, "transparency_bright", "Transparency of bright clouds" },
		{ GIMP_PDB_FLOAT, "transparency_dark", "Transparency of dark clouds" },
		{ GIMP_PDB_INT32, "noise", "Noise in % of sqrt(N)" },
		{ GIMP_PDB_FLOAT, "multiplier", "Multiplier to all pixel values" },
		{ GIMP_PDB_INT32, "split_layers", "Split objects to layers" },
		{ GIMP_PDB_INT32, "random_seed", "Seed number" },
		{ GIMP_PDB_INT32, "random_seed_bool", "Use a random seed number" },
		{ GIMP_PDB_INT32, "galaxy_type", "Galaxy type (0=spiral,1=spiral barred,2=elliptical)" },
		{ GIMP_PDB_INT32, "object_x", "Center of object (x)" },
		{ GIMP_PDB_INT32, "object_y", "Center of object (y)" },
		{ GIMP_PDB_INT32, "spiral_bright1_objects", "Number of spiral galaxy's bright objects (1)" },
		{ GIMP_PDB_INT32, "spiral_bright2_objects", "Number of spiral galaxy's bright objects (2)" },
		{ GIMP_PDB_INT32, "spiral_bright3_objects", "Number of spiral galaxy's bright objects (3)" },
		{ GIMP_PDB_INT32, "spiral_dark_objects", "Number of spiral galaxy's dark objects" },
		{ GIMP_PDB_INT32, "spiral_spiral_type", "Spiral galaxy's spiral type (0=log.,1=archim.,2=circles)" },
		{ GIMP_PDB_FLOAT, "spiral_bulge_radius", "Radius of spiral galaxy's bulge" },
		{ GIMP_PDB_FLOAT, "spiral_spiral_radius", "Radius of spiral galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiral_spiral_thickness", "Thickness of spiral galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiral_spiral_turns", "Turns of spiral galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiral_spiral_randomness", "Randomness of spiral galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiral_angle_randomness", "Randomness of spiral galaxy's angle" },
		{ GIMP_PDB_FLOAT, "spiral_radius_randomness", "Randomness of spiral galaxy's radius" },
		{ GIMP_PDB_FLOAT, "spiral_spiral_lumpiness", "Lumpiness of spiral galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiral_disc_lumpiness", "Lumpiness of spiral galaxy's disc" },
		{ GIMP_PDB_INT32, "spiral_bulge_color_mean", "Color mean value of spiral galaxy's bulge in K" },
		{ GIMP_PDB_COLOR, "spiral_color_bright1", "Color of spiral galaxy's bright objects (1)" },
		{ GIMP_PDB_INT32, "spiral_color_bright1_sigma", "Sigma of color of spiral galaxy's bright objects (1)" },
		{ GIMP_PDB_COLOR, "spiral_color_bright2", "Color of spiral galaxy's bright objects (2)" },
		{ GIMP_PDB_INT32, "spiral_color_bright2_sigma", "Sigma of color of spiral galaxy's bright objects (2)" },
		{ GIMP_PDB_COLOR, "spiral_color_bright3", "Color of spiral galaxy's bright objects (3)" },
		{ GIMP_PDB_INT32, "spiral_color_bright3_sigma", "Sigma of color of spiral galaxy's bright objects (3)" },
		{ GIMP_PDB_COLOR, "spiral_color_dark", "Color of spiral galaxy's dark objects" },
		{ GIMP_PDB_INT32, "spiral_color_dark_sigma", "Sigma of color of spiral galaxy's dark objects" },
		{ GIMP_PDB_INT32, "spiralbarred_bright1_objects", "Number of spiral barred galaxy's bright objects (1)" },
		{ GIMP_PDB_INT32, "spiralbarred_bright2_objects", "Number of spiral barred galaxy's bright objects (2)" },
		{ GIMP_PDB_INT32, "spiralbarred_bright3_objects", "Number of spiral barred galaxy's bright objects (3)" },
		{ GIMP_PDB_INT32, "spiralbarred_dark_objects", "Number of spiral barred galaxy's dark objects" },
		{ GIMP_PDB_INT32, "spiralbarred_spiral_type", "Spiral barred galaxy's spiral type (0=log.,1=archim.,2=circles)" },
		{ GIMP_PDB_FLOAT, "spiralbarred_bulge_radius", "Radius of spiral barred galaxy's bulge" },
		{ GIMP_PDB_FLOAT, "spiralbarred_bar_length", "Radius of spiral barred galaxy's bar" },
		{ GIMP_PDB_FLOAT, "spiralbarred_spiral_radius", "Radius of spiral barred galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiralbarred_spiral_thickness", "Thickness of spiral barred galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiralbarred_spiral_turns", "Turns of spiral barred galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiralbarred_spiral_randomness", "Randomness of spiral barred galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiralbarred_angle_randomness", "Randomness of spiral barred galaxy's angle" },
		{ GIMP_PDB_FLOAT, "spiralbarred_radius_randomness", "Randomness of spiral barred galaxy's radius" },
		{ GIMP_PDB_FLOAT, "spiralbarred_spiral_lumpiness", "Lumpiness of spiral barred galaxy's spiral" },
		{ GIMP_PDB_FLOAT, "spiralbarred_disc_lumpiness", "Lumpiness of spiral barred galaxy's disc" },
		{ GIMP_PDB_INT32, "spiralbarred_bulge_color_mean", "Color mean value of spiral barred galaxy's bulge in K" },
		{ GIMP_PDB_COLOR, "spiralbarred_color_bright1", "Color of spiral barred galaxy's bright objects (1)" },
		{ GIMP_PDB_INT32, "spiralbarred_color_bright1_sigma", "Sigma of color of spiral barred galaxy's bright objects (1)" },
		{ GIMP_PDB_COLOR, "spiralbarred_color_bright2", "Color of spiral barred galaxy's bright objects (2)" },
		{ GIMP_PDB_INT32, "spiralbarred_color_bright2_sigma", "Sigma of color of spiral barred galaxy's bright objects (2)" },
		{ GIMP_PDB_COLOR, "spiralbarred_color_bright3", "Color of spiral barred galaxy's bright objects (3)" },
		{ GIMP_PDB_INT32, "spiralbarred_color_bright3_sigma", "Sigma of color of spiral barred galaxy's bright objects (3)" },
		{ GIMP_PDB_COLOR, "spiralbarred_color_dark", "Color of spiral barred galaxy's dark objects" },
		{ GIMP_PDB_INT32, "spiralbarred_color_dark_sigma", "Sigma of color of spiral barred galaxy's dark objects" },
		{ GIMP_PDB_FLOAT, "elliptical_radius", "Radius of elliptical galaxy" },
		{ GIMP_PDB_FLOAT, "elliptical_excentricity", "Excentricity of elliptical galaxy" },
		{ GIMP_PDB_FLOAT, "elliptical_boxiness", "Boxiness of elliptical galaxy" },
		{ GIMP_PDB_INT32, "elliptical_color_mean", "Color mean value of elliptical galaxy in K" },
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
		_("Create an artificial galaxy"),
		_("This plug-in creates an artificial galaxy."),
		"Georg Hennig <georg.hennig@web.de>",
		"Georg Hennig <georg.hennig@web.de>",
		PLUG_IN_DATE,
		_("Artificial Galaxy"),
		"RGB*",
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

	gboolean added_alpha_channel = FALSE;
	gint32 active_layer = gimp_image_get_active_layer( image_id );

	switch( run_mode )
	{
		case GIMP_RUN_INTERACTIVE:
			gimp_get_data( PLUG_IN_NAME, &parameters );

			if ( gimp_drawable_type( param[2].data.d_drawable ) == GIMP_RGB_IMAGE )
			{
				gimp_image_undo_disable( image_id );
				gimp_layer_add_alpha( active_layer );
				added_alpha_channel = TRUE;
			}

			if ( !dialog( param[1].data.d_image, gimp_drawable_get( param[2].data.d_drawable ) ) )
			{
				if ( added_alpha_channel )
				{
					gimp_layer_flatten( active_layer );
					gimp_image_undo_enable( image_id );
				}
				return;
			}

			if ( added_alpha_channel )
			{
				gimp_layer_flatten( active_layer );
				gimp_image_undo_enable( image_id );
			}

			gimp_set_data( PLUG_IN_NAME, &parameters, sizeof( tparameter ) );
			break;
		case GIMP_RUN_NONINTERACTIVE:
			if ( nparams != 67 )
			{
				status = GIMP_PDB_CALLING_ERROR;
			}
			else
			{
				image_id = param[1].data.d_image;

				parameters.phi = param[2].data.d_float;
				parameters.psi = param[3].data.d_float;
				parameters.theta = param[4].data.d_float;

				parameters.transparency_bright = param[5].data.d_float;
				parameters.transparency_dark = param[6].data.d_float;

				parameters.noise = param[7].data.d_int32;
				parameters.multiplier = param[8].data.d_float;

				parameters.split_layers = param[9].data.d_int32;

				parameters.random_seed = param[10].data.d_int32;
				parameters.random_seed_bool = param[11].data.d_int32;

				parameters.galaxy_type = param[12].data.d_int32;

				parameters.object_x = param[13].data.d_int32;
				parameters.object_y = param[14].data.d_int32;


				parameters.spiral_bright1_objects = param[15].data.d_int32;
				parameters.spiral_bright2_objects = param[16].data.d_int32;
				parameters.spiral_bright3_objects = param[17].data.d_int32;
				parameters.spiral_dark_objects = param[18].data.d_int32;

				parameters.spiral_spiral_type = param[19].data.d_int32;

				parameters.spiral_bulge_radius = param[20].data.d_float;
				parameters.spiral_spiral_radius = param[21].data.d_float;
				parameters.spiral_spiral_thickness = param[22].data.d_float;
				parameters.spiral_spiral_turns = param[23].data.d_float;

				parameters.spiral_spiral_randomness = param[24].data.d_float;
				parameters.spiral_angle_randomness = param[25].data.d_float;
				parameters.spiral_radius_randomness = param[26].data.d_float;

				parameters.spiral_spiral_lumpiness = param[27].data.d_float;
				parameters.spiral_disc_lumpiness = param[28].data.d_float;

				parameters.spiral_bulge_color_mean = param[29].data.d_int32;

				parameters.spiral_color_bright1 = param[30].data.d_color;
				parameters.spiral_color_bright1_sigma = param[31].data.d_int32;
				parameters.spiral_color_bright2 = param[32].data.d_color;
				parameters.spiral_color_bright2_sigma = param[33].data.d_int32;
				parameters.spiral_color_bright3 = param[34].data.d_color;
				parameters.spiral_color_bright3_sigma = param[35].data.d_int32;
				parameters.spiral_color_dark = param[36].data.d_color;
				parameters.spiral_color_dark_sigma = param[37].data.d_int32;


				parameters.spiralbarred_bright1_objects = param[38].data.d_int32;
				parameters.spiralbarred_bright2_objects = param[39].data.d_int32;
				parameters.spiralbarred_bright3_objects = param[40].data.d_int32;
				parameters.spiralbarred_dark_objects = param[41].data.d_int32;

				parameters.spiralbarred_spiral_type = param[42].data.d_int32;

				parameters.spiralbarred_bulge_radius = param[43].data.d_float;
				parameters.spiralbarred_bar_length = param[44].data.d_float;
				parameters.spiralbarred_spiral_radius = param[45].data.d_float;
				parameters.spiralbarred_spiral_thickness = param[46].data.d_float;
				parameters.spiralbarred_spiral_turns = param[47].data.d_float;

				parameters.spiralbarred_spiral_randomness = param[48].data.d_float;
				parameters.spiralbarred_angle_randomness = param[49].data.d_float;
				parameters.spiralbarred_radius_randomness = param[50].data.d_float;

				parameters.spiralbarred_spiral_lumpiness = param[51].data.d_float;
				parameters.spiralbarred_disc_lumpiness = param[52].data.d_float;

				parameters.spiralbarred_bulge_color_mean = param[53].data.d_int32;

				parameters.spiralbarred_color_bright1 = param[54].data.d_color;
				parameters.spiralbarred_color_bright1_sigma = param[55].data.d_int32;
				parameters.spiralbarred_color_bright2 = param[56].data.d_color;
				parameters.spiralbarred_color_bright2_sigma = param[57].data.d_int32;
				parameters.spiralbarred_color_bright3 = param[58].data.d_color;
				parameters.spiralbarred_color_bright3_sigma = param[59].data.d_int32;
				parameters.spiralbarred_color_dark = param[60].data.d_color;
				parameters.spiralbarred_color_dark_sigma = param[61].data.d_int32;


				parameters.elliptical_radius = param[62].data.d_float;
				parameters.elliptical_excentricity = param[63].data.d_float;
				parameters.elliptical_boxiness = param[64].data.d_float;

				parameters.elliptical_color_mean = param[65].data.d_int32;


				parameters.show_preview = param[66].data.d_int32;
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
		create_galaxy();
	}

	values[0].data.d_status = status;
}

void object_new( gdouble **list, gint32 *list_number, const gint32 list_number_new )
{
	*list = malloc( list_number_new * sizeof(gdouble) * 10 );

	*list_number = list_number_new;
}

void object_free( gdouble **list, gint32 *list_number )
{
	if ( !*list ) return;

	if ( *list_number == 0 ) return;

	free( *list );
	*list = NULL;
	*list_number = 0;
}

void object_add( gdouble **list, const gint32 list_number, const gint32 list_position,
	const gdouble x, const gdouble y, const gdouble z, const gdouble brightness,
	const gdouble r, const gdouble g, const gdouble b,
	const gdouble sigma_x, const gdouble sigma_y, const gdouble rotation )
{
	if ( list_position >= list_number ) return;

	(*list)[list_position*10+0] = x;
	(*list)[list_position*10+1] = y;
	(*list)[list_position*10+2] = z;
	(*list)[list_position*10+3] = brightness;
	(*list)[list_position*10+4] = r;
	(*list)[list_position*10+5] = g;
	(*list)[list_position*10+6] = b;
	(*list)[list_position*10+7] = sigma_x;
	(*list)[list_position*10+8] = sigma_y;
	(*list)[list_position*10+9] = rotation;
}

void object_get( gdouble **list, const gint32 list_position,
	gdouble *x, gdouble *y, gdouble *z, gdouble *brightness,
	gdouble *r, gdouble *g, gdouble *b,
	gdouble *sigma_x, gdouble *sigma_y, gdouble *rotation )
{
	if ( x ) *x = (*list)[list_position*10+0];
	if ( y ) *y = (*list)[list_position*10+1];
	if ( z ) *z = (*list)[list_position*10+2];
	if ( brightness ) *brightness = (*list)[list_position*10+3];
	if ( r ) *r = (*list)[list_position*10+4];
	if ( g ) *g = (*list)[list_position*10+5];
	if ( b ) *b = (*list)[list_position*10+6];
	if ( sigma_x ) *sigma_x = (*list)[list_position*10+7];
	if ( sigma_y ) *sigma_y = (*list)[list_position*10+8];
	if ( rotation ) *rotation = (*list)[list_position*10+9];
}

inline gdouble gauss( const gdouble x, const gdouble y, const gdouble A, const gdouble x0, const gdouble y0,
	const gdouble sigma_x, const gdouble sigma_y )
{
	/* z = A * exp(-0.5*(x-x0)^2/sigma_x^2) * exp(-0.5*(y-y0)^2/sigma_y^2) + b */
	return A * exp( -0.5*(x-x0)*(x-x0)/(sigma_x*sigma_x) ) * exp( -0.5*(y-y0)*(y-y0)/(sigma_y*sigma_y) );
}

inline gdouble gauss_rotation( const gdouble x, const gdouble y, const gdouble A, const gdouble x0, const gdouble y0,
	const gdouble sigma_x, const gdouble sigma_y, const gdouble rotation )
{
	gdouble x_new =  cos(rotation)*(x-x0) + sin(rotation)*(y-y0) + x0;
	gdouble y_new = -sin(rotation)*(x-x0) + cos(rotation)*(y-y0) + y0;
	/* z = A * exp(-0.5*(x-x0)^2/sigma_x^2) * exp(-0.5*(y-y0)^2/sigma_y^2) + b */
	return A * exp( -0.5*(x_new-x0)*(x_new-x0)/(sigma_x*sigma_x) ) * exp( -0.5*(y_new-y0)*(y_new-y0)/(sigma_y*sigma_y) );
}

int comp( const void *ptr1, const void *ptr2 ) /* sort by brightness */
{
	const gdouble *p_1 = (gdouble*)ptr1;
	const gdouble *p_2 = (gdouble*)ptr2;

	if ( ( p_2[3] - p_1[3] ) > 0. ) return 1;
	if ( ( p_2[3] - p_1[3] ) < 0. ) return -1;
	return 0;
}

gdouble distance_spiral( const gdouble x, const gdouble y, const gint32 galaxy_type, const gdouble *spiral_radius_parameter,
	const gdouble *spiral_angle_parameter, const gdouble *spiral_rotation, const gdouble *spiral_min, const gdouble *spiral_max )
{
	gdouble dist = 100000.;
	gdouble dist_tmp;
	gdouble t, spiral;
	gdouble x_param, y_param;
	gint32 i, j;
	for ( i=0; i<1/*4*/; i++ )
	{
		for ( j=1; j<10000; j++ )
		{
			switch ( galaxy_type )
			{
				case SPIRAL_LOGARITHMIC:
					t = j * (spiral_max[i]-spiral_min[i])/10000;
					spiral = spiral_radius_parameter[i] * exp( spiral_angle_parameter[i] * t );
					x_param =  spiral*cos(t)*cos(spiral_rotation[i]) + spiral*sin(t)*sin(spiral_rotation[i]);
					y_param = -spiral*cos(t)*sin(spiral_rotation[i]) * spiral*sin(t)*cos(spiral_rotation[i]);
					dist_tmp = sqrt( (x-x_param)*(x-x_param) + (y-y_param)*(y-y_param) );
					if ( dist > dist_tmp ) dist = dist_tmp;
					break;
				case SPIRAL_ARCHIMEDIC:
					break;
				case SPIRAL_CIRCLES:
					break;
				default:
					break;
			}
		}
	}
printf("distance: %f\n",dist);
	return dist;
}

void create_object_distribution()
{
	object_free( &objects_bright1, &objects_bright1_number );
	object_free( &objects_bright2, &objects_bright2_number );
	object_free( &objects_bright3, &objects_bright3_number );
	object_free( &objects_dark, &objects_dark_number );

	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed+10 );

/* First two spiral arms are "main" spirals, the other ones random ones, number depends on spiral randomness */
	gdouble *spiral_radius_parameter;
	gdouble *spiral_angle_parameter;
	gdouble *spiral_rotation;
	gdouble *spiral_min;
	gdouble *spiral_max;

	gdouble x, y, z, brightness, color_r, color_g, color_b, sigma_x, sigma_y, rotation;
	gdouble rmax, rmin;
	gdouble radius_random, angle_random, x_random, y_random;

	gdouble random, rotation_random, sigma_random;
	gint32 i, j, k;

	switch( parameters.galaxy_type )
	{
		case GALAXY_SPIRAL:
		{
			gint number_of_spirals = 2 + (gint)(parameters.spiral_spiral_randomness/0.1);

			spiral_radius_parameter = malloc( number_of_spirals*sizeof(gdouble) );
			spiral_angle_parameter = malloc( number_of_spirals*sizeof(gdouble) );
			spiral_rotation = malloc( number_of_spirals*sizeof(gdouble) );
			spiral_min = malloc( number_of_spirals*sizeof(gdouble) );
			spiral_max = malloc( number_of_spirals*sizeof(gdouble) );

			object_new( &objects_bright1, &objects_bright1_number, parameters.spiral_bright1_objects );
			object_new( &objects_bright2, &objects_bright2_number, parameters.spiral_bright2_objects );
			object_new( &objects_bright3, &objects_bright3_number, parameters.spiral_bright3_objects );
			object_new( &objects_dark, &objects_dark_number, parameters.spiral_dark_objects );

			spiral_min[0] = 0.1;
			spiral_max[0] = parameters.spiral_spiral_turns * 2 * M_PI;

			spiral_min[1] = 0.1;
			spiral_max[1] = parameters.spiral_spiral_turns * 2 * M_PI;

			for ( j=2; j<number_of_spirals; j++ )
			{
				spiral_min[j] = ( 1. / ( 1. + parameters.spiral_spiral_randomness ) ) * 0.5 * gsl_rng_uniform( r ) * parameters.spiral_spiral_turns * 2 * M_PI;
				spiral_max[j] = parameters.spiral_spiral_randomness * ( 0.5 + 0.5 * gsl_rng_uniform( r ) ) * parameters.spiral_spiral_turns * 2 * M_PI;
			}

/*		r = a * exp( k * phi )
			k = ln( r_max / r_min ) / ( phi_max - phi_min )
			a = r_max / exp( k * phi_max )
*/

			rmin = 0.2 * parameters.spiral_bulge_radius * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100;
			rmax = parameters.spiral_spiral_radius * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100;

			spiral_angle_parameter[0] = log( rmax / rmin ) / ( spiral_max[0] - spiral_min[0] );
			spiral_radius_parameter[0] = rmax / exp( spiral_angle_parameter[0] * spiral_max[0] );
			spiral_rotation[0] = 0.;

			spiral_angle_parameter[1] = log( rmax / rmin ) / ( spiral_max[1] - spiral_min[1] );
			spiral_radius_parameter[1] = rmax / exp( spiral_angle_parameter[1] * spiral_max[1] );
			spiral_rotation[1] = M_PI;

			for ( j=2; j<number_of_spirals; j++ )
			{
				spiral_angle_parameter[j] = log( rmax / rmin ) / ( spiral_max[j] - spiral_min[j] ) * gsl_rng_uniform( r );
				spiral_radius_parameter[j] = rmax / exp( spiral_angle_parameter[j] * spiral_max[j] );
				spiral_rotation[j] = 2 * M_PI * gsl_rng_uniform( r );
			}

			for ( i=0; i<parameters.spiral_bright1_objects; i++ )
			{
				random = gsl_rng_uniform( r );
				for ( j=0; j<number_of_spirals; j++ )
				{
					k = j;
					if ( random < (gdouble)(j+1)/number_of_spirals ) break;
				}

				angle_random = ( spiral_max[k] - spiral_min[k] ) * gsl_rng_uniform( r );
				radius_random = spiral_radius_parameter[k] * exp( spiral_angle_parameter[k] * angle_random );
				radius_random += fabs( 0.5 * parameters.spiral_radius_randomness * gsl_ran_gaussian( r, ( rmax - rmin ) / 2 ) );
				rotation_random = spiral_rotation[k] + parameters.spiral_angle_randomness * M_PI * gsl_rng_uniform( r );
				x_random = radius_random * ( cos( angle_random )*cos( rotation_random ) + sin( angle_random )*sin( rotation_random ) );
				y_random = radius_random * (-cos( angle_random )*sin( rotation_random ) + sin( angle_random )*cos( rotation_random ) );

				x = x_random + parameters.object_x;
				y = y_random + parameters.object_y;

				z = CLAMP_DOUBLE( MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 2 +
					0.1 * gsl_ran_gaussian( r, parameters.spiral_spiral_thickness * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ), 0., MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );

				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 0.06 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				sigma_x = sigma_random;
				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 2*sigma_x/*0.08 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100*/ ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				if ( sigma_random > sigma_x )
				{
					sigma_y = sigma_random;
				}
				else
				{
					sigma_y = sigma_x;
					sigma_x = sigma_random;
				}
				CLAMP_DOUBLE( sigma_y, sigma_x, 4*sigma_x );

				rotation = atan2( 1., spiral_angle_parameter[k] ) - M_PI/2 + angle_random + gsl_ran_gaussian( r, M_PI/16 );

				brightness = CLAMP_DOUBLE( 255. * ( 3. / pow( ( sigma_x + sigma_y ), 0.8 ) + 0.2 * gsl_rng_uniform( r ) ), 0., 255. );

				color_r = parameters.spiral_color_bright1.r;
				color_g = parameters.spiral_color_bright1.g;
				color_b = parameters.spiral_color_bright1.b;
				color_r += gsl_ran_gaussian( r, parameters.spiral_color_bright1_sigma )/255;
				color_g += gsl_ran_gaussian( r, parameters.spiral_color_bright1_sigma )/255;
				color_b += gsl_ran_gaussian( r, parameters.spiral_color_bright1_sigma )/255;

				object_add( &objects_bright1, objects_bright1_number, i, x, y, z, brightness,
					color_r, color_g, color_b, sigma_x, sigma_y, rotation );

				if ( i%(1+parameters.spiral_bright1_objects/100) ) reporter( "", (gdouble)(i+1)/parameters.spiral_bright1_objects );
			}

			for ( i=0; i<parameters.spiral_bright2_objects; i++ )
			{
				random = gsl_rng_uniform( r );
				for ( j=0; j<number_of_spirals; j++ )
				{
					k = j;
					if ( random < (gdouble)(j+1)/number_of_spirals ) break;
				}

				angle_random = ( spiral_max[k] - spiral_min[k] ) * gsl_rng_uniform( r );
				radius_random = spiral_radius_parameter[k] * exp( spiral_angle_parameter[k] * angle_random );
				radius_random += fabs( 0.5 * parameters.spiral_radius_randomness * gsl_ran_gaussian( r, ( rmax - rmin ) / 2 ) );
				rotation_random = spiral_rotation[k] + parameters.spiral_angle_randomness * M_PI * gsl_rng_uniform( r );
				x_random = radius_random * ( cos( angle_random )*cos( rotation_random ) + sin( angle_random )*sin( rotation_random ) );
				y_random = radius_random * (-cos( angle_random )*sin( rotation_random ) + sin( angle_random )*cos( rotation_random ) );

				x = x_random + parameters.object_x;
				y = y_random + parameters.object_y;

				z = CLAMP_DOUBLE( MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 2 +
					0.1 * gsl_ran_gaussian( r, parameters.spiral_spiral_thickness * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ), 0., MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );

				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 0.06 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				sigma_x = sigma_random;
				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 2*sigma_x/*0.08 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100*/ ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				if ( sigma_random > sigma_x )
				{
					sigma_y = sigma_random;
				}
				else
				{
					sigma_y = sigma_x;
					sigma_x = sigma_random;
				}
				CLAMP_DOUBLE( sigma_y, sigma_x, 4*sigma_x );

				rotation = atan2( 1., spiral_angle_parameter[k] ) - M_PI/2 + angle_random + gsl_ran_gaussian( r, M_PI/16 );

				brightness = CLAMP_DOUBLE( 255. * ( 3. / pow( ( sigma_x + sigma_y ), 0.8 ) + 0.2 * gsl_rng_uniform( r ) ), 0., 255. );

				color_r = parameters.spiral_color_bright2.r;
				color_g = parameters.spiral_color_bright2.g;
				color_b = parameters.spiral_color_bright2.b;
				color_r += gsl_ran_gaussian( r, parameters.spiral_color_bright2_sigma )/255;
				color_g += gsl_ran_gaussian( r, parameters.spiral_color_bright2_sigma )/255;
				color_b += gsl_ran_gaussian( r, parameters.spiral_color_bright2_sigma )/255;

				object_add( &objects_bright2, objects_bright2_number, i, x, y, z, brightness,
					color_r, color_g, color_b, sigma_x, sigma_y, rotation );

				if ( i%(1+parameters.spiral_bright2_objects/100) ) reporter( "", (gdouble)(i+1)/parameters.spiral_bright2_objects );
			}

			for ( i=0; i<parameters.spiral_bright3_objects; i++ )
			{
				random = gsl_rng_uniform( r );
				for ( j=0; j<number_of_spirals; j++ )
				{
					k = j;
					if ( random < (gdouble)(j+1)/number_of_spirals ) break;
				}

				angle_random = ( spiral_max[k] - spiral_min[k] ) * gsl_rng_uniform( r );
				radius_random = spiral_radius_parameter[k] * exp( spiral_angle_parameter[k] * angle_random );
				radius_random += fabs( 0.5 * parameters.spiral_radius_randomness * gsl_ran_gaussian( r, ( rmax - rmin ) / 2 ) );
				rotation_random = spiral_rotation[k] + parameters.spiral_angle_randomness * M_PI * gsl_rng_uniform( r );
				x_random = radius_random * ( cos( angle_random )*cos( rotation_random ) + sin( angle_random )*sin( rotation_random ) );
				y_random = radius_random * (-cos( angle_random )*sin( rotation_random ) + sin( angle_random )*cos( rotation_random ) );

				x = x_random + parameters.object_x;
				y = y_random + parameters.object_y;

				z = CLAMP_DOUBLE( MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 2 +
					0.1 * gsl_ran_gaussian( r, parameters.spiral_spiral_thickness * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ), 0., MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );

				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 0.06 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				sigma_x = sigma_random;
				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 2*sigma_x/*0.08 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100*/ ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				if ( sigma_random > sigma_x )
				{
					sigma_y = sigma_random;
				}
				else
				{
					sigma_y = sigma_x;
					sigma_x = sigma_random;
				}
				CLAMP_DOUBLE( sigma_y, sigma_x, 4*sigma_x );

				rotation = atan2( 1., spiral_angle_parameter[k] ) - M_PI/2 + angle_random + gsl_ran_gaussian( r, M_PI/16 );

				brightness = CLAMP_DOUBLE( 255. * ( 3. / pow( ( sigma_x + sigma_y ), 0.8 ) + 0.2 * gsl_rng_uniform( r ) ), 0., 255. );

				color_r = parameters.spiral_color_bright3.r;
				color_g = parameters.spiral_color_bright3.g;
				color_b = parameters.spiral_color_bright3.b;
				color_r += gsl_ran_gaussian( r, parameters.spiral_color_bright3_sigma )/255;
				color_g += gsl_ran_gaussian( r, parameters.spiral_color_bright3_sigma )/255;
				color_b += gsl_ran_gaussian( r, parameters.spiral_color_bright3_sigma )/255;

				object_add( &objects_bright3, objects_bright3_number, i, x, y, z, brightness,
					color_r, color_g, color_b, sigma_x, sigma_y, rotation );

				if ( i%(1+parameters.spiral_bright3_objects/100) ) reporter( "", (gdouble)(i+1)/parameters.spiral_bright3_objects );
			}

			for ( i=0; i<parameters.spiral_dark_objects; i++ )
			{
				random = gsl_rng_uniform( r );
				for ( j=0; j<number_of_spirals; j++ )
				{
					k = j;
					if ( random < (gdouble)(j+1)/number_of_spirals ) break;
				}

				angle_random = ( spiral_max[k] - spiral_min[k] ) * gsl_rng_uniform( r );
				radius_random = spiral_radius_parameter[k] * exp( spiral_angle_parameter[k] * angle_random );
				radius_random += fabs( 0.5 * parameters.spiral_radius_randomness * gsl_ran_gaussian( r, ( rmax - rmin ) / 2 ) );
				rotation_random = spiral_rotation[k] + parameters.spiral_angle_randomness * M_PI * gsl_rng_uniform( r );
				x_random = radius_random * ( cos( angle_random )*cos( rotation_random ) + sin( angle_random )*sin( rotation_random ) );
				y_random = radius_random * (-cos( angle_random )*sin( rotation_random ) + sin( angle_random )*cos( rotation_random ) );

				x = x_random + parameters.object_x;
				y = y_random + parameters.object_y;

				z = CLAMP_DOUBLE( MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 2 +
					0.1 * gsl_ran_gaussian( r, parameters.spiral_spiral_thickness * MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) ), 0., MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );

				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 0.06 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				sigma_x = sigma_random;
				sigma_random = CLAMP_DOUBLE( fabs( 0.03 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100 +
					gsl_ran_gaussian( r, 2*sigma_x/*0.08 * parameters.spiral_spiral_radius *
					MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100*/ ) ), 1.,
						MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) );
				if ( sigma_random > sigma_x )
				{
					sigma_y = sigma_random;
				}
				else
				{
					sigma_y = sigma_x;
					sigma_x = sigma_random;
				}
				CLAMP_DOUBLE( sigma_y, sigma_x, 4*sigma_x );

				rotation = atan2( 1., spiral_angle_parameter[k] ) - M_PI/2 + angle_random + gsl_ran_gaussian( r, M_PI/16 );

				brightness = CLAMP_DOUBLE( 255. * ( 3. / pow( ( sigma_x + sigma_y ), 0.8 ) + 0.2 * gsl_rng_uniform( r ) ), 0., 255. );

				color_r = parameters.spiral_color_dark.r;
				color_g = parameters.spiral_color_dark.g;
				color_b = parameters.spiral_color_dark.b;
				color_r += gsl_ran_gaussian( r, parameters.spiral_color_dark_sigma )/255;
				color_g += gsl_ran_gaussian( r, parameters.spiral_color_dark_sigma )/255;
				color_b += gsl_ran_gaussian( r, parameters.spiral_color_dark_sigma )/255;

				object_add( &objects_dark, objects_dark_number, i, x, y, z, brightness,
					color_r, color_g, color_b, sigma_x, sigma_y, rotation );

				if ( i%(1+parameters.spiral_dark_objects/100) ) reporter( "", (gdouble)(i+1)/parameters.spiral_dark_objects );
			}

			free( spiral_radius_parameter );
			free( spiral_angle_parameter );
			free( spiral_rotation );
			free( spiral_min );
			free( spiral_max );

			break;
		}
		case GALAXY_SPIRALBARRED:
			object_new( &objects_bright1, &objects_bright1_number, parameters.spiralbarred_bright1_objects );
			object_new( &objects_bright2, &objects_bright2_number, parameters.spiralbarred_bright2_objects );
			object_new( &objects_bright3, &objects_bright3_number, parameters.spiralbarred_bright3_objects );
			object_new( &objects_dark, &objects_dark_number, parameters.spiralbarred_dark_objects );

			break;
		case GALAXY_ELLIPTICAL:
			break;
		default:
			break;
	}

	gsl_rng_free( r );

	reporter( "", 0. );

	if ( !cancel_recalculation )
	{
		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Recalculate distribution") );
		recalculating = FALSE;

		gimp_preview_invalidate( GIMP_PREVIEW( preview ) );
	}
}

void draw_bulge_or_ellipse( const gint galaxy_type, gdouble *image_buffer, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, const gboolean show_progress )
{
	gint32 x, y;
	gdouble color_r, color_g, color_b;

	switch ( galaxy_type )
	{
		case GALAXY_SPIRAL:
		case GALAXY_SPIRALBARRED:
		{
			temperature_to_rgb_relative( galaxy_type == GALAXY_SPIRAL ? parameters.spiral_bulge_color_mean :
				parameters.spiralbarred_bulge_color_mean, &color_r, &color_g, &color_b );

			gdouble radius_norm = 0.2 *
				( galaxy_type == GALAXY_SPIRAL ? parameters.spiral_bulge_radius : parameters.spiralbarred_bulge_radius ) *
				MAX( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100.;

			gdouble radius;
			gdouble brightness;

			for ( y=y_start; y<y_size+y_start; y++ )
			{
				for ( x=x_start; x<x_size+x_start; x++ )
				{
					radius = sqrt( ((double)(x-parameters.object_x))*((double)(x-parameters.object_x)) +
						((double)(y-parameters.object_y))*((double)(y-parameters.object_y)) );

					brightness = 255. / ( 1. + pow( radius/radius_norm, 2. ) );

					image_buffer[4*((y-y_start)*x_size+(x-x_start))+0] = CLAMP( brightness*color_r, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+1] = CLAMP( brightness*color_g, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+2] = CLAMP( brightness*color_b, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+3] = 127.5;
				}
			}

			break;
		}
		case GALAXY_ELLIPTICAL:
		{
			temperature_to_rgb_relative( parameters.elliptical_color_mean, &color_r, &color_g, &color_b );

			gdouble radius_norm = 0.1 * parameters.elliptical_radius * MAX( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 100.;

			gdouble b;
			gdouble radius, dradius, phi;
			gdouble brightness;

			for ( y=y_start; y<y_size+y_start; y++ )
			{
				for ( x=x_start; x<x_size+x_start; x++ )
				{
					radius = sqrt( ((double)(x-parameters.object_x))*((double)(x-parameters.object_x)) +
						((double)(y-parameters.object_y))*((double)(y-parameters.object_y)) );

					phi = atan2( (double)(y)-(double)(parameters.object_y), (double)(x)-(double)(parameters.object_x) );
					phi += M_PI * parameters.phi / 180.;

					dradius = -parameters.elliptical_boxiness*cos(4.*phi);
					radius *= 1. + dradius;

					b = radius*sqrt( 1. - parameters.elliptical_excentricity*parameters.elliptical_excentricity*cos(phi)*cos(phi) );

					brightness = 255. / ( 1. + pow( b/radius_norm, 2. ) );

					image_buffer[4*((y-y_start)*x_size+(x-x_start))+0] = CLAMP( brightness*color_r, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+1] = CLAMP( brightness*color_g, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+2] = CLAMP( brightness*color_b, 0, 255 );
					image_buffer[4*((y-y_start)*x_size+(x-x_start))+3] = CLAMP( brightness, 0, 255 );
				}
			}
			break;
		}
		default:
			break;
	}
}

void draw_objects( gdouble *image_buffer, gint32 x_start, gint32 y_start,
	gint32 x_size, gint32 y_size, gdouble **list, const gint32 list_number, const gboolean show_progress )
{
	gdouble x, y, z;
	gdouble x_rotated, y_rotated, z_rotated;
	gdouble x_center, y_center, z_center;
	gdouble phi, psi, theta;
	gdouble brightness, color_r, color_g, color_b;
	gdouble sigma_x, sigma_y, rotation;

	gdouble draw_width;
	gdouble brightness_new;
	gdouble r_new, g_new, b_new, transparency_new;

	gint32 i, x_it, y_it;

	x_center = parameters.object_x;
	y_center = parameters.object_y;
	z_center = MIN( gimp_image_width( image_id ), gimp_image_height( image_id ) ) / 2;

	phi = parameters.phi * M_PI / 180;
	psi = parameters.psi * M_PI / 180;
	theta = parameters.theta * M_PI / 180;

	for ( i=0; i<list_number; i++ )
	{
		object_get( list, i, &x, &y, &z, &brightness, &color_r, &color_g, &color_b, &sigma_x, &sigma_y, &rotation );

		draw_width = 4 * MAX( sigma_x, sigma_y );

/* rotation by euler angles, x convention */
/*
( cos(psi)cos(phi)-sin(psi)cos(theta)sin(phi)  -cos(psi)sin(phi)-sin(psi)cos(theta)cos(phi)  sin(psi)sin(theta) )
( sin(psi)cos(phi)+cos(psi)cos(theta)sin(phi)   cos(psi)cos(theta)cos(phi)-sin(psi)sin(phi) -cos(psi)sin(theta) )
(           sin(theta)sin(phi)                                 sin(theta)cos(phi)               cos(theta)      )
*/

		x_rotated =
			  (cos(psi)*cos(phi)-sin(psi)*cos(theta)*sin(phi))*(x-x_center)
			+ (-cos(psi)*sin(phi)-sin(psi)*cos(theta)*cos(phi))*(y-y_center)
			+ (sin(psi)*sin(theta))*(z-z_center)
			+ x_center;
		y_rotated =
			  (sin(psi)*cos(phi)+cos(psi)*cos(theta)*sin(phi))*(x-x_center)
			+ (cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi))*(y-y_center)
			+ (-cos(psi)*sin(theta))*(z-z_center)
			+ y_center;

		z_rotated =
			  (sin(theta)*sin(phi))*(x-x_center)
			+ (sin(theta)*cos(phi))*(y-y_center)
			+ (cos(theta))*(z-z_center)
			+ z_center;


		if ( (gint32)(x_rotated+draw_width) >= x_start &&
			(gint32)(y_rotated+draw_width) >= y_start &&
			(gint32)(x_rotated-draw_width) < x_start+x_size &&
			(gint32)(y_rotated-draw_width) < y_start+y_size )
		{
			for ( x_it=(gint32)(x_rotated-draw_width); x_it<=(gint32)(x_rotated+draw_width); x_it++ )
			{
				for ( y_it=(gint32)(y_rotated-draw_width); y_it<=(gint32)(y_rotated+draw_width); y_it++ )
				{
					if ( x_it >= x_start &&
						y_it >= y_start &&
						x_it < x_start+x_size &&
						y_it < y_start+y_size )
					{
						brightness_new = gauss_rotation( x_it, y_it, brightness, x_rotated, y_rotated, sigma_x, sigma_y, rotation );

						transparency_new = MAX( z_rotated, image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+3] );

						r_new = MAX( brightness_new*color_r, image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+0] );
						g_new = MAX( brightness_new*color_g, image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+1] );
						b_new = MAX( brightness_new*color_b, image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+2] );

						image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+0] = r_new;
						image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+1] = g_new;
						image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+2] = b_new;
						image_buffer[4*((y_it-y_start)*x_size+(x_it-x_start))+3] = transparency_new;
					}
				}
			}
		}

		if ( show_progress && i % (1+list_number/100) ) gimp_progress_update( (gdouble)(i) / list_number );
	}
}

void normalize( gdouble *image_buffer, gint32 x_size, gint32 y_size, const gboolean show_progress )
{
	gint32 x, y;
	for ( y=0; y<y_size; y++ )
	{
		for ( x=0; x<x_size; x++ )
		{
			image_buffer[4*(y*x_size+x)+0] *= parameters.multiplier;
			image_buffer[4*(y*x_size+x)+1] *= parameters.multiplier;
			image_buffer[4*(y*x_size+x)+2] *= parameters.multiplier;

			image_buffer[4*(y*x_size+x)+3] = 255.;
		}

		if ( show_progress && y % (1+4*y_size/100) )
			gimp_progress_update( (gdouble)(y) / 4*y_size );
	}
}

void create_noise( gdouble *image_buffer, gint32 x_size, gint32 y_size, const gboolean show_progress )
{
	const gsl_rng_type *T;
	gsl_rng *r;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc( T );

	gsl_rng_set( r, parameters.random_seed+30/*time( NULL )*/ );

	gint32 x, y;
	for ( y=0; y<y_size; y++ )
	{
		for ( x=0; x<x_size; x++ )
		{
			image_buffer[4*(y*x_size+x)+0] += gsl_ran_gaussian( r, parameters.noise*sqrt(image_buffer[4*(y*x_size+x)+0])/100 );
			image_buffer[4*(y*x_size+x)+1] += gsl_ran_gaussian( r, parameters.noise*sqrt(image_buffer[4*(y*x_size+x)+1])/100 );
			image_buffer[4*(y*x_size+x)+2] += gsl_ran_gaussian( r, parameters.noise*sqrt(image_buffer[4*(y*x_size+x)+2])/100 );
		}

		if ( show_progress && y % (1+4*y_size/100) )
			gimp_progress_update( (gdouble)(y) / 4*y_size );
	}

	gsl_rng_free( r );
}

/*
Create galaxy into destination region
*/
static void create_galaxy_selection( GimpPixelRgn *region_destination, const gboolean show_progress )
{
	gint32 x, y;
	gint32 i;

	/* write the galaxy into a buffer before actually drawing it */
	gdouble *image_buffer;
	image_buffer = g_malloc( sizeof(gdouble)*4*region_destination->w*region_destination->h );

	reporter( "artificial_galaxy: Creating background...", -1.0 );
	if ( show_progress ) gimp_progress_init( _("Creating background...") );

	/* Background: transparent and black */
	for( i=0; i<4*region_destination->w*region_destination->h; i++ )
	{
		image_buffer[i] = 0.;

		if ( show_progress && i % (1+4*region_destination->w*region_destination->h/100) )
			gimp_progress_update( (gdouble)(i) / 4*region_destination->w*region_destination->h );
	}

	/* bulge / ellipse */
	if ( ( parameters.galaxy_type == GALAXY_SPIRAL && parameters.spiral_bulge_radius > 0. ) ||
		( parameters.galaxy_type == GALAXY_SPIRALBARRED && parameters.spiralbarred_bulge_radius > 0. ) ||
		( parameters.galaxy_type == GALAXY_ELLIPTICAL && parameters.elliptical_radius > 0. ) )
	{
		reporter( "artificial_galaxy: Creating bulge / ellipse...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating bulge / ellipse...") );

		draw_bulge_or_ellipse( parameters.galaxy_type, image_buffer, region_destination->x, region_destination->y,
			region_destination->w, region_destination->h, show_progress );
	}

	/* spiral, bar */
	switch( parameters.galaxy_type )
	{
		case GALAXY_SPIRALBARRED:
		case GALAXY_SPIRAL:
			reporter( "artificial_galaxy: Creating bright objects (1)...", -1.0 );
			if ( show_progress ) gimp_progress_init( _("Creating bright objects (1)...") );
			draw_objects( image_buffer, region_destination->x, region_destination->y,
				region_destination->w, region_destination->h, &objects_bright1, objects_bright1_number, show_progress );

			reporter( "artificial_galaxy: Creating bright objects (2)...", -1.0 );
			if ( show_progress ) gimp_progress_init( _("Creating bright objects (2)...") );
			draw_objects( image_buffer, region_destination->x, region_destination->y,
				region_destination->w, region_destination->h, &objects_bright2, objects_bright2_number, show_progress );

			reporter( "artificial_galaxy: Creating bright objects (3)...", -1.0 );
			if ( show_progress ) gimp_progress_init( _("Creating bright objects (3)...") );
			draw_objects( image_buffer, region_destination->x, region_destination->y,
				region_destination->w, region_destination->h, &objects_bright3, objects_bright3_number, show_progress );

			reporter( "artificial_galaxy: Creating dark objects...", -1.0 );
			if ( show_progress ) gimp_progress_init( _("Creating dark objects...") );
			draw_objects( image_buffer, region_destination->x, region_destination->y,
				region_destination->w, region_destination->h, &objects_dark, objects_dark_number, show_progress );
			break;
		case GALAXY_ELLIPTICAL:
			break;
	}

	/* Normalizing image */
	reporter( "artificial_galaxy: Normalizing image...", -1.0 );
	if ( show_progress ) gimp_progress_init( _("Normalizing image...") );

	normalize( image_buffer, region_destination->w, region_destination->h, show_progress );

	/* Noise */
	if ( parameters.noise > 0 )
	{
		reporter( "artificial_galaxy: Creating noise...", -1.0 );
		if ( show_progress ) gimp_progress_init( _("Creating noise...") );

		create_noise( image_buffer, region_destination->w, region_destination->h, show_progress );
	}

	reporter( "artificial_galaxy: Writing normalized values to buffer...", -1.0 );
	if ( show_progress ) gimp_progress_init( _("Writing normalized values to buffer...") );

	guchar *buffer;
	buffer = g_malloc( sizeof(guchar)*region_destination->bpp*region_destination->w*region_destination->h );
	for ( y=0; y<region_destination->h; y++ )
	{
		for ( x=0; x<region_destination->w; x++ )
		{
			for ( i=0; i<4; i++ )
			{
				buffer[4*(y*region_destination->w+x)+i] = CLAMP( image_buffer[4*(y*region_destination->w+x)+i], 0, 255 );
			}
		}

		if ( show_progress && y % (1+region_destination->h/100) )
			gimp_progress_update( (gdouble)(y) / region_destination->h );
	}

	reporter( "artificial_galaxy: Setting buffer to image...", -1.0 );
	gimp_pixel_rgn_set_rect( region_destination, buffer, region_destination->x, region_destination->y,
		region_destination->w, region_destination->h );

	if ( show_progress ) gimp_progress_end();

	g_free( image_buffer );
	g_free( buffer );
}

/*
create galaxy
*/
static void create_galaxy()
{
	gimp_image_undo_group_start( image_id );

	gint32 layer_destination;
	GimpPixelRgn region_destination;

	layer_destination = gimp_layer_new( image_id, _("Artificial galaxy"), gimp_image_width( image_id ),
		gimp_image_height( image_id ), GIMP_RGBA_IMAGE, 100, GIMP_NORMAL_MODE );

	gimp_image_add_layer( image_id, layer_destination, 0 );

	gimp_pixel_rgn_init( &region_destination, gimp_drawable_get( layer_destination ), 0, 0,
		gimp_drawable_width( layer_destination ), gimp_drawable_height( layer_destination ), TRUE, TRUE );

	create_galaxy_selection( &region_destination, TRUE );

	gimp_drawable_flush( gimp_drawable_get( layer_destination ) );
	gimp_drawable_merge_shadow( layer_destination, TRUE );
	gimp_drawable_update( layer_destination, 0, 0, gimp_drawable_width( layer_destination ),
	gimp_drawable_height( layer_destination ) );

	object_free( &objects_bright1, &objects_bright1_number );
	object_free( &objects_bright2, &objects_bright2_number );
	object_free( &objects_bright3, &objects_bright3_number );
	object_free( &objects_dark, &objects_dark_number );
	object_free( &objects_all, &objects_all_number );

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

	create_galaxy_selection( &region_destination, FALSE );

	gimp_drawable_flush( drawable );
	gimp_drawable_update( drawable->drawable_id, 0, 0, drawable->width, drawable->height );

	/* Read from shadow tiles ("FALSE, TRUE") */
	gimp_pixel_rgn_init( &region_destination, drawable, 0, 0, drawable->width, drawable->height, FALSE, TRUE );
	gimp_drawable_preview_draw_region( GIMP_DRAWABLE_PREVIEW( preview ), &region_destination );
}

static void spin_spiral_bright1_objects( GtkWidget *spin )
{
	parameters.spiral_bright1_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_bright2_objects( GtkWidget *spin )
{
	parameters.spiral_bright2_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_bright3_objects( GtkWidget *spin )
{
	parameters.spiral_bright3_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_dark_objects( GtkWidget *spin )
{
	parameters.spiral_dark_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_bright1_objects( GtkWidget *spin )
{
	parameters.spiralbarred_bright1_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_bright2_objects( GtkWidget *spin )
{
	parameters.spiralbarred_bright2_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_bright3_objects( GtkWidget *spin )
{
	parameters.spiralbarred_bright3_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_dark_objects( GtkWidget *spin )
{
	parameters.spiralbarred_dark_objects = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_color_bright1_sigma( GtkWidget *spin )
{
	parameters.spiral_color_bright1_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_color_bright2_sigma( GtkWidget *spin )
{
	parameters.spiral_color_bright2_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_color_bright3_sigma( GtkWidget *spin )
{
	parameters.spiral_color_bright3_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_color_dark_sigma( GtkWidget *spin )
{
	parameters.spiral_color_dark_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_color_bright1_sigma( GtkWidget *spin )
{
	parameters.spiralbarred_color_bright1_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_color_bright2_sigma( GtkWidget *spin )
{
	parameters.spiralbarred_color_bright2_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_color_bright3_sigma( GtkWidget *spin )
{
	parameters.spiralbarred_color_bright3_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_color_dark_sigma( GtkWidget *spin )
{
	parameters.spiralbarred_color_dark_sigma = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiral_bulge_color_mean( GtkWidget *spin )
{
	parameters.spiral_bulge_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_spiralbarred_bulge_color_mean( GtkWidget *spin )
{
	parameters.spiralbarred_bulge_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_elliptical_color_mean( GtkWidget *spin )
{
	parameters.elliptical_color_mean = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void create_spiral_table()
{
	galaxy_spiral_table = gtk_table_new( 16, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( galaxy_spiral_table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( galaxy_spiral_table ), 2 );
	gtk_widget_show( galaxy_spiral_table );

	GtkWidget *label;
	GtkWidget *button;
	GtkWidget *spin;
	GtkObject *adj;
	GtkWidget *combo;

	label = gtk_label_new( _("Bright objects (1):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_bright1_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_bright1_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (1) color"), 40, 20,
		&parameters.spiral_color_bright1, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiral_color_bright1 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), button, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_color_bright1_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_color_bright1_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Bright objects (2):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_bright2_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_bright2_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 1, 2, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (2) color"), 40, 20,
		&parameters.spiral_color_bright2, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiral_color_bright2 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), button, 0, 1, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 1, 2, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_color_bright2_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_color_bright2_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 2, 3, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Bright objects (3):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 0, 1, 4, 5, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_bright3_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_bright3_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 1, 2, 4, 5, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (3) color"), 40, 20,
		&parameters.spiral_color_bright3, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiral_color_bright3 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), button, 0, 1, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 1, 2, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_color_bright3_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_color_bright3_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 2, 3, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Dark objects:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 0, 1, 6, 7, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_dark_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_dark_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 1, 2, 6, 7, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Dark objects color"), 40, 20,
		&parameters.spiral_color_dark, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiral_color_dark );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), button, 0, 1, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 1, 2, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_color_dark_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_color_dark_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 2, 3, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	combo = gimp_int_combo_box_new(
		_("Logarithmic spiral"), SPIRAL_LOGARITHMIC,
		_("Archimedic spiral"), SPIRAL_ARCHIMEDIC,
		_("Circles"), SPIRAL_CIRCLES,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( combo ), parameters.spiral_spiral_type );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( combo ), parameters.spiral_spiral_type,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.spiral_spiral_type );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), combo, 0, 3, 8, 9, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( combo );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 9,
		_("Bulge radius:"), 185, 75,
		parameters.spiral_bulge_radius, 0.0, 200.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Radius of the bulge, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_bulge_radius );

	label = gtk_label_new( _("Bulge color:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 0, 1, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 1000, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiral_bulge_color_mean );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiral_bulge_color_mean ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), spin, 1, 2, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("K") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiral_table ), label, 2, 3, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 11,
		_("Spiral radius:"), 185, 75,
		parameters.spiral_spiral_radius, 0.0, 200.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Radius of the spiral, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_spiral_radius );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 12,
		_("Spiral thickness:"), 185, 75,
		parameters.spiral_spiral_thickness, 0.0, 50.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Thickness sigma of the spiral, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_spiral_thickness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 13,
		_("Spiral turns:"), 185, 75,
		parameters.spiral_spiral_turns, 0.1, 3.0, 0.1, 5, 1,
		TRUE, 0, 0, _("How often the spirals turn around"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_spiral_turns );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 14,
		_("Spiral randomness:"), 185, 75,
		parameters.spiral_spiral_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the spiral"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_spiral_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 15,
		_("Angle randomness:"), 185, 75,
		parameters.spiral_angle_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the angle"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_angle_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 16,
		_("Radius randomness:"), 185, 75,
		parameters.spiral_radius_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the radius"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_radius_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 17,
		_("Spiral lumpiness:"), 185, 75,
		parameters.spiral_spiral_lumpiness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Lumpiness of the spiral"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_spiral_lumpiness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiral_table ), 0, 18,
		_("Disc lumpiness:"), 185, 75,
		parameters.spiral_disc_lumpiness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Lumpiness of the disc"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiral_disc_lumpiness );
}

static void create_spiralbarred_table()
{
	galaxy_spiralbarred_table = gtk_table_new( 3, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( galaxy_spiralbarred_table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( galaxy_spiralbarred_table ), 2 );
	gtk_widget_show( galaxy_spiralbarred_table );

	GtkWidget *label;
	GtkWidget *button;
	GtkWidget *spin;
	GtkObject *adj;
	GtkWidget *combo;

	label = gtk_label_new( _("Bright objects (1):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_bright1_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_bright1_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (1) color"), 40, 20,
		&parameters.spiralbarred_color_bright1, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiralbarred_color_bright1 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), button, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_color_bright1_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_color_bright1_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Bright objects (2):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_bright2_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_bright2_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 1, 2, 2, 3, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (2) color"), 40, 20,
		&parameters.spiralbarred_color_bright2, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiralbarred_color_bright2 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), button, 0, 1, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 1, 2, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_color_bright2_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_color_bright2_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 2, 3, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Bright objects (3):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 0, 1, 4, 5, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_bright3_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_bright3_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 1, 2, 4, 5, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Bright objects (3) color"), 40, 20,
		&parameters.spiralbarred_color_bright3, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiralbarred_color_bright3 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), button, 0, 1, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 1, 2, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_color_bright3_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_color_bright3_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 2, 3, 5, 6, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	label = gtk_label_new( _("Dark objects:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 0, 1, 6, 7, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 100000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_dark_objects );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_dark_objects ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 1, 2, 6, 7, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	button = gimp_color_button_new( _("Dark objects color"), 40, 20,
		&parameters.spiralbarred_color_dark, GIMP_COLOR_AREA_FLAT );
	g_signal_connect( button, "color_changed", G_CALLBACK( gimp_color_button_get_color ),
		&parameters.spiralbarred_color_dark );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), button, 0, 1, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( button );
	label = gtk_label_new( _("Sigma:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 1, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 1, 2, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 0, 255, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_color_dark_sigma );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_color_dark_sigma ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 2, 3, 7, 8, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );

	combo = gimp_int_combo_box_new(
		_("Logarithmic spiral"), SPIRAL_LOGARITHMIC,
		_("Archimedic spiral"), SPIRAL_ARCHIMEDIC,
		_("Circles"), SPIRAL_CIRCLES,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( combo ), parameters.spiralbarred_spiral_type );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( combo ), parameters.spiralbarred_spiral_type,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.spiralbarred_spiral_type );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), combo, 0, 3, 8, 9, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( combo );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 9,
		_("Bulge radius:"), 185, 75,
		parameters.spiralbarred_bulge_radius, 0.0, 200.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Radius of the bulge, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_bulge_radius );

	label = gtk_label_new( _("Bulge color:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 0, 1, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 1000, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.spiralbarred_bulge_color_mean );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_spiralbarred_bulge_color_mean ), NULL );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), spin, 1, 2, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("K") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_spiralbarred_table ), label, 2, 3, 10, 11, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 11,
		_("Spiral radius:"), 185, 75,
		parameters.spiralbarred_spiral_radius, 0.0, 200.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Radius of the spiral, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_spiral_radius );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 12,
		_("Spiral thickness:"), 185, 75,
		parameters.spiralbarred_spiral_thickness, 0.0, 50.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Thickness sigma of the spiral, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_spiral_thickness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 13,
		_("Spiral turns:"), 185, 75,
		parameters.spiralbarred_spiral_turns, 0.1, 3.0, 0.1, 5, 1,
		TRUE, 0, 0, _("How often the spirals turn around"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_spiral_turns );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 14,
		_("Spiral randomness:"), 185, 75,
		parameters.spiralbarred_spiral_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the spiral"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_spiral_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 15,
		_("Angle randomness:"), 185, 75,
		parameters.spiralbarred_angle_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the angle"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_angle_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 16,
		_("Radius randomness:"), 185, 75,
		parameters.spiralbarred_radius_randomness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Randomness of the radius"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_radius_randomness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 17,
		_("Spiral lumpiness:"), 185, 75,
		parameters.spiralbarred_spiral_lumpiness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Lumpiness of the spiral"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_spiral_lumpiness );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_spiralbarred_table ), 0, 18,
		_("Disc lumpiness:"), 185, 75,
		parameters.spiralbarred_disc_lumpiness, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Lumpiness of the disc"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.spiralbarred_disc_lumpiness );
}

static void create_elliptical_table()
{
	galaxy_elliptical_table = gtk_table_new( 3, 4, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( galaxy_elliptical_table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( galaxy_elliptical_table ), 2 );
	gtk_widget_show( galaxy_elliptical_table );

	GtkWidget *label;
	GtkWidget *spin;
	GtkObject *adj;

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_elliptical_table ), 0, 0,
		_("Radius:"), 185, 75,
		parameters.elliptical_radius, 0.0, 200.0, 0.1, 5, 1,
		TRUE, 0, 0, _("Radius of the elliptical galaxy, measured in percents of the image size"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.elliptical_radius );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_elliptical_table ), 0, 1,
		_("Excentricity:"), 185, 75,
		parameters.elliptical_excentricity, 0.0, 1.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Numerical excentricity of the elliptical galaxy"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.elliptical_excentricity );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( galaxy_elliptical_table ), 0, 2,
		_("Boxy - Disky:"), 185, 75,
		parameters.elliptical_boxiness, -0.1, 0.2, 0.01, 5, 2,
		TRUE, 0, 0, _("a4 parameter, dr=a4*cos(4*t). a4<0 is boxy, a4>0 is disky"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.elliptical_boxiness );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	label = gtk_label_new( _("Color:") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_elliptical_table ), label, 0, 1, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	spin = gtk_spin_button_new_with_range( 1000, 40000, 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.elliptical_color_mean );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_elliptical_color_mean ), NULL );
	g_signal_connect_swapped( spin, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );
	gtk_table_attach( GTK_TABLE( galaxy_elliptical_table ), spin, 1, 2, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( spin );
	label = gtk_label_new( _("K") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_table_attach( GTK_TABLE( galaxy_elliptical_table ), label, 2, 3, 3, 4, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
}

static void show_spiral_settings()
{
	gtk_widget_hide( galaxy_current_table );
	g_object_ref( galaxy_current_table );
	gtk_container_remove( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );

	galaxy_current_table = galaxy_spiral_table;

	gtk_container_add( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );
	gtk_widget_show( galaxy_current_table );
}

static void show_spiralbarred_settings()
{
	gtk_widget_hide( galaxy_current_table );
	g_object_ref( galaxy_current_table );
	gtk_container_remove( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );

	galaxy_current_table = galaxy_spiralbarred_table;

	gtk_container_add( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );
	gtk_widget_show( galaxy_current_table );
}

static void show_elliptic_settings()
{
	gtk_widget_hide( galaxy_current_table );
	g_object_ref( galaxy_current_table );
	gtk_container_remove( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );

	galaxy_current_table = galaxy_elliptical_table;

	gtk_container_add( GTK_CONTAINER( galaxy_specific_frame ), galaxy_current_table );
	gtk_widget_show( galaxy_current_table );
}

static void galaxy_type_selected( GtkWidget *combo_box )
{
	gint active;
	gimp_int_combo_box_get_active( GIMP_INT_COMBO_BOX( combo_box ), &active );
	switch ( active )
	{
		case GALAXY_SPIRAL:
			show_spiral_settings();
			break;
		case GALAXY_SPIRALBARRED:
			show_spiralbarred_settings();
			break;
		case GALAXY_ELLIPTICAL:
			show_elliptic_settings();
			break;
		default:
			break;
	}
}

static void spin_object_x( GtkWidget *spin )
{
	parameters.object_x = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_object_y( GtkWidget *spin )
{
	parameters.object_y = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void spin_random_seed( GtkWidget *spin )
{
	parameters.random_seed = gtk_spin_button_get_value( GTK_SPIN_BUTTON( spin ) );
}

static void split_layers_toggled( GtkWidget *check )
{
	parameters.split_layers = gtk_toggle_button_get_active( GTK_TOGGLE_BUTTON( check ) );
}

static void recalculate_distribution_clicked( GtkWidget *button )
{
	if ( !recalculating )
	{
		recalculating = TRUE;
		cancel_recalculation = FALSE;

		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Cancel recalculation") );

		create_object_distribution();
	}
	else
	{
		recalculating = FALSE;
		cancel_recalculation = TRUE;

		gtk_button_set_label( GTK_BUTTON( recalculate_button ), _("Recalculate distribution") );
	}
}

static void a_clicked( GtkWidget *button )
{
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_number_spin ), 2500 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_spin ), 35. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_brightness_sigma_spin ), 30. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_spin ), 7000 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( background_color_sigma_spin ), 1000 );
// 
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_number_spin ), 10 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_spin ), 195. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_brightness_sigma_spin ), 30. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_spin ), 7800 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_color_sigma_spin ), 500 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_x_spin ), gimp_image_width( image_id )/2 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( object_center_y_spin ), gimp_image_height( image_id )/2 );
// 	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( object_density_combo ), DENSITY_GAUSS );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( GIMP_SCALE_ENTRY_SPINBUTTON( object_radius_adj ) ), 60 );
// 
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_number_spin ), 3 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_spin ), 110. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_brightness_sigma_spin ), 60. );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_spin ), 7500 );
// 	gtk_spin_button_set_value( GTK_SPIN_BUTTON( foreground_color_sigma_spin ), 1000 );
}

static void b_clicked( GtkWidget *button )
{
}

static void c_clicked( GtkWidget *button )
{
}

static void d_clicked( GtkWidget *button )
{
}

static void e_clicked( GtkWidget *button )
{
}

static void f_clicked( GtkWidget *button )
{
}

void reporter( const gchar *string, const gdouble fraction )
{
	if ( string == "" )
	{
		if ( fraction == 0. )
		{
			gtk_progress_bar_set_text( GTK_PROGRESS_BAR( progress_bar ), string );
			gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
			gtk_main_iteration_do( FALSE );
		}
		else if ( fraction > 0. && fraction <= 1. )
		{
			gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
			gtk_main_iteration_do( FALSE );
		}
	}
	else
	{
		if ( fraction >= 0. && fraction <= 1. )
		{
			gtk_progress_bar_set_text( GTK_PROGRESS_BAR( progress_bar ), string );
			gtk_progress_bar_set_fraction( GTK_PROGRESS_BAR( progress_bar ), fraction );
			gtk_main_iteration_do( FALSE );
		}
		else
		{
			printf( "artificial_galaxy: %s\n", string );
		}
	}
}

/*
main dialog
*/
static gint dialog( gint32 image_id, GimpDrawable *drawable )
{
	GtkWidget *dlg;
	GtkWidget *main_hbox;
	GtkWidget *left_vbox;
	GtkWidget *right_vbox;
	GtkWidget *label;
	GtkWidget *button;
	GtkWidget *star_frame;
	GtkWidget *check_box;
	GtkWidget *random_seed;
	GtkWidget *frame;
	GtkWidget *table;
	GtkObject *adj;
	GtkWidget *spin;

  gimp_ui_init( PLUG_IN_NAME, TRUE );

	dlg = gimp_dialog_new( _("Artificial Galaxy"), "astro_artificial_galaxy", NULL, 0,
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

	// GtkWidget *vertical_line = gtk_vseparator_new();
	// gtk_box_pack_start( GTK_BOX( main_hbox ), vertical_line, FALSE, FALSE, 0 );
	// gtk_widget_show( vertical_line );

	right_vbox = gtk_vbox_new( FALSE, 8 );
	gtk_container_set_border_width( GTK_CONTAINER( right_vbox ), 8 );
	gtk_box_pack_start( GTK_BOX( main_hbox ), right_vbox, FALSE, FALSE, 0 );


/* Draw options (left side) */
	preview = gimp_drawable_preview_new( drawable, &parameters.show_preview );
	gtk_box_pack_start( GTK_BOX( left_vbox ), preview, FALSE, FALSE, 0 );
	gtk_widget_show( preview );

	g_signal_connect( preview, "invalidated", G_CALLBACK( preview_callback ), NULL );

	frame = gimp_frame_new( _("Rendering Options") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	gtk_widget_show( frame );

	table = gtk_table_new( 13, 3, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	parameters.object_x = gimp_image_width( image_id )/2;
	parameters.object_y = gimp_image_height( image_id )/2;
	label = gtk_label_new( _("Object center (X, Y):") );
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	spin = gtk_spin_button_new_with_range( 0, gimp_image_width( image_id ), 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.object_x );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_object_x ), NULL );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
	gtk_widget_show( label );
	gtk_table_attach( GTK_TABLE( table ), spin, 1, 2, 0, 1, 0/*GTK_FILL*/, 0, 0, 0 );
	gtk_widget_show( spin );
	spin = gtk_spin_button_new_with_range( 0, gimp_image_height( image_id ), 1 );
	gtk_spin_button_set_value( GTK_SPIN_BUTTON( spin ), parameters.object_y );
	g_signal_connect( spin, "value_changed", G_CALLBACK( spin_object_y ), NULL );
	gtk_table_attach( GTK_TABLE( table ), spin, 2, 3, 0, 1, 0/*GTK_FILL*/, 0, 0, 0 );
	gtk_widget_show( spin );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 2,
		_("Phi:"), 185, 75,
		parameters.phi, 0., 360., 0.1, 5, 1,
		TRUE, 0, 0, _("Rotation inside galaxy plane"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.phi );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

// 	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 3,
// 		_("Psi:"), 185, 75,
// 		parameters.psi, 0., 360., 0.1, 5, 1,
// 		TRUE, 0, 0, _("Rotation towards observer"), NULL );
// 	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
// 		&parameters.psi );
// 	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );
// 
// 	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 4,
// 		_("Theta:"), 185, 75,
// 		parameters.theta, 0., 360., 0.1, 5, 1,
// 		TRUE, 0, 0, _("Rotation around axis towards observer"), NULL );
// 	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
// 		&parameters.theta );
// 	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );
// 
// 	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 5,
// 		_("Transparency of bright clouds:"), 185, 75,
// 		parameters.transparency_bright, 0., 100., 0.1, 5, 1,
// 		TRUE, 0, 0, _("Transparency of bright clouds"), NULL );
// 	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
// 		&parameters.transparency_bright );
// 	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );
// 
// 	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 6,
// 		_("Transparency of dark clouds:"), 185, 75,
// 		parameters.transparency_dark, 0., 100., 0.1, 5, 1,
// 		TRUE, 0, 0, _("Transparency of dark clouds"), NULL );
// 	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
// 		&parameters.transparency_dark );
// 	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 7,
		_("Noise:"), 185, 75,
		parameters.noise, 0, 200, 1, 5, 0,
		TRUE, 0, 0, _("Noise in % of sqrt(N) = photon noise"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_int_adjustment_update ),
		&parameters.noise );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	adj = gimp_scale_entry_new( GTK_TABLE( table ), 0, 8,
		_("Multiplier:"), 185, 75,
		parameters.multiplier, 0.01, 10.0, 0.01, 5, 2,
		TRUE, 0, 0, _("Multiply with value"), NULL );
	g_signal_connect( adj, "value_changed", G_CALLBACK( gimp_double_adjustment_update ),
		&parameters.multiplier );
	g_signal_connect_swapped( adj, "value_changed", G_CALLBACK( gimp_preview_invalidate ), preview );

// 	check_box = gtk_check_button_new_with_label( _("Split bulge, spirals and dark clouds to layers") );
// 	gtk_table_attach_defaults( GTK_TABLE( table ), check_box, 0, 2, 9, 10 );
// 	gtk_toggle_button_set_active( GTK_TOGGLE_BUTTON( check_box ), parameters.split_layers );
// 	g_signal_connect( check_box, "toggled", G_CALLBACK( split_layers_toggled ), NULL );
// 	g_signal_connect_swapped( check_box, "toggled", G_CALLBACK( gimp_preview_invalidate ), preview );
// 	gtk_widget_show( check_box );

// 	star_frame = gtk_frame_new( _("Sample distributions:") );
// 	gtk_box_pack_start( GTK_BOX( left_vbox ), star_frame, FALSE, FALSE, 0 );
// 	gtk_widget_show( star_frame );
// 
// 	table = gtk_table_new( 2, 3, FALSE );
// 	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
// 	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
// 	gtk_container_add( GTK_CONTAINER( star_frame ), table );
// 	gtk_widget_show( table );
// 
// 	button = gtk_button_new_with_label( _("Sample A") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 0, 1, 0, 1, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( a_clicked ), NULL );
// 
// 	button = gtk_button_new_with_label( _("Sample B") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( b_clicked ), NULL );
// 
// 	button = gtk_button_new_with_label( _("Sample C") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 0, 1, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( c_clicked ), NULL );
// 
// 	button = gtk_button_new_with_label( _("Sample D") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 1, 2, 1, 2, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( d_clicked ), NULL );
// 
// 	button = gtk_button_new_with_label( _("Sample E") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 2, 3, 0, 1, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( e_clicked ), NULL );
// 
// 	button = gtk_button_new_with_label( _("Sample F") );
// 	gtk_table_attach( GTK_TABLE( table ), button, 2, 3, 1, 2, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( button );
// 	g_signal_connect( button, "clicked", G_CALLBACK( f_clicked ), NULL );


/* Distribution options (right side) */
	// frame = gimp_frame_new( _("Object distribution:") );
	// gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	// gtk_widget_show( frame );

	table = gtk_table_new( 3, 1, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

// 	recalculate_button = gtk_button_new_with_label( _("Recalculate distribution") );
// 	gtk_table_attach( GTK_TABLE( table ), recalculate_button, 0, 1, 0, 1, GTK_FILL | GTK_EXPAND, 0, 0, 0 );
// 	gtk_widget_show( recalculate_button );
// 	g_signal_connect( recalculate_button, "clicked", G_CALLBACK( recalculate_distribution_clicked ), NULL );
// 
// 	progress_bar = gtk_progress_bar_new();
// 	gtk_widget_set_size_request( progress_bar, -1, 24 );
// 	gtk_table_attach( GTK_TABLE( table ), progress_bar, 0, 1, 1, 2, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( progress_bar );
// 	gtk_progress_bar_set_orientation( GTK_PROGRESS_BAR( progress_bar ), GTK_PROGRESS_LEFT_TO_RIGHT );
// 
// 	random_seed = gimp_random_seed_new( &parameters.random_seed, &parameters.random_seed_bool );
// 	gtk_table_attach( GTK_TABLE( table ), random_seed, 0, 1, 2, 3, GTK_FILL, 0, 0, 0 );
// 	gtk_widget_show( random_seed );
// 	g_signal_connect( GIMP_RANDOM_SEED_SPINBUTTON( random_seed ), "value_changed", G_CALLBACK( spin_random_seed ), NULL );


	// frame = gimp_frame_new( "" );
	// gtk_box_pack_start( GTK_BOX( right_vbox ), frame, FALSE, FALSE, 0 );
	// gtk_widget_show( frame );

	table = gtk_table_new( 2, 1, FALSE );
	gtk_table_set_col_spacings( GTK_TABLE( table ), 6 );
	gtk_table_set_row_spacings( GTK_TABLE( table ), 2 );
	gtk_container_add( GTK_CONTAINER( frame ), table );
	gtk_widget_show( table );

	label = gtk_label_new( _("Galaxy type:") );
	gtk_table_attach( GTK_TABLE( table ), label, 0, 1, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_misc_set_alignment( GTK_MISC( label ), 0, 0.5 );
	gtk_widget_show( label );

	galaxy_type = gimp_int_combo_box_new(
// 		_("Spiral galaxy"), GALAXY_SPIRAL,
// 		_("Spiral barred galaxy"), GALAXY_SPIRALBARRED,
		_("Elliptical galaxy"), GALAXY_ELLIPTICAL,
		NULL );
	gimp_int_combo_box_set_active( GIMP_INT_COMBO_BOX( galaxy_type ), parameters.galaxy_type );
	gimp_int_combo_box_connect( GIMP_INT_COMBO_BOX( galaxy_type ), parameters.galaxy_type,
		G_CALLBACK( gimp_int_combo_box_get_active ), &parameters.galaxy_type );
	g_signal_connect( galaxy_type, "changed", G_CALLBACK( galaxy_type_selected ), NULL );
	g_signal_connect_swapped( galaxy_type, "changed", G_CALLBACK( gimp_preview_invalidate ), preview );

	gtk_table_attach( GTK_TABLE( table ), galaxy_type, 1, 3, 0, 1, GTK_EXPAND | GTK_FILL, 0, 0, 0);
	gtk_widget_show( galaxy_type );

	galaxy_specific_frame = gimp_frame_new( _("Galaxy Settings") );
	gtk_box_pack_start( GTK_BOX( right_vbox ), galaxy_specific_frame, FALSE, FALSE, 0 );
	gtk_widget_show( galaxy_specific_frame );

	create_spiral_table();
	create_spiralbarred_table();
	create_elliptical_table();

	switch ( parameters.galaxy_type )
	{
		case GALAXY_SPIRAL:
			show_spiral_settings();
			break;
		case GALAXY_SPIRALBARRED:
			show_spiralbarred_settings();
			break;
		case GALAXY_ELLIPTICAL:
			show_elliptic_settings();
			break;
		default:
			break;
	}

	gtk_widget_show( left_vbox );
	gtk_widget_show( right_vbox );
	gtk_widget_show( main_hbox );
	gtk_widget_show( dlg );

	gboolean run = ( gimp_dialog_run( GIMP_DIALOG( dlg ) ) == GTK_RESPONSE_OK );

	gtk_widget_destroy( dlg );

	return run;
}
