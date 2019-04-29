## Process this file with automake to produce Makefile.in

bin_PROGRAMS = astronomy-alignment astronomy-merge astronomy-background-gradient astronomy-artificial-stars astronomy-artificial-galaxy astronomy-star-rounding

bindir = $(GIMP_LIBDIR)/plug-ins

astronomy_alignment_SOURCES = \
	alignment.c

astronomy_merge_SOURCES = \
	merge.c

astronomy_background_gradient_SOURCES = \
	background_gradient.c

astronomy_artificial_stars_SOURCES = \
	artificial_stars.c

astronomy_artificial_galaxy_SOURCES = \
	artificial_galaxy.c

astronomy_star_rounding_SOURCES = \
	star_rounding.c

AM_CPPFLAGS = \
	-DLOCALEDIR=\""$(LOCALEDIR)"\"		\
	-DDATADIR=\""$(DATADIR)"\"

INCLUDES =\
	-I$(top_srcdir)		\
	@GIMP_CFLAGS@		\
	-I$(includedir)

LDADD = $(GIMP_LIBS) -lgsl -lgslcblas -lfftw3 -lm

EXTRA_DIST = \
	plugin-intl.h \
	artificial_stars_temperature.h
