## Process this file with automake to produce Makefile.in

documentationdir = $(GIMP_DATADIR)/doc

documentation_DATA = \
	documentation.tex \
	documentation.pdf

documentation.pdf: documentation.tex
	pdflatex documentation.tex

EXTRA_DIST = \
	documentation.pdf \
	documentation.tex \
	align.jpg \
	border_information.jpg \
	star_rounding_drawing.png \
	artificial_galaxy.jpg \
	star_rounding.jpg \
	artificial_stars.jpg \
	star_rounding_long.png \
	background_gradient.jpg \
	merge.jpg \
	star_rounding_rounded.png