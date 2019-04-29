;; /***************************************************************************
;;  *   Copyright (C) 2008-2014 by Georg Hennig                               *
;;  *   georg.hennig@web.de                                                   *
;;  *                                                                         *
;;  *   This program is free software; you can redistribute it and/or modify  *
;;  *   it under the terms of the GNU General Public License as published by  *
;;  *   the Free Software Foundation; either version 2 of the License, or     *
;;  *   (at your option) any later version.                                   *
;;  *                                                                         *
;;  *   This program is distributed in the hope that it will be useful,       *
;;  *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
;;  *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
;;  *   GNU General Public License for more details.                          *
;;  *                                                                         *
;;  *   You should have received a copy of the GNU General Public License     *
;;  *   along with this program; if not, write to the                         *
;;  *   Free Software Foundation, Inc.,                                       *
;;  *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
;;  ***************************************************************************/

(define (background-gradient-batch pattern boxa boxb iterations sigma minimum weighted_values percentage)
	(let* ((filelist (cadr (file-glob pattern 1))))
		(while (not (null? filelist))
			(let* ((filename (car filelist))
						(image (car (gimp-file-load RUN-NONINTERACTIVE filename filename)))
						(drawable (car (gimp-image-get-active-layer image))))
				(plug-in-astro-background-gradient RUN-NONINTERACTIVE image drawable boxa boxb iterations sigma minimum weighted_values percentage)
				(set! drawable (car (gimp-image-merge-visible-layers image CLIP-TO-IMAGE)))
				(gimp-file-save RUN-NONINTERACTIVE image drawable filename filename)
				(gimp-image-delete image)
			)
			(set! filelist (cdr filelist))
		)
	)
)
