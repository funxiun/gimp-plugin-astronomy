;; /***************************************************************************
;;  *   Copyright (C) 2006-2014 by Georg Hennig                               *
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

(define (astronomy-brightness-contrast-batch image drawable visible_only brightness contrast)
	(gimp-undo-push-group-start image)

	(let* ((layers (gimp-image-get-layers image))
		(layers_num (car layers))
		(layers_array (cadr layers))
		(i (- layers_num 1)))
		(while (>= i 0)
			(if (or (= visible_only FALSE) (= (car (gimp-drawable-get-visible (aref layers_array i))) TRUE))
				(gimp-brightness-contrast (aref layers_array i) brightness contrast))
			(set! i (- i 1))
		)
	)

	(gimp-undo-push-group-end image)

	(gimp-displays-flush)
)

(script-fu-register "astronomy-brightness-contrast-batch"
	_"Set Brightness and Contrast for all Layers"
	_"Apply the Brightness/Contrast tool to all layers in an image"
	"Georg Hennig"
	"Georg Hennig"
	"12-2014"
	""
	SF-IMAGE    "SF-IMAGE" 0
	SF-DRAWABLE "SF-DRAWABLE" 0
	SF-TOGGLE   _"Visible layers only" FALSE
	SF-ADJUSTMENT _"Brightness" (list 0 -127 127 1 10 0 0)
	SF-ADJUSTMENT _"Contrast" (list 0 -127 127 1 10 0 0)
)

(script-fu-menu-register "astronomy-brightness-contrast-batch"
                         "<Image>/Filters/Astronomy")
