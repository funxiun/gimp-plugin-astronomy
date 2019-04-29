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

(define mode_name_list (list _"Normal" _"Addition"
 _"Subtract" _"Difference"
 _"Multiply" _"Divide"
 _"Dodge" _"Burn" _"Screen" _"Overlay"
 _"Hard Light" _"Soft Light"
 _"Darken Only" _"Lighten Only"
 _"Grain extract" _"Grain merge"
 _"Hue" _"Color" _"Saturation" _"Value"))

(define mode_list (list NORMAL-MODE ADDITION-MODE
 SUBTRACT-MODE DIFFERENCE-MODE
 MULTIPLY-MODE DIVIDE-MODE
 DODGE-MODE BURN-MODE SCREEN-MODE OVERLAY-MODE
 HARDLIGHT-MODE SOFTLIGHT-MODE
 DARKEN-ONLY-MODE LIGHTEN-ONLY-MODE
 GRAIN-EXTRACT-MODE GRAIN-MERGE-MODE
 HUE-MODE COLOR-MODE SATURATION-MODE VALUE-MODE))

(define (astronomy-mode-batch image drawable visible_only mode_number opacity)
	(gimp-undo-push-group-start image)

	(let* ((layers (gimp-image-get-layers image))
		(layers_num (car layers))
		(layers_array (cadr layers))
		(i (- layers_num 1))
		(mode_current (nth mode_number mode_list)))
		(while (>= i 0)
			(if (or (= visible_only FALSE) (= (car (gimp-drawable-get-visible (aref layers_array i))) TRUE))
				(begin (gimp-layer-set-opacity (aref layers_array i) opacity)
				(gimp-layer-set-mode (aref layers_array i) mode_current))
			)
			(set! i (- i 1))
		)
	)

	(gimp-undo-push-group-end image)

	(gimp-displays-flush)
)

(script-fu-register "astronomy-mode-batch"
	_"Set Mode for all Layers"
	_"Set a particular blending mode and opacity to all layers"
	"Georg Hennig"
	"Georg Hennig"
	"12-2014"
	""
	SF-IMAGE    "SF-IMAGE" 0
	SF-DRAWABLE "SF-DRAWABLE" 0
	SF-TOGGLE   _"Visible layers only" FALSE
	SF-OPTION   _"Mode" mode_name_list
	SF-ADJUSTMENT _"Opacity" (list 100.0 0.0 100.0 0.1 10 1 0)
)

(script-fu-menu-register "astronomy-mode-batch"
                         "<Image>/Filters/Astronomy")
