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

(define i)
(define mode_tmp)
(define opacity_tmp)
(define layer_new)

(define (astronomy-flat-division image drawable visible_only filename)
	(gimp-undo-push-group-start image)

	(let* ((layers (gimp-image-get-layers image))
		(layers_num (car layers))
		(layers_array (cadr layers))
		(i (- layers_num 1))
		(darkframe (car (gimp-file-load RUN-NONINTERACTIVE filename filename))))
		(while (>= i 0)
			(if (or (= visible_only FALSE) (= (car (gimp-drawable-get-visible (aref layers_array i))) TRUE))
				(let* ((darklayer (car (gimp-layer-new-from-drawable (car (gimp-image-get-active-layer darkframe)) image))))
					(set! mode_tmp (car (gimp-layer-get-mode (aref layers_array i))))
					(set! opacity_tmp (car (gimp-layer-get-opacity (aref layers_array i))))
					(gimp-image-add-layer image darklayer i)
					(gimp-layer-set-mode darklayer DIVIDE-MODE)
					(set! layer_new (car (gimp-image-merge-down image darklayer CLIP-TO-BOTTOM-LAYER)))
					(gimp-layer-set-mode layer_new mode_tmp)
					(gimp-layer-set-opacity layer_new opacity_tmp)
				)
			)
			(set! i (- i 1))
		)
	)

	(gimp-undo-push-group-end image)

	(gimp-displays-flush)
)


(script-fu-register "astronomy-flat-division"
	_"Divide all Layers by Flat Field"
	_"Divides all layers by a flat field"
	"Georg Hennig"
	"Georg Hennig"
	"12-2014"
	""
	SF-IMAGE    "SF-IMAGE" 0
	SF-DRAWABLE "SF-DRAWABLE" 0
	SF-TOGGLE   _"Visible layers only" FALSE
	SF-FILENAME "Darkframe to subtract" "./"
)

(script-fu-menu-register "astronomy-flat-division"
                         "<Image>/Filters/Astronomy")