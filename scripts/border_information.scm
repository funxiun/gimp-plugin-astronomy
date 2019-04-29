;; /***************************************************************************
;;  *   Copyright (C) 2007-2014 by Georg Hennig                               *
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

;; Set your defaults here, sizes in pixels
;; Empty strings as "" (or 0 size).
(define LEFT1 "")
(define LEFT1_SIZE 0)

(define LEFT2 "Object")
(define LEFT2_SIZE 32)

(define MIDDLE1 "Description 1")
(define MIDDLE1_SIZE 13)

(define MIDDLE2 "Description 2")
(define MIDDLE2_SIZE 13)

(define RIGHT1 "Author's Name")
(define RIGHT1_SIZE 12)

(define RIGHT2 "01.12.2014")
(define RIGHT2_SIZE 16)

;; If BORDERx is set to "0", it simply isn't drawn
(define BORDER1 2)
(define BORDER1_COLOR '(255 255 255))

(define BORDER2 12)
(define BORDER2_COLOR '(0 0 0))

(define GAUSS_RADIUS 0)

;; Should be at least max{LEFT1_SIZE+LEFT2_SIZE; MIDDLE1_SIZE+MIDDLE2_SIZE; RIGHT1_SIZE+RIGHT2_SIZE}
;; Otherwise text will be drawn on the BORDER2, over each other or even outside the image.
(define INFORMATION_HEIGHT 42)

(define FONT "Century Gothic")
(define FONT_COLOR '(255 255 255))

(define (astronomy-border-information image drawable inner_size inner_color
					outer_size outer_color gauss_radius information_height font font_color
					text_left1 font_size_left1 text_left2 font_size_left2
					text_middle1 font_size_middle1 text_middle2 font_size_middle2
					text_right1 font_size_right1 text_right2 font_size_right2)
	(gimp-undo-push-group-start image)

	;; font_size_** may not be zero; assume that user doesn't want to display any text in this case
	(if (<= font_size_left1 0) (begin (set! font_size_left1 1) (set! text_left1 "")))
	(if (<= font_size_left2 0) (begin (set! font_size_left2 1) (set! text_left2 "")))
	(if (<= font_size_middle1 0) (begin (set! font_size_middle1 1) (set! text_middle1 "")))
	(if (<= font_size_middle2 0) (begin (set! font_size_middle2 1) (set! text_middle2 "")))
	(if (<= font_size_right1 0) (begin (set! font_size_right1 1) (set! text_right1 "")))
	(if (<= font_size_right2 0) (begin (set! font_size_right2 1) (set! text_right2 "")))

	;; Add alpha channel to all layers, so we can move our information layer to the bottom
	(define layers_number 0)
	(let* ((layers (gimp-image-get-layers image))
		(layers_num (car layers))
		(layers_array (cadr layers))
		(i (- layers_num 1)))
		(set! layers_number layers_num)
		(while (>= i 0)
			(gimp-layer-add-alpha (aref layers_array i))
			(set! i (- i 1))
		)
	)

	;; Backup foreground color
	(define color_tmp (car (gimp-context-get-foreground)))

	;; Original Image width and height
	(define image_width (car (gimp-image-width image)))
	(define image_height (car (gimp-image-height image)))

	;; Add border 1 as layer
	(gimp-image-resize
		image
		(+ image_width (* 2 inner_size))
		(+ image_height (* 2 inner_size))
		inner_size
		inner_size
	)
	(define border1_layer
		(car (gimp-layer-new
			image
			(car (gimp-image-width image))
			(car (gimp-image-height image))
			RGBA-IMAGE
			"border 1"
			100.0
			NORMAL-MODE
		))
	)
	(gimp-image-add-layer image border1_layer layers_number)
	(gimp-image-set-active-layer image border1_layer)
	(gimp-context-set-foreground inner_color)
	(gimp-drawable-fill (car (gimp-image-get-active-drawable image)) FOREGROUND-FILL)

	;; Add border 2 as layer
	(gimp-image-resize
		image
		(+ (car (gimp-image-width image)) (* 2 outer_size))
		(+ (+ (car (gimp-image-height image)) (* 2 outer_size)) information_height)
		outer_size
		outer_size
	)
	(define border2_layer
		(car (gimp-layer-new
			image
			(car (gimp-image-width image))
			(car (gimp-image-height image))
			RGBA-IMAGE
			"Border + Information"
			100.0
			NORMAL-MODE
		))
	)
	(gimp-image-add-layer image border2_layer (+ 1 layers_number))
	(gimp-image-set-active-layer image border2_layer)
	(gimp-context-set-foreground outer_color)
	(gimp-drawable-fill (car (gimp-image-get-active-drawable image)) FOREGROUND-FILL)

	;; Text
	(gimp-context-set-foreground font_color)

	(define middle_width (/ (car (gimp-image-width image)) 2))

	;; Left text (aligned left on the left of the inner border)
	(define text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_left1 ;; text
			font_size_left1 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(define text_x (max 0 outer_size))
	;; inner size as spacer between inner border and text (therefore "* 3", not "* 2").
	(define text_y (max 0 (+ image_height (+ outer_size (* 3 inner_size)))))
	(define text1_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_left1 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_left1 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))
	(set! text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_left2 ;; text
			font_size_left2 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(set! text_x (max 0 outer_size))
	(set! text_y (max 0 (- (car (gimp-image-height image)) (+ outer_size (cadr text_extents)))))
	(define text2_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_left2 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_left2 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))

	;; Text in the middle
	(set! text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_middle1 ;; text
			font_size_middle1 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(set! text_x (max 0 (- middle_width (/ (car text_extents) 2))))
	(set! text_y (max 0 (+ image_height (+ outer_size (* 3 inner_size)))))
	(define text3_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_middle1 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_middle1 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))
	(set! text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_middle2 ;; text
			font_size_middle2 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(set! text_x (max 0 (- middle_width (/ (car text_extents) 2))))
	(set! text_y (max 0 (- (car (gimp-image-height image)) (+ outer_size (cadr text_extents)))))
	(define text4_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_middle2 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_middle2 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))

	;; Right text
	(set! text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_right1 ;; text
			font_size_right1 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(set! text_x (max (- (- (car (gimp-image-width image)) outer_size) (car text_extents))))
	(set! text_y (max 0 (+ image_height (+ outer_size (* 3 inner_size)))))
	(define text5_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_right1 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_right1 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))
	(set! text_extents ;; contains width, height, ascent, descent afterwards.
		(gimp-text-get-extents-fontname
			text_right2 ;; text
			font_size_right2 ;; size
			PIXELS ;; alternatively: POINTS
			font ;; font name
		)
	)
	(set! text_x (max (- (- (car (gimp-image-width image)) outer_size) (car text_extents))))
	(set! text_y (max 0 (- (car (gimp-image-height image)) (+ outer_size (cadr text_extents)))))
	(define text6_layer (car (gimp-text-fontname
		image ;; image
		-1 ;; drawable; -1 creates new text layer
		text_x ;; x position
		text_y ;; y position
		text_right2 ;; text
		-1 ;; border
		TRUE ;; antialias
		font_size_right2 ;; font size
		PIXELS ;; alternatively: POINTS
		font ;; font name
	)))

	;; Move it to the bottom
	(gimp-image-lower-layer-to-bottom image border1_layer)
	(gimp-image-lower-layer-to-bottom image border2_layer)

	;; Merge all layers into one (only merge text layers if they were created;
	;; empty text layers are not created. merge the resulting layer (might be border1_layer) down.
	(define merge_layer (car (gimp-image-merge-down image border1_layer EXPAND-AS-NECESSARY)))

	;; smooth inner border -> outer border
	(if (> gauss_radius 0) (plug-in-gauss 1 image merge_layer gauss_radius gauss_radius 0))

	(if (> text1_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text1_layer EXPAND-AS-NECESSARY))))
	(if (> text2_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text2_layer EXPAND-AS-NECESSARY))))
	(if (> text3_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text3_layer EXPAND-AS-NECESSARY))))
	(if (> text4_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text4_layer EXPAND-AS-NECESSARY))))
	(if (> text5_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text5_layer EXPAND-AS-NECESSARY))))
	(if (> text6_layer -1)
		(set! merge_layer (car (gimp-image-merge-down image text6_layer EXPAND-AS-NECESSARY))))

	;; Set foreground color to previous value
	(gimp-context-set-foreground color_tmp)

	(gimp-undo-push-group-end image)

	(gimp-displays-flush)
)

(script-fu-register "astronomy-border-information"
	_"Draw Border with Image Information"
	_"Draw a border around an image and some information about it image in the bottom"
	"Georg Hennig"
	"Georg Hennig"
	"12-2014"
	"RGB*"
	SF-IMAGE    "SF-IMAGE" 0
	SF-DRAWABLE "SF-DRAWABLE" 0
	SF-ADJUSTMENT _"Inner border size" (list BORDER1 0 100 1 10 0 1)
	SF-COLOR _"Inner border color" BORDER1_COLOR
	SF-ADJUSTMENT _"Outer border size" (list BORDER2 0 100 1 10 0 1)
	SF-COLOR _"Outer border color" BORDER2_COLOR
	SF-ADJUSTMENT _"Gauss Radius" (list GAUSS_RADIUS 0 20 0.1 1 0 1)
	SF-ADJUSTMENT _"Bottom information height" (list INFORMATION_HEIGHT 0 200 1 10 0 1)
	SF-FONT _"Font" FONT
	SF-COLOR _"Font Color" FONT_COLOR
	SF-STRING _"Text (left, first line)" LEFT1
	SF-ADJUSTMENT _"Font size (left, first line)" (list LEFT1_SIZE 0 100 1 10 0 1)
	SF-STRING _"Text (left, second line)" LEFT2
	SF-ADJUSTMENT _"Font size (left, second line)" (list LEFT2_SIZE 0 100 1 10 0 1)
	SF-STRING _"Text (middle, first line)" MIDDLE1
	SF-ADJUSTMENT _"Font size (middle, first line)" (list MIDDLE1_SIZE 0 100 1 10 0 1)
	SF-STRING _"Text (middle, second line)" MIDDLE2
	SF-ADJUSTMENT _"Font size (middle, second line)" (list MIDDLE2_SIZE 0 100 1 10 0 1)
	SF-STRING _"Text (right, first line)" RIGHT1
	SF-ADJUSTMENT _"Font size (right, first line)" (list RIGHT1_SIZE 0 100 1 10 0 1)
	SF-STRING _"Text (right, second line)" RIGHT2
	SF-ADJUSTMENT _"Font size (right, second line)" (list RIGHT2_SIZE 0 100 1 10 0 1)
)

(script-fu-menu-register "astronomy-border-information"
                         "<Image>/Filters/Astronomy")