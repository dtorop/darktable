/*
    This file is part of darktable,
    Copyright (C) 2021 darktable developers.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once

#include <glib.h>
#include <inttypes.h>
#include "common/pdf.h"
#include "common/cups_print.h"
#include "common/image.h"
#include "common/math.h"

// FIXME: increase this #
#define MAX_IMAGE_PER_PAGE 20

// FIXME: should these be double?
typedef struct _image_pos
{
  float x, y, width, height;
} dt_image_pos;

typedef struct _image_box
{
  dt_imgid_t imgid;
  // FIXME: should this be float?
  int32_t max_width, max_height; // max size for the export (in page pixels)
  int32_t exp_width, exp_height; // final exported size (in pixels)
  int32_t dis_width, dis_height; // image size on screen (in pixels)
  int32_t img_width, img_height; // the final image size as it will be exported
  dt_alignment_t alignment;
  // FIXME: do we need all this dimensions data?
  // FIXME: is this true?
  // while we track the image in screen/print location, the ground
  // truth for its location is relative to the page dimensions, so
  // that if page dimensions change, the layout will remain plausible
  // FIXME: s/pos/rel_pos/
  dt_image_pos pos;              // relative pos from screen.page
  // FIXME: if screen position is always calculated from pos, then just use that
  dt_image_pos screen;           // current screen pos (in pixels)
  dt_image_pos print;            // current print pos (in pixels) depending on paper size + DPI
  uint16_t *buf;
} dt_image_box;

typedef struct dt_screen_pos
{
  dt_image_pos page;       // this is for reference and is the box of the
                           // white page (in pixels) in print module.
                           // it is the full page.

  dt_image_pos print_area; // this is for reference and is the box of the
                           // grey area in the white page (in pixels) in print
                           // module. it is the print area (without margins).
  gboolean borderless;     // whether the print is borderless (user's margins below
                           // hardware margins.
} dt_screen_pos;

typedef struct dt_images_box
{
  dt_imgid_t imgid_to_load;
  int32_t motion_over;
  int count;
  // FIXME: this should be a list so that there isn't a hard limit
  dt_image_box box[MAX_IMAGE_PER_PAGE];
  // FIXME: We could only keep track of measurements in page pixels as this is the ground truth for printing, and from that (and page dpi which could be kept in dt_print_t) we can then on the fly calculate dimensions in mm/cm/in for overlay and right panel. Then the other ground truth is proportional measurements of page in GUI from which we figure out page position in pixels.
  // FIXME: this can be calculated in layout_page, not needed here
  float page_width_px, page_height_px;  // full print page in page pixels
  // FIXME: this overlaps with dt_print_info_t->dt_paper_info_t!
  float page_width_mm, page_height_mm; // full print page in mm
  dt_screen_pos screen;
} dt_images_box;

// return the box index or -1 if (x, y) coordinate is not over an image
int32_t dt_printing_get_image_box(const dt_images_box *imgs, const int x, const int y);

void dt_printing_clear_box(dt_image_box *img);
void dt_printing_clear_boxes(dt_images_box *imgs);

/* (x, y) -> (width, height) are in pixels (on screen position) */
void dt_printing_setup_display(dt_images_box *imgs,
                               const float px, const float py, const float pwidth, const float pheight,
                               const float ax, const float ay, const float awidth, const float aheight, gboolean borderless);

void dt_printing_setup_box(dt_images_box *imgs, const int idx,
                           const float x, const float y,
                           const float width, const float height);

/* page_width page_height in mm, compute the max_width and max_height in
   pixels for the image */
void dt_printing_setup_page(dt_images_box *imgs,
                            const float page_width_mm, const float page_height_mm,
                            const int resolution);

/* setup the image id and exported width x height */
void dt_printing_setup_image(dt_images_box *imgs, const int idx,
                             const dt_imgid_t imgid, const int32_t width, const int32_t height,
                             const dt_alignment_t alignment);

/* return the on screen pos with alignment */
void dt_printing_get_screen_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos);
void dt_printing_get_screen_rel_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos);
void dt_printing_get_image_pos_mm(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos);
void dt_printing_get_image_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos);

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on

