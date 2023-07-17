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

#include "common/printing.h"

void _clear_pos(dt_image_pos *pos)
{
  pos->x = pos->y = pos->width = pos->height = 0.0f;
}

void dt_printing_clear_box(dt_image_box *img)
{
  img->imgid = NO_IMGID;
  img->max_width = img->max_height = 0;
  img->exp_width = img->exp_height = 0;
  img->dis_width = img->dis_height = 0;
  img->img_width = img->img_height = 0;
  img->alignment = ALIGNMENT_CENTER;
  img->buf = NULL;
  img->w_box = NULL;

  _clear_pos(&img->screen);
  _clear_pos(&img->pos);
  _clear_pos(&img->print);
}

void dt_printing_clear_boxes(dt_images_box *imgs)
{
  for(int k=0; k<MAX_IMAGE_PER_PAGE; k++)
    dt_printing_clear_box(&imgs->box[k]);

  _clear_pos(&imgs->screen.print_area);

  imgs->count = 0;
  imgs->motion_over = -1;
  imgs->page_width = imgs->page_height = 0;
  imgs->page_width_mm = imgs->page_height_mm = 0;
}

int32_t dt_printing_get_image_box(const dt_images_box *imgs, const int x, const int y)
{
  int box = -1;
  float dist = FLT_MAX;

  for(int k=0; k<imgs->count; k++)
  {
    const dt_image_box *b = &imgs->box[k];

    const float x1 = b->screen.x;
    const float x2 = b->screen.x + b->screen.width;
    const float y1 = b->screen.y;
    const float y2 = b->screen.y + b->screen.height;

    // check if over a box
    if(x > x1 && x < x2 && y > y1 && y < y2)
    {
      // compute min dist
      float dd = sqrf(x1 - x);
      dd = fminf(dd, sqrf(x2 - x));
      dd = fminf(dd, sqrf(y1 - y));
      dd = fminf(dd, sqrf(y2 - y));

      if(dd < dist)
      {
        box = k;
        dist = dd;
      }
    }
  }

  return box;
}

void _compute_rel_pos(const dt_images_box *imgs, const dt_image_pos *ref, dt_image_pos *pos)
{
  // compute the printing position & width as % of the page

  const float page_width  = imgs->screen.page_width;
  const float page_height = imgs->screen.page_height;

  pos->x      = ref->x / page_width;
  pos->y      = ref->y / page_height;
  pos->width  = ref->width / page_width;
  pos->height = ref->height / page_height;
}

void dt_printing_setup_display(dt_images_box *imgs,
                               const float ax, const float ay, const float awidth, const float aheight,
                               gboolean borderless)
{
  const float pwidth = imgs->screen.page_width;
  const float pheight = imgs->screen.page_height;

  imgs->screen.print_area.x      = ax;
  imgs->screen.print_area.y      = ay;
  imgs->screen.print_area.width  = awidth;
  imgs->screen.print_area.height = aheight;

  dt_print(DT_DEBUG_PRINT, "[printing] screen/page  (0.0, 0.0) -> (%3.1f, %3.1f)\n",
           pwidth, pheight);
  dt_print(DT_DEBUG_PRINT, "[printing] screen/parea (%3.1f, %3.1f) -> (%3.1f, %3.1f)\n",
           ax, ay, awidth, aheight);

  // FIXME: store this in pinfo->page.borderless or such if it is actually used by cups_print
  imgs->screen.borderless = borderless;

  // and now reset the box to be resized accordingly if needed
  for(int k=0; k<imgs->count; k++)
  {
    dt_image_box *box = &imgs->box[k];

    if(box->pos.x > 0)
    {
      box->screen.x      = pwidth * box->pos.x;
      box->screen.y      = pheight * box->pos.y;
      box->screen.width  = pwidth * box->pos.width;
      box->screen.height = pheight * box->pos.height;
    }
  }
}

void _contain_to_bounds(const float request_low, const float request_size,
                        const float bound_low, const float bound_size,
                        float *out_low, float *out_size)
{
  // apply bounds, with a sane miminum
  *out_low = fmaxf(request_low, bound_low);
  *out_size = CLAMPF(request_size, 100.0f, bound_size);

  // if past bounds, shift back in, best case without shrinking
  // FIXME: would it be better if the widget just wouldn't let itself be dragged outside of print area
  const float out_high = *out_low + *out_size;
  const float bound_high = bound_low + bound_size;
  if(out_high > bound_high)
  {
    const float off = (out_high - bound_high);
    *out_low = fmaxf(bound_low, *out_low - off);
  }
}

// for the given layout box at idx, set its dimensions in screen and
// relative dimensions, while not exceeding page body
void dt_printing_setup_box(dt_images_box *imgs, const int idx,
                           const float screen_x, const float screen_y,
                           const float screen_width, const float screen_height)
{
  dt_image_box *box = &imgs->box[idx];
  // FIXME: would it be better if the widget just wouldn't let itself be dragged outside of print area?
  _contain_to_bounds(screen_x, screen_width,
                     imgs->screen.print_area.x, imgs->screen.print_area.width,
                     &box->screen.x, &box->screen.width);
  _contain_to_bounds(screen_y, screen_height,
                     imgs->screen.print_area.y, imgs->screen.print_area.height,
                     &box->screen.y, &box->screen.height);

  // FIXME: should this be relative to printable area rather than page dimensions, so image boxes never go outside of printable area?
  _compute_rel_pos(imgs, &box->screen, &box->pos);

  // add box if doesn't exist
  // FIXME: if boxed are a GList will need to allocate before setting bounds
  if(idx == imgs->count) imgs->count++;
}

void dt_printing_setup_page(dt_images_box *imgs,
                            const float page_width, const float page_height,
                            const int resolution)
{
  imgs->page_width_mm = page_width;
  imgs->page_height_mm = page_height;
  imgs->page_width  = dt_pdf_point_to_pixel(dt_pdf_mm_to_point(page_width), resolution);
  imgs->page_height = dt_pdf_point_to_pixel(dt_pdf_mm_to_point(page_height), resolution);

  for(int k=0; k<imgs->count; k++)
  {
    dt_image_box *box = &imgs->box[k];

    box->max_width = box->pos.width * imgs->page_width;
    box->max_height = box->pos.height * imgs->page_height;
  }
}

void _align_pos(const dt_image_pos *ref, const dt_alignment_t alignment,
                const int32_t width, const int32_t height, dt_image_pos *pos)
{
  pos->width  = width;
  pos->height = height;

  switch(alignment)
  {
    case ALIGNMENT_TOP_LEFT:
      pos->x = ref->x;
      pos->y = ref->y;
      break;
    case ALIGNMENT_TOP:
      pos->x = ref->x + (ref->width - width) / 2;
      pos->y = ref->y;
      break;
    case ALIGNMENT_TOP_RIGHT:
      pos->x = ref->x + (ref->width - width);
      pos->y = ref->y;
      break;
    case ALIGNMENT_LEFT:
      pos->x = ref->x;
      pos->y = ref->y + (ref->height - height) / 2;
      break;
    case ALIGNMENT_CENTER:
      pos->x = ref->x + (ref->width - width) / 2;
      pos->y = ref->y + (ref->height - height) / 2;
      break;
    case ALIGNMENT_RIGHT:
      pos->x = ref->x + (ref->width - width);
      pos->y = ref->y + (ref->height - height) / 2;
      break;
    case ALIGNMENT_BOTTOM_LEFT:
      pos->x = ref->x;
      pos->y = ref->y + (ref->height - height);
      break;
    case ALIGNMENT_BOTTOM:
      pos->x = ref->x + (ref->width - width) / 2;
      pos->y = ref->y + (ref->height - height);
      break;
    case ALIGNMENT_BOTTOM_RIGHT:
      pos->x = ref->x + (ref->width - width);
      pos->y = ref->y + (ref->height - height);
      break;
  }
}

void dt_printing_get_screen_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos)
{
  _clear_pos(pos);

  _align_pos(&img->screen, img->alignment, img->dis_width, img->dis_height, pos);
}

void dt_printing_get_screen_rel_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos)
{
  dt_image_pos screen_pos;

  dt_printing_get_screen_pos(imgs, img, &screen_pos);

  _compute_rel_pos(imgs, &screen_pos, pos);
}

void dt_printing_get_image_pos_mm(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos)
{
  dt_image_pos rpos;

  dt_printing_get_screen_rel_pos(imgs, img, &rpos);

  pos->x      = rpos.x * imgs->page_width_mm;
  pos->y      = rpos.y * imgs->page_height_mm;
  pos->width  = rpos.width * imgs->page_width_mm;
  pos->height = rpos.height * imgs->page_height_mm;
}

void dt_printing_get_image_pos(const dt_images_box *imgs, const dt_image_box *img, dt_image_pos *pos)
{
  dt_image_pos rpos;

  dt_printing_get_screen_rel_pos(imgs, img, &rpos);

  pos->x      = rpos.x * imgs->page_width;
  pos->y      = rpos.y * imgs->page_height;
  pos->width  = rpos.width * imgs->page_width;
  pos->height = rpos.height * imgs->page_height;
}

void dt_printing_setup_image(dt_images_box *imgs, const int idx,
                             const dt_imgid_t imgid, const int32_t width, const int32_t height,
                             const dt_alignment_t alignment)
{
  dt_image_box *box = &imgs->box[idx];

  if(box->imgid != imgid)
    dt_image_get_final_size(imgid, &box->img_width, &box->img_height);

  box->imgid      = imgid;
  box->exp_width  = width;
  box->exp_height = height;
  box->alignment  = alignment;

  // for the print (pdf) the origin is bottom/left, so y must be inverted compared to
  // screen coordinates.
  // FIXME: this is awkward that this is set then set again
  box->print.x      = box->pos.x * imgs->page_width;
  box->print.y      = box->pos.y * imgs->page_height;
  box->print.width  = box->pos.width * imgs->page_width;
  box->print.height = box->pos.height * imgs->page_height;

  dt_image_pos pos;
  _align_pos(&box->print, box->alignment, box->exp_width, box->exp_height, &pos);

  box->print.x = pos.x;
  box->print.y = imgs->page_height - (pos.y + pos.height);
  box->print.width = pos.width;
  box->print.height = pos.height;

  // compute image size on display

  box->dis_width  = box->img_width;
  box->dis_height = box->img_height;

  if(box->dis_width > box->screen.width)
  {
    const float scale =  box->screen.width / (float)box->dis_width;
    box->dis_width = box->screen.width;
    box->dis_height = (int32_t)(((float)box->dis_height + 0.5f) * scale);
  }

  if(box->dis_height > box->screen.height)
  {
    const float scale = box->screen.height / (float)box->dis_height;
    box->dis_height = box->screen.height;
    box->dis_width = (int32_t)(((float)box->dis_width + 0.5f) * scale);
  }
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on

