/*
    This file is part of darktable,
    Copyright (C) 2023 darktable developers.

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
/** a class to draw the area of a page within its margins, and hold layout boxes */

#pragma once

#include "print/page_dsc.h"

typedef struct _dt_layout_body_t dt_layout_body_t;

dt_layout_body_t *dt_layout_body_new(DTPageDsc *page_dsc);
void dt_layout_body_destroy(dt_layout_body_t *const body);

GtkWidget *dt_layout_body_get_widget(dt_layout_body_t *body);

#if 0
void dt_layout_body_add_box(dt_layout_body_t *const body,
                            const dt_image_box *const box,
                            const dt_imgid_t imgid);
#endif


// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
