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
/** a class to draw a page, margins, and printable area, and hold the layout body */

#pragma once

#include "print/page_dsc.h"


typedef struct _dt_layout_page_t dt_layout_page_t;

// FIXME: once set it up via gobject, don't have to pass in printer
dt_layout_page_t *dt_layout_page_new(DTPageDsc *page_dsc);
void dt_layout_page_destroy(dt_layout_page_t *const page);

GtkWidget *dt_layout_page_get_widget(dt_layout_page_t *page);

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
