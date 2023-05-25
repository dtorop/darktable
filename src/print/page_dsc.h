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

#pragma once

#include <glib-object.h>

#include "common/cups_print.h"

G_BEGIN_DECLS

typedef enum page_orientation_t
{
  ORIENT_PORTRAIT = 0,
  ORIENT_LANDSCAPE,
  ORIENT_N // needs to be the last one
} page_orientation_t;

#define DT_TYPE_PAGE_DSC (dt_page_dsc_get_type())
G_DECLARE_FINAL_TYPE(DTPageDsc, dt_page_dsc, DT, PAGE_DSC, GObject)

DTPageDsc *dt_page_dsc_new();
void dt_page_dsc_set_paper(DTPageDsc *self, const gchar *const paper_name,
                              const dt_paper_info_t *const paper);

// FIXME: make getter/setter functions for commonly used attributes and to mass-set margins and other box-ish attributes (unless we use a struct to contain these?)
//gboolean dt_page_dsc_get_value(DTPageDsc *d, double *value);
//void dt_page_dsc_set_value(DTPageDsc *d, double value);

G_END_DECLS

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
