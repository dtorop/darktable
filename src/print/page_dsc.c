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

#include "common/cups_print.h"
#include "common/darktable.h"
#include "control/conf.h"
#include "print/page_dsc.h"


// FIXME: if GTK >= 2.74 could use G_DEFINE_ENUM_TYPE and G_DEFINE_ENUM_VALUE
// FIXME: move this to header file, share with print_settings?
typedef enum _unit_t
{
  UNIT_MM = 0,
  UNIT_CM,
  UNIT_IN,
  UNIT_N // needs to be the last one
} _unit_t;

static const float units[UNIT_N] = { 1.0f, 0.1f, 1.0f/25.4f };
static const gchar *_unit_names[] = { N_("mm"), N_("cm"), N_("inch"), NULL };


// FIXME: keep the hierarchical definition in cups_print.h or flat style?
// FIXME: retain the dt_print_info_t and dt_images_box types from extant print code?

struct _DTPageDsc {
  GObject parent;
  _unit_t unit;
  page_orientation_t orientation;
  // units: mm
  gdouble margin_top, margin_bottom, margin_left, margin_right;
  gboolean margin_lock;
  gdouble paper_width, paper_height;
#if 0
  // FIXME: this should go in a separte model which handles printer operations, including interacting with cups to set printer/paper
  gchar *printer_name;
  gboolean printer_is_turboprint;
#endif
  // FIXME: HW margins info should come from that separate printer model, rather than being stored/maintained here
  double hw_margin_top, hw_margin_bottom, hw_margin_left, hw_margin_right;
  guint resolution;
  // FIXME: we may not need all this data! but just name and/or common_name?
  dt_paper_info_t paper;
};

// FIXME: should be G_DEFINE_FINAL_TYPE?
G_DEFINE_TYPE(DTPageDsc, dt_page_dsc, G_TYPE_OBJECT)

enum {
  PROP_0,
  PROP_UNITS,
  PROP_ORIENTATION,
  // FIXME: s/PROP_MARGIN_/PROP_PAGE_MARGIN_/?
  PROP_MARGIN_TOP,
  PROP_MARGIN_BOTTOM,
  PROP_MARGIN_LEFT,
  PROP_MARGIN_RIGHT,
  PROP_MARGIN_LOCK,
  PROP_PAPER_WIDTH,
  PROP_PAPER_HEIGHT,
#if 0
  PROP_PRINTER_NAME,
  PROP_PRINTER_IS_TURBOPRINT,
#endif
  PROP_HW_MARGIN_TOP,
  PROP_HW_MARGIN_BOTTOM,
  PROP_HW_MARGIN_LEFT,
  PROP_HW_MARGIN_RIGHT,
  PROP_RESOLUTION,
  NUM_PROPERTIES
};

static GParamSpec *page_dsc_props[NUM_PROPERTIES] = { NULL, };


static void dt_page_dsc_set_property(GObject *object, guint prop_id, const GValue *value, GParamSpec *pspec) {
  DTPageDsc *self = DT_PAGE_DSC(object);

  const double from_mm = units[self->unit];
  //printf("from_mm is %f\n", from_mm);

  switch(prop_id)
  {
    case PROP_ORIENTATION:
      self->orientation = g_value_get_uint(value);
      break;
    case PROP_UNITS:
      self->unit = g_value_get_uint(value);
      dt_conf_set_string("plugins/print/print/unit", _unit_names[self->unit]);
      printf("dt_page_dsc_set_property: unit changed to %d '%s'\n", self->unit, _unit_names[self->unit]);
      break;
    case PROP_MARGIN_TOP:
      self->margin_top = g_value_get_double(value);
      printf("dt_page_dsc_set_property: setting top_margin to %f mm = %f %s\n", self->margin_top, self->margin_top * from_mm, _unit_names[self->unit]);
      dt_conf_set_float("plugins/print/print/top_margin",
                        self->margin_top * from_mm);
      break;
    case PROP_MARGIN_BOTTOM:
      self->margin_bottom = g_value_get_double(value);
      printf("dt_page_dsc_set_property: setting bottom_margin to %f mm = %f %s\n", self->margin_bottom, self->margin_bottom * from_mm, _unit_names[self->unit]);
      dt_conf_set_float("plugins/print/print/bottom_margin",
                        self->margin_bottom * from_mm);
      break;
    case PROP_MARGIN_LEFT:
      self->margin_left = g_value_get_double(value);
      dt_conf_set_float("plugins/print/print/left_margin",
                        self->margin_left * from_mm);
      break;
    case PROP_MARGIN_RIGHT:
      self->margin_right = g_value_get_double(value);
      dt_conf_set_float("plugins/print/print/right_margin",
                        self->margin_right * from_mm);
      break;
    case PROP_MARGIN_LOCK:
      self->margin_lock = g_value_get_boolean(value);
      dt_conf_set_float("plugins/print/print/lock_borders", self->margin_lock);
      break;
    case PROP_PAPER_WIDTH:
      self->paper_width = g_value_get_double(value);
      printf("set paper width to %f\n", self->paper_width);
      break;
    case PROP_PAPER_HEIGHT:
      self->paper_height = g_value_get_double(value);
      printf("set paper width to %f\n", self->paper_height);
      break;
#if 0
    case PROP_PRINTER_NAME:
      g_free(self->printer_name);
      self->printer_name = g_value_dup_string(g_value_get_string(value));
      dt_conf_set_string("plugins/print/print/printer", self->printer_name);
      break;
    case PROP_PRINTER_IS_TURBOPRINT:
      self->printer_is_turboprint = g_value_get_boolean(value);
#endif
    case PROP_HW_MARGIN_TOP:
      self->hw_margin_top = g_value_get_double(value);
      break;
    case PROP_HW_MARGIN_BOTTOM:
      self->hw_margin_bottom = g_value_get_double(value);
      break;
    case PROP_HW_MARGIN_LEFT:
      self->hw_margin_left = g_value_get_double(value);
      break;
    case PROP_HW_MARGIN_RIGHT:
      self->hw_margin_right = g_value_get_double(value);
      break;
    case PROP_RESOLUTION:
      self->resolution = g_value_get_uint(value);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID(object, prop_id, pspec);
      break;
  }
}

static void dt_page_dsc_get_property(GObject *object, guint prop_id, GValue *value, GParamSpec *pspec) {
  DTPageDsc *self = DT_PAGE_DSC(object);

  switch(prop_id)
  {
    case PROP_ORIENTATION:
      g_value_set_uint(value, self->orientation);
      break;
    case PROP_UNITS:
      g_value_set_uint(value, self->unit);
      break;
    case PROP_MARGIN_TOP:
      g_value_set_double(value, self->margin_top);
      break;
    case PROP_MARGIN_BOTTOM:
      g_value_set_double(value, self->margin_bottom);
      break;
    case PROP_MARGIN_LEFT:
      g_value_set_double(value, self->margin_left);
      break;
    case PROP_MARGIN_RIGHT:
      g_value_set_double(value, self->margin_right);
      break;
    case PROP_MARGIN_LOCK:
      g_value_set_boolean(value, self->margin_lock);
      break;
    case PROP_PAPER_WIDTH:
      g_value_set_double(value, self->paper_width);
      break;
    case PROP_PAPER_HEIGHT:
      g_value_set_double(value, self->paper_height);
      break;
#if 0
    case PROP_PRINTER_NAME:
      g_value_set_string(value, self->printer_name);
      break;
    case PROP_PRINTER_IS_TURBOPRINT:
      g_value_set_boolean(value, self->printer_is_turboprint);
      break;
#endif
    case PROP_HW_MARGIN_TOP:
      g_value_set_double(value, self->hw_margin_top);
      break;
    case PROP_HW_MARGIN_BOTTOM:
      g_value_set_double(value, self->hw_margin_bottom);
      break;
    case PROP_HW_MARGIN_LEFT:
      g_value_set_double(value, self->hw_margin_left);
      break;
    case PROP_HW_MARGIN_RIGHT:
      g_value_set_double(value, self->hw_margin_right);
      break;
    case PROP_RESOLUTION:
      g_value_set_uint(value, self->resolution);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID(object, prop_id, pspec);
      break;
  }
}

static void dt_page_dsc_class_init(DTPageDscClass *class)
{
  GObjectClass *gobject_class = G_OBJECT_CLASS(class);

  gobject_class->set_property = dt_page_dsc_set_property;
  gobject_class->get_property = dt_page_dsc_get_property;

  page_dsc_props[PROP_UNITS] =
    g_param_spec_uint("units", "units", "Measurement units",
                      UNIT_MM, UNIT_N-1, UNIT_MM, G_PARAM_READWRITE);
  // FIXME: should be an enum
  page_dsc_props[PROP_ORIENTATION] =
    g_param_spec_uint("orientation", "orient", "Page orientation",
                      ORIENT_PORTRAIT, ORIENT_N-1, ORIENT_PORTRAIT, G_PARAM_READWRITE);
  page_dsc_props[PROP_MARGIN_TOP] =
    g_param_spec_double("margin-top", "mtop", "Page top margin in mm",
                        0, 1000.0, 17.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_MARGIN_BOTTOM] =
    g_param_spec_double("margin-bottom", "mbottom", "Page bottom margin in mm",
                        0, 1000.0, 17.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_MARGIN_LEFT] =
    g_param_spec_double("margin-left", "mleft", "Page left margin in mm",
                        0, 1000.0, 17.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_MARGIN_RIGHT] =
    g_param_spec_double("margin-right", "mright", "Page right margin in mm",
                        0, 1000.0, 17.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_MARGIN_LOCK] =
    g_param_spec_boolean("margin-lock", "mlock", "Lock all margins to top margin value",
                         TRUE, G_PARAM_READWRITE);
  page_dsc_props[PROP_PAPER_WIDTH] =
    g_param_spec_double("paper-width", "pwidth", "Page width in mm",
                        1, 10000.0, 210.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_PAPER_HEIGHT] =
    g_param_spec_double("paper-height", "pheight", "Page height in mm",
                        1, 10000.0, 297.0, G_PARAM_READWRITE);
#if 0
  page_dsc_props[PROP_PRINTER_NAME] =
    g_param_spec_string("printer-name", "pname", "Printer name",
                        NULL, G_PARAM_READWRITE);
  page_dsc_props[PROP_PRINTER_IS_TURBOPRINT] =
    g_param_spec_boolean("printer-is-turboprint", "turboprint",
                         FALSE, G_PARAM_READWRITE);
#endif
  page_dsc_props[PROP_HW_MARGIN_TOP] =
    g_param_spec_double("hw-margin-top", "hwtop",
                        "Printer hardware top margin in mm",
                        0, 1000.0, 0.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_HW_MARGIN_BOTTOM] =
    g_param_spec_double("hw-margin-bottom", "hwbottom",
                        "Printer hardware bottom margin in mm",
                        0, 1000.0, 0.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_HW_MARGIN_LEFT] =
    g_param_spec_double("hw-margin-left", "hwleft",
                        "Printer hardware left margin in mm",
                        0, 1000.0, 0.0, G_PARAM_READWRITE);
  page_dsc_props[PROP_HW_MARGIN_RIGHT] =
    g_param_spec_double("hw-margin-right", "hwright",
                        "Printer hardware right margin in mm",
                        0, 1000.0, 0.0, G_PARAM_READWRITE);
  // FIXME: do we max out resolution at 360?
  page_dsc_props[PROP_RESOLUTION] =
    g_param_spec_uint("resolution", "res", "Printer resolution",
                      1, 2880, 360, G_PARAM_READWRITE);

  g_object_class_install_properties(gobject_class, NUM_PROPERTIES, page_dsc_props);  
}

static void dt_page_dsc_init(DTPageDsc *self)
{
  const char *str = dt_conf_get_string_const("plugins/print/print/unit");
  const char **names = _unit_names;
  for(_unit_t i=0; *names; names++, i++)
    if(g_strcmp0(str, *names) == 0)
      self->unit = i;

  const double to_mm = 1.0 / units[self->unit];
  // FIXME: wouldn't it be better to save these in conf as a known unit (mm), and convert them as needed?
  self->margin_top = dt_conf_get_float("plugins/print/print/top_margin") * to_mm;
  self->margin_bottom = dt_conf_get_float("plugins/print/print/bottom_margin") * to_mm;
  self->margin_left =
    dt_conf_get_float("plugins/print/print/left_margin") * to_mm;
  self->margin_right = dt_conf_get_float("plugins/print/print/right_margin") * to_mm;
#if 0
  // we don't load this from conf, as that printer may not be
  // connected -- we don't know until the printer detect stage
  self->printer_name = NULL;
#endif
  // FIXME: should load paper name from conf and when instantiate print settings set paper name dropdown based on that which will load in paper description (once CUPS printers is loaded)
}

DTPageDsc *dt_page_dsc_new(void) {
  DTPageDsc *self;

  self = g_object_new(DT_TYPE_PAGE_DSC, NULL);

  return self;
}

#if 0
void dt_page_dsc_finalize(DTPageDsc *self)
{
  g_free(self->printer_name);
  G_OBJECT_CLASS(dt_page_dsc_parent_class)->finalize(self);
}
#endif

void dt_page_dsc_set_paper(DTPageDsc *self, const gchar *const paper_name,
                              const dt_paper_info_t *const paper)
{
  // FIXME: do we need the whole dt_paper_info_t structure? or should just g_object_set width/height in print_settings?
  memcpy(&self->paper, paper, sizeof(dt_paper_info_t));
  // FIXME: the print settings code sets this to the paper name from combobox, which could be name or common_name -- should we replicate that? or just set to paper->name?
  printf("dt_page_dsc_set_paper: to '%s' %fx%f\n", paper_name, paper->width, paper->height);
  dt_conf_set_string("plugins/print/print/paper", paper_name);
  // FIXME: some papers in "portrait" mode are actually wider than tall -- in that case should we manually flip them to landscape, then reflip to portrait when running print job
  // FIXME: use g_object_set or set struct directly and raise a notify signal?
  g_object_set(self,
               "paper-width", paper->width,
               "paper-height", paper->height,
               NULL);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
