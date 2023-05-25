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

#include "dtgtk/layout_page.h"

#include "control/control.h"
#include "dtgtk/layout_body.h"
#include "gui/drag_and_drop.h"
#include "gui/gtk.h"

// FIXME: make the print code in print_settings work with the gobject data
// FIXME: make the get_params/set_params in print_settings work with the gobject data
// FIXME: should use bindings to keep paper background & layout box areas in sync (aspect ratio + margins)?
// FIXME: add grid overlay, and make snap to grid work
// FIXME: for debugging make select area have an X inside ito
// FIXME: make a drag-drop handler on print_layout->w_main which dims the images in each layout box?
// FIXME: split the right panel into print settings, image layout (x/y/width/height, image dimensions & ppi, alignment), and a print button with no expander at the bottom

struct _dt_layout_page_t
{
  // FIXME: reference these data structures, or just hold locally the relevant #'s?
  //dt_printer_info_t printer;
  //dt_paper_info_t paper;
  //dt_page_setup_t page_setup;

  DTPageDsc *page_dsc;

  // FIXME: if this data is encapsulated in this class/widget, then we don't have to keep it in other data structures
  // FIXME: units should be in printer pixels (determined from resolution)?
  // FIXME: if units are in mm, then do we need printer_resolution here?
  // from dt_printer_info_t
  //int printer_resolution;
  double hw_margin_top, hw_margin_bottom, hw_margin_left, hw_margin_right;
  // from dt_paper_info_t, in mm
  gdouble paper_width, paper_height;
  // from dt_page_setup_t, in mm
  gdouble margin_top, margin_bottom, margin_left, margin_right;

  // FIXME: there's a bunch of other information from cups_print.h we could encapsulate -- the current medium, the current profile and rendering intent for the paper, whether to use turboprint -- should we store this and take it off the hands of external data strucutres? or have gobject hold it?

  // page is a sandwich of (bottom-to-top):
  //   1) background of center area
  //   2) aspectframe for page which contains via an overlay
  //      a) paper margins
  //      b) paper area inside margins
  //   3) hardware margins
  //   4) widgets to draw body w/margins (matches body in aspectframe above)
  // FIXME: how many of these widgets need access via struct outside of init? can some be local?
  GtkWidget *w_main;         // GtkOverlay -- contains all other widgets
  GtkWidget *w_background;   // GtkDrawingArea -- center view background
  GtkWidget *w_aspect1;      // GtkAspectFrame -- contains margins/content background
  GtkWidget *w_margins;      // GtkDrawingArea -- paper outside margins
  GtkWidget *w_content;      // GtkDrawingArea -- paper inside margins;
  GtkWidget *w_hw_margins;   // GtkDrawingArea -- tick-marks showing hardware margins of printer
  GtkWidget *w_aspect2;      // GtkAspectFrame -- contains layout boxes
  // FIXME: an early version of print module allowed for negative margins to expand a photo, do we want to bring that back?

  dt_layout_body_t *body;
};

static void _recalc_body_margins(dt_layout_page_t *page)
{
  const int body_width_px = gtk_widget_get_allocated_width(page->w_margins);
  const int body_height_px = gtk_widget_get_allocated_height(page->w_margins);

  // FIXME: this boilerplate is repeated, make getter method for orientation or is_landscape
  page_orientation_t orientation;
  g_object_get(page->page_dsc, "orientation", &orientation, NULL);
  const gboolean landscape = (orientation == ORIENT_LANDSCAPE);

  const double pg_width = landscape ? page->paper_height : page->paper_width;
  const double pg_height = landscape ? page->paper_width : page->paper_height;
  printf("_recalc_body_margins called, w_page %p, alloc %d x %d, landscape %d page %f x %f\n", page->w_margins, body_width_px, body_height_px, landscape, pg_width, pg_height);

  // page margins are set by the user from the GUI, and the top margin
  // is *always* at the top of the screen (doesn't change with
  // page orientation)
  const double margin_top_px = page->margin_top / pg_height * body_height_px;
  const double margin_bottom_px = page->margin_bottom / pg_height * body_height_px;
  const double margin_left_px = page->margin_left / pg_width * body_width_px;
  const double margin_right_px = page->margin_right / pg_width * body_width_px;
  printf("  margins %f %f %f %f mm, %f %f %f %f px\n", page->margin_top, page->margin_bottom, page->margin_left, page->margin_right, margin_top_px, margin_bottom_px, margin_left_px, margin_right_px);

  gtk_widget_set_margin_top(page->w_content, round(margin_top_px));
  gtk_widget_set_margin_bottom(page->w_content, round(margin_bottom_px));
  gtk_widget_set_margin_start(page->w_content, round(margin_left_px));
  gtk_widget_set_margin_end(page->w_content, round(margin_right_px));

  // FIXME: do need this? if so do it via bind?
  GtkWidget *w_body = dt_layout_body_get_widget(page->body);
  gtk_widget_set_margin_top(w_body, round(margin_top_px));
  gtk_widget_set_margin_bottom(w_body, round(margin_bottom_px));
  gtk_widget_set_margin_start(w_body, round(margin_left_px));
  gtk_widget_set_margin_end(w_body, round(margin_right_px));
}

static gboolean _event_configure_margins(GtkWidget *widget,
                                         const GdkEventConfigure *const event,
                                         dt_layout_page_t *page)
{
  printf("_event_configure_margins called\n");
  _recalc_body_margins(page);
  return FALSE;
}

static gboolean _event_draw_rect(GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
  // Background-color is CSS dependent, hence same draw code can
  // handle center area background, paper margins and content. The
  // page content widget is simply drawn over the margins widget.
  GtkStyleContext *context = gtk_widget_get_style_context(widget);
  guint width = gtk_widget_get_allocated_width (widget);
  guint height = gtk_widget_get_allocated_height (widget);
  gtk_render_background(context, cr, 0, 0, width, height);
  gtk_render_frame(context, cr, 0, 0, width, height);
  return FALSE;
}

// FIXME: instead draw these over background, margins, and page widgets so don't need to deal with a separate overlay
static gboolean _event_draw_hw_margins(GtkWidget *widget, cairo_t *cr,
                                       dt_layout_page_t *page)
{
  gint page_x, page_y;
  if(!gtk_widget_translate_coordinates(page->w_margins, widget, 0, 0, &page_x, &page_y))
  {
    // FIXME: presumably the page widget is not yet realized
    dt_print(DT_DEBUG_ALWAYS, "_event_draw_hw_margins: could not translate from page to hw coordinates\n");
    return FALSE;
  }

  GtkStyleContext *context = gtk_widget_get_style_context(widget);
  // it would be more succinct to use gtk_render_line(), but it always
  // draws 1 px wide lines, and the tick-marks should be 2 px
  GdkRGBA color;
  gtk_style_context_get_color(context,
                              gtk_style_context_get_state (context),
                              &color);
  gdk_cairo_set_source_rgba(cr, &color);

  //printf("_event_draw_hw_margins: page coords %dx%d color %f,%f,%f,%f\n", page_x, page_y, color.red, color.green, color.blue, color.alpha);
  cairo_translate(cr, page_x, page_y);

  // display non-printable area
  page_orientation_t orientation;
  g_object_get(page->page_dsc, "orientation", &orientation, NULL);
  const gboolean landscape = (orientation == ORIENT_LANDSCAPE);
  // FIXME: better to handle landscape mode by rotating the cairo drawing area by 90 degrees
  const double pg_width_mm = landscape ? page->paper_height : page->paper_width;
  const double pg_height_mm = landscape ? page->paper_width : page->paper_height;

  // all measurements in mm
  const double np_top = landscape ? page->hw_margin_right : page->hw_margin_top;
  const double np_left = landscape ? page->hw_margin_top : page->hw_margin_left;
  const double np_right = landscape ? page->hw_margin_bottom : page->hw_margin_right;
  const double np_bottom = landscape ? page->hw_margin_left : page->hw_margin_bottom;

  const int paper_width_px = gtk_widget_get_allocated_width(page->w_margins);
  const int paper_height_px = gtk_widget_get_allocated_height(page->w_margins);

  const double np1x = (np_left / pg_width_mm) * paper_width_px;
  const double np1y = (np_top / pg_height_mm) * paper_height_px;
  const double np2x = ((pg_width_mm - np_right) / pg_width_mm) * paper_width_px;
  const double np2y = ((pg_height_mm - np_bottom) / pg_height_mm) * paper_height_px;
  const double tick_width = DT_PIXEL_APPLY_DPI(10.0);

  // top-left
  cairo_move_to(cr, np1x-tick_width, np1y);
  cairo_line_to(cr, np1x, np1y); cairo_line_to(cr, np1x, np1y-tick_width);
  cairo_stroke(cr);

  // top-right
  cairo_move_to(cr, np2x+tick_width, np1y);
  cairo_line_to(cr, np2x, np1y); cairo_line_to(cr, np2x, np1y-tick_width);
  cairo_stroke(cr);

  // bottom-left
  cairo_move_to(cr, np1x-tick_width, np2y);
  cairo_line_to(cr, np1x, np2y); cairo_line_to (cr, np1x, np2y+tick_width);
  cairo_stroke(cr);

  // bottom-right
  cairo_move_to(cr, np2x+tick_width, np2y);
  cairo_line_to(cr, np2x, np2y); cairo_line_to (cr, np2x, np2y+tick_width);
  cairo_stroke(cr);

  return FALSE;
}

void _update_aspect(dt_layout_page_t *const page)
{
  gdouble aspect;
  if(page->paper_width == 0.0 || page->paper_height == 0.0)
  {
    // the paper somehow not initialized
    // FIXME: should probably de-activate print button
    gtk_widget_hide(page->w_aspect1);
    gtk_widget_hide(page->w_aspect2);
    dt_control_log(_("please select a printer and paper"));
  }
  else
  {
    page_orientation_t orientation;
    g_object_get(page->page_dsc, "orientation", &orientation, NULL);
    aspect = (orientation == ORIENT_LANDSCAPE) ?
      page->paper_height / page->paper_width :
      page->paper_width / page->paper_height;

    // two matching aspect frames, one for background of margins/page,
    // one for layout boxes
    // FIXME: could also bind the ratio parameter?
    gtk_aspect_frame_set(GTK_ASPECT_FRAME(page->w_aspect1), 0.5f, 0.5f, aspect, FALSE);
    gtk_aspect_frame_set(GTK_ASPECT_FRAME(page->w_aspect2), 0.5f, 0.5f, aspect, FALSE);
    gtk_widget_show(page->w_aspect1);
    gtk_widget_show(page->w_aspect2);
  }
}

// FIXME: _orientation_changed and paper width/height changed all are similar, can combine just to an _update_aspect() call?
void _orientation_changed(GObject *self, GParamSpec *pspec,
                          dt_layout_page_t *const page)
{
  _update_aspect(page);
}

void _paper_width_changed(GObject *self, GParamSpec *pspec,
                          dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->paper_width, NULL);
  printf("paper width changed to %f\n", page->paper_width);
  _update_aspect(page);
}

void _paper_height_changed(GObject *self, GParamSpec *pspec,
                          dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->paper_height, NULL);
  printf("paper height changed to %f\n", page->paper_height);
  _update_aspect(page);
}

// FIXME: can make these one function? should have a margin type which consists of a struct of 4 doubles to compress these updates and unify these functions?
// FIXME: can use bind and a conversion function for this?
void _top_margin_changed(GObject *self, GParamSpec *pspec,
                         dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->margin_top, NULL);
  printf("_top_margin_changed %s -> %f\n", pspec->name, page->margin_top);
  _recalc_body_margins(page);
}

void _bottom_margin_changed(GObject *self, GParamSpec *pspec,
                            dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->margin_bottom, NULL);
  printf("_bottom_margin_changed %s -> %f\n", pspec->name, page->margin_bottom);
  _recalc_body_margins(page);
}

void _left_margin_changed(GObject *self, GParamSpec *pspec,
                          dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->margin_left, NULL);
  printf("_left_margin_changed %s -> %f\n", pspec->name, page->margin_left);
  _recalc_body_margins(page);
}

void _right_margin_changed(GObject *self, GParamSpec *pspec,
                           dt_layout_page_t *const page)
{
  g_object_get(self, pspec->name, &page->margin_right, NULL);
  printf("_right_margin_changed %s -> %f\n", pspec->name, page->margin_right);
  _recalc_body_margins(page);
}

void _hw_margin_changed(GObject *self, GParamSpec *pspec,
                        dt_layout_page_t *const page)
{
  // FIXME: instead should get these from gobject in _event_draw_hw_margins()?
  // FIXME: can we store a struct in GObject with all four dimensions so this isn't called x4?
  g_object_get(page->page_dsc,
               "hw-margin-top", &page->hw_margin_top,
               "hw-margin-bottom", &page->hw_margin_bottom,
               "hw-margin-left", &page->hw_margin_left,
               "hw-margin-right", &page->hw_margin_right,
               NULL);
  gtk_widget_queue_draw(page->w_hw_margins);
}

GtkWidget *dt_layout_page_get_widget(dt_layout_page_t *page)
{
  return page->w_main;
}

dt_layout_page_t *dt_layout_page_new(DTPageDsc *page_dsc)
{
  dt_layout_page_t *page =
    (dt_layout_page_t *)calloc(1, sizeof(dt_layout_page_t));

  page->page_dsc = g_object_ref(page_dsc);

  // FIXME: should get these here and store them internally, or just get them as needed from gobject?
  // FIXME: margins come in via mm, put them in screen pixels or percentages of page?
  g_object_get(page->page_dsc,
               "margin-top", &page->margin_top,
               "margin-bottom", &page->margin_bottom,
               "margin-left", &page->margin_left,
               "margin-right", &page->margin_right,
               "paper-width", &page->paper_width,
               "paper-height", &page->paper_height,
               "hw-margin-top", &page->hw_margin_top,
               "hw-margin-bottom", &page->hw_margin_bottom,
               "hw-margin-left", &page->hw_margin_left,
               "hw-margin-right", &page->hw_margin_right,
               NULL);
  printf("new page %p, paper %fx%f hw margins %f,%f,%f,%f\n", page, page->paper_width, page->paper_height, page->hw_margin_top, page->hw_margin_bottom, page->hw_margin_left, page->hw_margin_right);

  // view background outside of page
  page->w_background = gtk_drawing_area_new();
  gtk_widget_set_name(page->w_background, "print-view-bkgd");

  // page margins
  page->w_margins = gtk_drawing_area_new();
  gtk_widget_set_name(page->w_margins, "print-paper-margins");

  // inner content (Note that this widget draws the background of the
  // paper containing the content, but page->body, above it, actually
  // contains the layout boxes. This is so that the hardware margins
  // can be above the paper and not block events to layout boxes.
  page->w_content = gtk_drawing_area_new();
  gtk_widget_set_name(page->w_content, "print-paper-content");

  GtkWidget *w_overlay = gtk_overlay_new();
  gtk_widget_set_name(w_overlay, "print-paper-overlay");
  gtk_container_add(GTK_CONTAINER(w_overlay), page->w_margins);
  gtk_overlay_add_overlay(GTK_OVERLAY(w_overlay), page->w_content);

  page->w_aspect1 = gtk_aspect_frame_new(NULL, 0.5f, 0.5f, 1.0f, FALSE);
  gtk_widget_set_name(page->w_aspect1, "print-paper-frame");
  gtk_container_add(GTK_CONTAINER(page->w_aspect1), w_overlay);
  gtk_widget_set_size_request(page->w_aspect1, DT_PIXEL_APPLY_DPI(200),
                              DT_PIXEL_APPLY_DPI(200));

  // FIXME: for call-outs on drawing do we use a GtkLabel?
  // TODO: if there are call-out overlays, render ends via gtk_render_arrow() and gtk_render_line()

  // tick marks to show hardware margins -- this may extend beyond paper widget
  page->w_hw_margins = gtk_drawing_area_new();
  gtk_widget_set_name(page->w_hw_margins, "print-hw-margins");

  // body of page (area w/in margins)
  page->body = dt_layout_body_new(page_dsc);
  GtkWidget *w_body = dt_layout_body_get_widget(page->body);
  // FIXME: set direction in dt_layout_body_new()?
  gtk_widget_set_direction(w_body, GTK_TEXT_DIR_LTR);

  // container for page layout boxes
  // FIXME: do need to bind margins to the content widget above
  page->w_aspect2 = gtk_aspect_frame_new(NULL, 0.5f, 0.5f, 1.0f, FALSE);
  gtk_widget_set_name(page->w_aspect2, "print-layout-frame");
  gtk_container_add(GTK_CONTAINER(page->w_aspect2), w_body);

  // FIXME: needed or is/should be called on configure event?
  _update_aspect(page);

  page->w_main = gtk_overlay_new();
  gtk_widget_set_name(page->w_main, "print-layout-main");
  gtk_container_add(GTK_CONTAINER(page->w_main), page->w_background);
  gtk_overlay_add_overlay(GTK_OVERLAY(page->w_main), page->w_aspect1);
  gtk_overlay_add_overlay(GTK_OVERLAY(page->w_main), page->w_hw_margins);
  gtk_overlay_add_overlay(GTK_OVERLAY(page->w_main), page->w_aspect2);

  // FIXME: we should dt_gui_add_help_link()

  // FIXME: either of these useful?
  //gtk_widget_set_app_paintable(layout->widget, TRUE);
  //gtk_widget_set_can_focus(layout->widget, TRUE);

  // FIXME: if use CSS for center area background, can use _event_draw_paper code here too, will just be _event_draw_bkgd_rect or such for all three
  g_signal_connect(G_OBJECT(page->w_background), "draw",
                   G_CALLBACK(_event_draw_rect), NULL);
  g_signal_connect(G_OBJECT(page->w_margins), "draw",
                   G_CALLBACK(_event_draw_rect), NULL);
  g_signal_connect(G_OBJECT(page->w_content), "draw",
                   G_CALLBACK(_event_draw_rect), NULL);
  g_signal_connect(G_OBJECT(page->w_margins), "configure-event",
                   G_CALLBACK(_event_configure_margins), page);
  g_signal_connect(G_OBJECT(page->w_hw_margins), "draw",
                   G_CALLBACK(_event_draw_hw_margins), page);

  g_signal_connect(G_OBJECT(page->page_dsc), "notify::orientation",
                   G_CALLBACK(_orientation_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::paper-width",
                   G_CALLBACK(_paper_width_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::paper-height",
                   G_CALLBACK(_paper_height_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::margin-top",
                   G_CALLBACK(_top_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::margin-bottom",
                   G_CALLBACK(_bottom_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::margin-left",
                   G_CALLBACK(_left_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::margin-right",
                   G_CALLBACK(_right_margin_changed), page);
  // FIXME: as these all probably change at once, instead notify on printer name changed?
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::hw-margin-top",
                   G_CALLBACK(_hw_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::hw-margin-bottom",
                   G_CALLBACK(_hw_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::hw-margin-left",
                   G_CALLBACK(_hw_margin_changed), page);
  g_signal_connect(G_OBJECT(page->page_dsc), "notify::hw-margin-right",
                   G_CALLBACK(_hw_margin_changed), page);

  gtk_widget_show(page->w_background);
  gtk_widget_show(page->w_margins);
  gtk_widget_show(page->w_content);
  gtk_widget_show(w_overlay);
  gtk_widget_show(page->w_aspect1);
  gtk_widget_show(page->w_hw_margins);
  gtk_widget_show(page->w_aspect2);
  gtk_widget_show(page->w_main);

#if 0
  gtk_widget_add_events(layout->w_fixed, GDK_POINTER_MOTION_MASK);
  g_signal_connect(G_OBJECT(layout->w_fixed), "motion-notify-event", G_CALLBACK(_mouse_moved_fixed), layout);
#endif

  return page;
}

void dt_layout_page_destroy(dt_layout_page_t *const page)
{
  // page_dsc is maintained when enter/leave print view, hence we
  // must disconnect these handlers before we destroy page
  g_signal_handlers_disconnect_by_data(G_OBJECT(page->page_dsc), page);

  dt_layout_body_destroy(page->body);

  // this will remove the main widget from its container
  // FIXME: ???
  if(page->w_aspect1) gtk_widget_destroy(page->w_main);

  g_object_unref(page->page_dsc);

  free(page);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
