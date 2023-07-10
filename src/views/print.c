/*
    This file is part of darktable,
    Copyright (C) 2014-2023 darktable developers.

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

/** this is the view for the print module.  */
#include "common/collection.h"
#include "common/cups_print.h"
#include "common/printing.h"
#include "common/darktable.h"
#include "common/debug.h"
#include "common/image_cache.h"
#include "common/selection.h"
#include "control/conf.h"
#include "control/control.h"
#include "develop/develop.h"
#include "dtgtk/thumbtable.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include "views/view.h"
#include "views/view_api.h"

#include <gdk/gdkkeysyms.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

DT_MODULE(1)

typedef struct dt_print_t
{
  dt_print_info_t *pinfo;
  dt_images_box *imgs;

  GtkWidget *w_main;         // GtkOverlay -- contains all other widgets
  GtkWidget *w_aspect1;      // GtkAspectFrame -- contains margins/content background
  GtkWidget *w_margins;      // GtkDrawingArea -- paper outside margins
  GtkWidget *w_content;      // GtkDrawingArea -- paper inside margins;
  GtkWidget *w_hw_margins;   // GtkDrawingArea -- tick-marks showing hardware margins of printer
  GtkWidget *w_page;         // GtkDrawingArea -- temp catchall
}
dt_print_t;

const char *name(const dt_view_t *self)
{
  return C_("view", "print");
}

uint32_t view(const dt_view_t *self)
{
  return DT_VIEW_PRINT;
}

static void _view_print_filmstrip_activate_callback(gpointer instance,
                                                    const dt_imgid_t imgid,
                                                    const dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t *)self->data;

  if(!dt_is_valid_imgid(imgid))
    return;

  // if the previous shown image is selected and the selection is unique
  // then we change the selected image to the new one
  if(dt_is_valid_imgid(prt->imgs->imgid_to_load))
  {
    sqlite3_stmt *stmt;
    DT_DEBUG_SQLITE3_PREPARE_V2(dt_database_get(darktable.db),
                                "SELECT m.imgid"
                                " FROM memory.collected_images as m, main.selected_images as s"
                                " WHERE m.imgid=s.imgid",
                                -1, &stmt, NULL);
    gboolean follow = FALSE;
    if(sqlite3_step(stmt) == SQLITE_ROW)
    {
      if(sqlite3_column_int(stmt, 0) == prt->imgs->imgid_to_load
         && sqlite3_step(stmt) != SQLITE_ROW)
      {
        follow = TRUE;
      }
    }
    sqlite3_finalize(stmt);
    if(follow)
    {
      dt_selection_select_single(darktable.selection, imgid);
    }
  }

  prt->imgs->imgid_to_load = imgid;

  // update the active images list
  dt_view_active_images_reset(FALSE);
  // this triggers print_settings callback which can load the activated image
  dt_view_active_images_add(imgid, TRUE);
}

static void _update_display_coords(dt_print_t *prt, int view_width, int view_height)
{
  // FIXME: if use aspect frame we can skip a lot of this work
  float px=.0f, py=.0f, pwidth=.0f, pheight=.0f;
  float ax=.0f, ay=.0f, awidth=.0f, aheight=.0f;
  gboolean borderless = FALSE;

  dt_get_print_layout(prt->pinfo, view_width, view_height,
                      &px, &py, &pwidth, &pheight,
                      &ax, &ay, &awidth, &aheight, &borderless);

  // record the screen page dimension. this will be used to draw the
  // page and to compute the actual layout of the areas placed over
  // the page.
  dt_printing_setup_display(prt->imgs,
                            px, py, pwidth, pheight,
                            ax, ay, awidth, aheight,
                            borderless);

  // FIXME: if some of the pixel dimensions aren't used outside of here, don't store theme! if we can put them in a form to make them easier here, do that
  if(prt->w_content)
  {
    gtk_widget_set_margin_start(prt->w_content, round(ax - px));
    gtk_widget_set_margin_end(prt->w_content, round((pwidth - awidth) - (ax - px)));
    gtk_widget_set_margin_top(prt->w_content, round(ay - py));
    gtk_widget_set_margin_bottom(prt->w_content, round((pheight - aheight) - (ay - py)));
  }

  if(prt->w_aspect1)
  {
    const dt_paper_info_t *paper = &prt->pinfo->paper;
    if(paper->width > 0.0 && paper->height > 0.0)
    {
      gdouble aspect = (prt->pinfo->page.landscape) ?
        paper->height / paper->width : paper->width / paper->height;

      gtk_aspect_frame_set(GTK_ASPECT_FRAME(prt->w_aspect1), 0.5f, 0.5f, aspect, FALSE);
    }
  }
}

static void _view_print_settings(const dt_view_t *view,
                                 dt_print_info_t *pinfo,
                                 dt_images_box *imgs)
{
  dt_print_t *prt = (dt_print_t *)view->data;

  prt->pinfo = pinfo;
  prt->imgs = imgs;

  _update_display_coords(prt, view->width, view->height);
  dt_control_queue_redraw();
}

void
init(dt_view_t *self)
{
  self->data = calloc(1, sizeof(dt_print_t));

  /* initialize CB to get the print settings from corresponding lib module */
  darktable.view_manager->proxy.print.view = self;
  darktable.view_manager->proxy.print.print_settings = _view_print_settings;
}

void cleanup(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t *)self->data;
  g_object_unref(prt->w_main);
  free(prt);
}

void configure(dt_view_t *self, int width, int height)
{
  dt_print_t *prt = (dt_print_t *)self->data;
  if(prt)
    _update_display_coords(prt, width, height);
}

void expose(dt_view_t *self, cairo_t *cr, int32_t width, int32_t height,
            int32_t pointerx, int32_t pointery)
{
  // FIXME: somehow is necessary?
}

static gboolean _event_draw_bkgd(GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
  dt_gui_gtk_set_source_rgb(cr, DT_GUI_COLOR_PRINT_BG);
  cairo_paint(cr);
  return FALSE;
}

static gboolean _event_draw_rect(GtkWidget *widget, cairo_t *cr, gpointer user_data)
{
  // Background-color is CSS dependent, hence same draw code can paper
  // margins and content. The page content widget is simply drawn over
  // the margins widget.
  GtkStyleContext *context = gtk_widget_get_style_context(widget);
  guint width = gtk_widget_get_allocated_width(widget);
  guint height = gtk_widget_get_allocated_height(widget);
  gtk_render_background(context, cr, 0, 0, width, height);
  gtk_render_frame(context, cr, 0, 0, width, height);
  return FALSE;
}

static gboolean _event_draw_hw_margins(GtkWidget *widget, cairo_t *cr,
                                       dt_print_t *prt)
{
  gint page_x, page_y;
  if(!gtk_widget_translate_coordinates(prt->w_margins, widget, 0, 0, &page_x, &page_y))
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
                              gtk_style_context_get_state(context),
                              &color);
  gdk_cairo_set_source_rgba(cr, &color);

  cairo_translate(cr, page_x, page_y);

  // page w/h in px
  const float pwidth = prt->imgs->screen.page.width;
  const float pheight = prt->imgs->screen.page.height;

  // page w/h in mm
  const float pg_width  = prt->pinfo->paper.width;
  const float pg_height = prt->pinfo->paper.height;

  // display non-printable area
  const dt_printer_info_t *prntr = &prt->pinfo->printer;
  const gboolean lndscp = prt->pinfo->page.landscape;
  const float np_top = lndscp ? prntr->hw_margin_right : prntr->hw_margin_top;
  const float np_left = lndscp ? prntr->hw_margin_top : prntr->hw_margin_left;
  const float np_right = lndscp ? prntr->hw_margin_bottom : prntr->hw_margin_right;
  const float np_bottom = lndscp ? prntr->hw_margin_left : prntr->hw_margin_bottom;

  const float np1x = (np_left / pg_width) * pwidth;
  const float np1y = (np_top / pg_height) * pheight;
  const float np2x = pwidth * (1.0f - np_right / pg_width);
  const float np2y = pheight * (1.0f - np_bottom / pg_height);

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

static gboolean _event_draw_page(GtkWidget *self, cairo_t *cr, dt_print_t *prt)
{
  // print page & borders only. Images are displayed in
  // gui_post_expose in print_settings module.
  if(!prt->pinfo)
    return FALSE;

  const float px = prt->imgs->screen.page.x;
  const float py = prt->imgs->screen.page.y;
  const float pwidth = prt->imgs->screen.page.width;
  const float pheight = prt->imgs->screen.page.height;
  const float ax = prt->imgs->screen.print_area.x;
  const float ay = prt->imgs->screen.print_area.y;
  const float awidth = prt->imgs->screen.print_area.width;
  const float aheight = prt->imgs->screen.print_area.height;

  // page w/h
  float pg_width  = prt->pinfo->paper.width;
  float pg_height = prt->pinfo->paper.height;

  // non-printable
  float np_top = prt->pinfo->printer.hw_margin_top;
  float np_left = prt->pinfo->printer.hw_margin_left;
  float np_right = prt->pinfo->printer.hw_margin_right;
  float np_bottom = prt->pinfo->printer.hw_margin_bottom;

  // handle the landscape mode if needed
  if(prt->pinfo->page.landscape)
  {
    float tmp = pg_width;
    pg_width = pg_height;
    pg_height = tmp;

    // rotate the non-printable margins
    tmp       = np_top;
    np_top    = np_right;
    np_right  = np_bottom;
    np_bottom = np_left;
    np_left   = tmp;
  }

  const float pright = px + pwidth;
  const float pbottom = py + pheight;

  // x page -> x display
  // (x / pg_width) * p_width + p_x
  cairo_set_source_rgb (cr, 0.9, 0.9, 0.9);
  cairo_rectangle (cr, px, py, pwidth, pheight);
  cairo_fill (cr);

  // display non-printable area
  cairo_set_source_rgb (cr, 0.1, 0.1, 0.1);

  const float np1x = px + (np_left / pg_width) * pwidth;
  const float np1y = py + (np_top / pg_height) * pheight;
  const float np2x = pright - (np_right / pg_width) * pwidth;
  const float np2y = pbottom - (np_bottom / pg_height) * pheight;

  // top-left
  cairo_move_to (cr, np1x-10, np1y);
  cairo_line_to (cr, np1x, np1y); cairo_line_to (cr, np1x, np1y-10);
  cairo_stroke (cr);

  // top-right
  // npy = p_y + (np_top / pg_height) * p_height;
  cairo_move_to (cr, np2x+10, np1y);
  cairo_line_to (cr, np2x, np1y); cairo_line_to (cr, np2x, np1y-10);
  cairo_stroke (cr);

  // bottom-left
  cairo_move_to (cr, np1x-10, np2y);
  cairo_line_to (cr, np1x, np2y); cairo_line_to (cr, np1x, np2y+10);
  cairo_stroke (cr);

  // bottom-right
  cairo_move_to (cr, np2x+10, np2y);
  cairo_line_to (cr, np2x, np2y); cairo_line_to (cr, np2x, np2y+10);
  cairo_stroke (cr);

  // clip to this area to ensure that the image won't be larger,
  // this is needed when using negative margin to enlarge the print

  cairo_rectangle (cr, np1x, np1y, np2x-np1x, np2y-np1y);
  cairo_clip(cr);

  cairo_set_source_rgb (cr, 0.77, 0.77, 0.77);
  cairo_rectangle (cr, ax, ay, awidth, aheight);
  cairo_fill (cr);

  return FALSE;
}

void mouse_moved(dt_view_t *self,
                 double x,
                 double y,
                 double pressure,
                 int which)
{
  const dt_print_t *prt = (dt_print_t *)self->data;

  // if we are not hovering over a thumbnail in the filmstrip -> show
  // metadata of first opened image.

  const dt_imgid_t mouse_over_id = dt_control_get_mouse_over_id();

  if(prt->imgs->count == 1
     && mouse_over_id != prt->imgs->box[0].imgid)
  {
    dt_control_set_mouse_over_id(prt->imgs->box[0].imgid);
  }
  else if(prt->imgs->count > 1)
  {
    const int bidx = dt_printing_get_image_box(prt->imgs, x, y);
    if(bidx == -1)
      dt_control_set_mouse_over_id(NO_IMGID);
    else if(mouse_over_id != prt->imgs->box[bidx].imgid)
    {
      dt_control_set_mouse_over_id(prt->imgs->box[bidx].imgid);
    }
  }
}

gboolean try_enter(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  //  now check that there is at least one selected image

  const dt_imgid_t imgid = dt_act_on_get_main_image();

  if(!dt_is_valid_imgid(imgid))
  {
    // fail :(
    dt_control_log(_("no image to open!"));
    return TRUE;
  }

  // this loads the image from db if needed:
  const dt_image_t *img = dt_image_cache_get(darktable.image_cache, imgid, 'r');
  // get image and check if it has been deleted from disk first!

  char imgfilename[PATH_MAX] = { 0 };
  gboolean from_cache = TRUE;
  dt_image_full_path(img->id, imgfilename, sizeof(imgfilename), &from_cache);
  if(!g_file_test(imgfilename, G_FILE_TEST_IS_REGULAR))
  {
    dt_control_log(_("image `%s' is currently unavailable"), img->filename);
    dt_image_cache_read_release(darktable.image_cache, img);
    return 1;
  }
  // and drop the lock again.
  dt_image_cache_read_release(darktable.image_cache, img);

  // set up load for the selected image when print_settings gets
  // signal that we've entered print view
  prt->imgs->imgid_to_load = imgid;

  return FALSE;
}

void enter(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  /* scroll filmstrip to the first selected image */
  if(dt_is_valid_imgid(prt->imgs->imgid_to_load))
  {
    // change active image
    dt_thumbtable_set_offset_image(dt_ui_thumbtable(darktable.gui->ui),
                                   prt->imgs->box[0].imgid, TRUE);
    // but no need to raise signal, as print settings is already will
    // load this image on view enter
    dt_view_active_images_reset(FALSE);
    dt_view_active_images_add(prt->imgs->imgid_to_load, FALSE);
  }

  gtk_overlay_add_overlay(GTK_OVERLAY(dt_ui_center_base(darktable.gui->ui)),
                          prt->w_main);

  // FIXME: needed?
  // be sure that log msg is always placed on top
  gtk_overlay_reorder_overlay
    (GTK_OVERLAY(dt_ui_center_base(darktable.gui->ui)),
     gtk_widget_get_parent(dt_ui_log_msg(darktable.gui->ui)), -1);
  gtk_overlay_reorder_overlay
    (GTK_OVERLAY(dt_ui_center_base(darktable.gui->ui)),
     gtk_widget_get_parent(dt_ui_toast_msg(darktable.gui->ui)), -1);

  DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals, DT_SIGNAL_VIEWMANAGER_THUMBTABLE_ACTIVATE,
                            G_CALLBACK(_view_print_filmstrip_activate_callback), self);

  gtk_widget_grab_focus(dt_ui_center(darktable.gui->ui));

  dt_control_set_mouse_over_id(prt->imgs->imgid_to_load);
}

void leave(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  /* disconnect from filmstrip image activate */
  DT_DEBUG_CONTROL_SIGNAL_DISCONNECT(darktable.signals,
                                     G_CALLBACK(_view_print_filmstrip_activate_callback),
                                     (gpointer)self);

  gtk_container_remove(GTK_CONTAINER(dt_ui_center_base(darktable.gui->ui)),
                       prt->w_main);

  dt_printing_clear_boxes(prt->imgs);
}

void gui_init(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  // view background outside of page
  GtkWidget* w_background = gtk_drawing_area_new();
  gtk_widget_set_name(w_background, "print-view-bkgd");

  // page margins
  prt->w_margins = gtk_drawing_area_new();
  gtk_widget_set_name(prt->w_margins, "print-paper-margins");

  // inner content (this widget draws the background of the paper
  // within margins, but prt->body, in a higher layer, contains the
  // layout boxes)
  prt->w_content = gtk_drawing_area_new();
  gtk_widget_set_name(prt->w_content, "print-paper-content");

  // page content & margins
  GtkWidget *w_overlay = gtk_overlay_new();
  gtk_widget_set_name(w_overlay, "print-paper-overlay");
  gtk_container_add(GTK_CONTAINER(w_overlay), prt->w_margins);
  gtk_overlay_add_overlay(GTK_OVERLAY(w_overlay), prt->w_content);

  // fits page in center area
  prt->w_aspect1 = gtk_aspect_frame_new(NULL, 0.5f, 0.5f, 1.0f, FALSE);
  gtk_widget_set_name(prt->w_aspect1, "print-paper-frame");
  gtk_container_add(GTK_CONTAINER(prt->w_aspect1), w_overlay);
#if 0
  gtk_widget_set_size_request(prt->w_aspect1, DT_PIXEL_APPLY_DPI(200),
                              DT_PIXEL_APPLY_DPI(200));
#endif

  // tick marks to show hardware margins -- this may extend beyond paper widget
  prt->w_hw_margins = gtk_drawing_area_new();
  gtk_widget_set_name(prt->w_hw_margins, "print-hw-margins");

  prt->w_page = gtk_drawing_area_new();
  gtk_widget_set_name(prt->w_page, "print-paper");

  prt->w_main = gtk_overlay_new();
  gtk_widget_set_name(prt->w_main, "print-main");
  gtk_container_add(GTK_CONTAINER(prt->w_main), w_background);
  gtk_overlay_add_overlay(GTK_OVERLAY(prt->w_main), prt->w_aspect1);
  gtk_overlay_add_overlay(GTK_OVERLAY(prt->w_main), prt->w_hw_margins);
  gtk_overlay_add_overlay(GTK_OVERLAY(prt->w_main), prt->w_page);
  gtk_overlay_add_overlay(GTK_OVERLAY(prt->w_main),
                          darktable.lib->proxy.print.w_settings_main);
  // so that margin start will be left, end will be right
  gtk_widget_set_direction(prt->w_main, GTK_TEXT_DIR_LTR);

  g_signal_connect(G_OBJECT(w_background), "draw", G_CALLBACK(_event_draw_bkgd), NULL);
  g_signal_connect(G_OBJECT(prt->w_margins), "draw", G_CALLBACK(_event_draw_rect), NULL);
  g_signal_connect(G_OBJECT(prt->w_content), "draw", G_CALLBACK(_event_draw_rect), NULL);
  g_signal_connect(G_OBJECT(prt->w_hw_margins), "draw",
                   G_CALLBACK(_event_draw_hw_margins), prt);
  g_signal_connect(G_OBJECT(prt->w_page), "draw", G_CALLBACK(_event_draw_page), prt);

  gtk_widget_show(w_background);
  gtk_widget_show(prt->w_margins);
  gtk_widget_show(prt->w_content);
  gtk_widget_show(w_overlay);
  gtk_widget_show(prt->w_aspect1);
  gtk_widget_show(prt->w_hw_margins);
  //gtk_widget_show(prt->w_page);
  gtk_widget_show(prt->w_main);

  // so that can remove it from ui center overlay when leave print
  // view, and but add it back in when re-enter print view
  g_object_ref(prt->w_main);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
