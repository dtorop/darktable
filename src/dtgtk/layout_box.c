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
/** a class to boxes containing images within a print layout */

#include "dtgtk/layout_box.h"

#include "dtgtk/thumbnail.h"
#include "gui/gtk.h"

// FIXME: add overlay with drag handles to edge of the box
// FIXME: add call-outs for box dimensions and margins (maybe means putting whole thing in a grid)
// FIXME: box should dims/hides image when dragging another one over, or alternately previews that image in the layout box

struct _dt_layout_box_t
{
  // FIXME: do we refer again to all of these? if not don't need pointers
  GtkWidget *w_draw;         // GtkDrawingArea -- layout box area
  GtkWidget *w_over;         // GtkOverlay -- contains box area, thumbnail, and eventbox
  GtkWidget *w_event;        // GtkEventBox -- handles dragging and drag-and-drop
  dt_thumbnail_t *thumb;
};


static gboolean _event_draw_box(GtkWidget *widget,
                                cairo_t *cr,
                                gpointer user_data)
{
  GtkStyleContext *context = gtk_widget_get_style_context(widget);

  guint width = gtk_widget_get_allocated_width (widget);
  guint height = gtk_widget_get_allocated_height (widget);

  printf("draw box event %dx%d\n", width, height);
  gtk_render_background(context, cr, 0, 0, width, height);
  gtk_render_frame(context, cr, 0, 0, width, height);

  // test X across page
  cairo_set_source_rgb (cr, 0.8, 0.3, 0.3);
  cairo_move_to(cr, 0, 0);
  cairo_line_to(cr, width, height);
  cairo_stroke(cr);
  cairo_move_to(cr, width, 0);
  cairo_line_to(cr, 0, height);
  cairo_stroke(cr);

  // should propagate?
  return FALSE;
}

#if 0
static gboolean _mouse_moved_box(GtkWidget *w, GdkEventMotion *event, dt_print_layout_t *layout)
{
  printf("mouse moved box\n");
  return FALSE;
}
#endif

#if 0
// FIXME: should this get pos_screen or a dt_image_box?
// FIXME: should this be dt_print_layout_box_full_page() and the rest happens via gui?
dt_layout_box_t *dt_layout_box_new(const dt_image_box *const box,
                                   const dt_imgid_t imgid)
{
  //printf("dt_print_layout_new_box screen %f x %f + %f x %f imgid %d\n", box->screen.x, box->screen.y, box->screen.width, box->screen.height, imgid);
  //printf(" rel %f x %f + %f x %f\n", box->pos.x, box->pos.y, box->pos.width, box->pos.height);

  const int paper_width_px = gtk_widget_get_allocated_width(layout->w_paper);
  const int paper_height_px = gtk_widget_get_allocated_height(layout->w_paper);

  // screen coordinates of bounding box
  // FIXME: if these are only used once, put inline below
  const int box_x = paper_width_px * box->pos.x;
  const int box_y = paper_height_px * box->pos.y;
  const int box_width = paper_width_px * box->pos.width;
  const int box_height = paper_height_px * box->pos.height;
  //printf("box in screen coords from rel %d x %d + %d x %d\n", box_x, box_y, box_width, box_height);

  dt_layout_box_t *b =
    _new_empty_box(layout, box_x, box_y, box_width, box_height);
  // FIXME: need this?
  memcpy(&b->box, box, sizeof(dt_image_box));

  // FIXME: have a way to make a thumbnail active, and have the active one be controlled from right panel
  // FIXME: should return a dt_thumbnail_t or something less transparent, e.g. an image # of abstract pointer to image list?

  // FIXME: make new container type DT_THUMBNAIL_CONTAINER_PRINT or DT_THUMBNAIL_CONTAINER_PRINT_LAYOUT which is like _PREVIEW but not zoomable
  // FIXME: should rowid be 0?
  b->thumb = dt_thumbnail_new(box_width, box_height, IMG_TO_FIT, imgid, -1,
                              DT_THUMBNAIL_OVERLAYS_NONE,
                              DT_THUMBNAIL_CONTAINER_PREVIEW, FALSE);
  gtk_overlay_add_overlay(GTK_OVERLAY(b->w_over), b->thumb->w_main);
  // FIXME: handle alignment of thumbnail within the box

  return b;
}
#endif

void dt_layout_box_resize(dt_layout_box_t *box, gint width, gint height)
{
  gtk_widget_set_size_request(box->w_over, width, height);
  // FIXME: do want to set force to TRUE?
  // FIXME: do need to call thumbnail resize?
  if(box->thumb)
    dt_thumbnail_resize(box->thumb, width, height, FALSE, IMG_TO_FIT);
}

void dt_layout_box_set_imgid(dt_layout_box_t *box, const dt_imgid_t imgid)
{
  printf("dt_layout_box_set_imgid: to img %d\n", imgid);
  if(box->thumb)
  {
    printf("dt_layout_box_set_imgid: destroying thumbnail\n");
    dt_thumbnail_destroy(box->thumb);
  }
  box->thumb = dt_thumbnail_new(gtk_widget_get_allocated_width(box->w_over),
                                gtk_widget_get_allocated_height(box->w_over),
                                IMG_TO_FIT, imgid, -1,
                                DT_THUMBNAIL_OVERLAYS_NONE,
                                DT_THUMBNAIL_CONTAINER_PREVIEW, FALSE);
  gtk_overlay_add_overlay(GTK_OVERLAY(box->w_over), box->thumb->w_main);

  // event box needs to be in front of thumbnail to get priority for drag-and-drop
  // FIXME: this won't matter if we can just change image on extant thumbnail
  gtk_overlay_reorder_overlay(GTK_OVERLAY(box->w_over), box->thumb->w_main, 0);
}

GtkWidget *dt_layout_box_widget(dt_layout_box_t *box)
{
  return box->w_over;
}

dt_layout_box_t *dt_layout_box_new(double x, double y, double width, double height,
                                   const dt_imgid_t imgid)
{
  // FIXME: use some of code above to handle if imgid is not NO_IMGID

  // FIXME: it would be nice if we could just attach box dimensions as user data connected to each widget, or send a configure event to each child and have the configure handler already know about the child's box data
  dt_layout_box_t *b = (dt_layout_box_t *)calloc(1, sizeof(dt_layout_box_t));

  //GtkWidget *layout_box = _new_layout_box(dt_print_layout_t *const layout, box_x, box_y, box_width, box_height);
  b->w_draw = gtk_drawing_area_new();
  dt_gui_add_class(b->w_draw, "print-layout-box");
  gtk_widget_show(b->w_draw);

  b->w_event = gtk_event_box_new();
  gtk_widget_show(b->w_event);

  // FIXME: even if use this, event box should wrap overlay so gets drag-and-drop
  // FIXME: do need overlay once get rid of overall box? or need it for handles for resizing?
  b->w_over = gtk_overlay_new();
  gtk_widget_set_size_request(b->w_over, width, height);
  gtk_widget_show(b->w_over);
  gtk_container_add(GTK_CONTAINER(b->w_over), b->w_draw);
  // FIXME: if we new thumb can to NO_IMGID then we can just create a blank thumbnail here
  gtk_overlay_add_overlay(GTK_OVERLAY(b->w_over), b->w_event);

  // FIXME: pass in layout or just info on the box?
  g_signal_connect(G_OBJECT(b->w_draw), "draw", G_CALLBACK(_event_draw_box), NULL);

  // FIXME: really want to make mouse move also resize, so need to detect where the click is, or set up cells which catch different events
  // FIXME: when enter layout box (or click on it?) should dim other boxes, show drag handles, and show measurements overlay

#if 0
  g_signal_connect(G_OBJECT(b->w_over), "motion-notify-event", G_CALLBACK(_mouse_moved_box), layout);
#endif

  return b;
}

void dt_layout_box_destroy(dt_layout_box_t *box)
{
  free(box);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
