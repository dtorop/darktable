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

#include "dtgtk/layout_body.h"

#include "dtgtk/layout_box.h"
#include "gui/drag_and_drop.h"
#include "gui/gtk.h"

// FIXME: when flip from portrait to landscape and back the boxes drift
// FIXME: experiment with making individual images active/inactive when they're clicked or moused over
// FIXME: when enter view with one image active, make a box containing that image at the full size of the body
// FIXME: add a spacing paramter, when enter view with >1 image active, tile layout boxes for each image with spacing gap between each, filling the body
// FIXME: similarly, when drag images into an empty body, make layout boxes for each, tile/space as necessary
// FIXME: similarly, when drag > 1 image into a layout box, split it into multiple layout boxes, tiling/spacing appropriately
// FIXME: when double-click on the only layout box, expand it to fit the page body
// FIXME: when double-click on one of several layout boxes, tile them all evenly according to the spacing paramter

#define MIN_BOX_SIZE 20

typedef struct _box_info_t
{
  double rel_x, rel_y, rel_width, rel_height;
  // s/data/box/
  dt_layout_box_t *data;
  // FIXME: needed?
  dt_layout_body_t *parent;
} _box_info_t;

struct _dt_layout_body_t
{
  DTPageDsc *page_dsc;

  // FIXME: s/w_overlay/w_body/ or maybe w_main
  GtkWidget *w_overlay;      // GtkOverlay -- main widget, contains paper body,
                             //   event box, and select box overlay
  GtkWidget *w_fixed;        // GtkFixed -- inside w_overlay, contains layout boxes
  GtkWidget *w_select_box;   // GtkDrawingArea -- box when creating new layout

  GSList *boxes;

  _box_info_t *drag_box;                // box being dragged, NULL if none
  double box_initial_x, box_initial_y;  // screen coords at drag start

  double select_x, select_y, select_width, select_height;

  // FIXME: if we can get drag data we don't need to remember this -- and then we can also drag in multiple images -- and reject drag if dragging multiple images to one box
  dt_imgid_t last_activated_img;
};

static gboolean _event_size_allocate_fixed(GtkWidget *w_fixed,
                                           GtkAllocation *alloc,
                                           dt_layout_body_t *body)
{
  // FIXME: doc of GtkWindow says that the event allocation may include client decorations, gtk_window_get_size() is more reliable, but unfortunately fixed or event box don't seem to have widgets, and configure-event doesn't seem to be called on any of these widgets
  // FIXME: create a GtkDrawingArea which does produce configure events, so even if it is empty this could be more reliable than size-allocate?
  printf("fixed size-allocate event x %d y %d width %d height %d\n", alloc->x, alloc->y, alloc->width, alloc->height);

  for(GSList *boxes = body->boxes; boxes; boxes = g_slist_next(boxes))
  {
    _box_info_t *binfo = boxes->data;
    const gint box_x = alloc->width * binfo->rel_x;
    const gint box_y = alloc->height * binfo->rel_y;
    const gint box_width = alloc->width * binfo->rel_width;
    const gint box_height = alloc->height * binfo->rel_height;
    printf("box update to %d x %d + %d x %d\n", box_x, box_y, box_width, box_height);
    GtkWidget *w = dt_layout_box_widget(binfo->data);
    gtk_fixed_move(GTK_FIXED(w_fixed), w, box_x, box_y);
    dt_layout_box_resize(binfo->data, box_width, box_height);
  }

  return FALSE;
}

// FIXME: be more specific naming these
static void _box_drag_gesture_begin(GtkGestureDrag *gesture,
                                    gdouble offset_x, gdouble offset_y,
                                    dt_layout_body_t *body)
{
  printf("_box_drag_gesture_begin drag_location %p\n", body->drag_box);
  if(body->drag_box)
  {
    gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);

    const int content_width = gtk_widget_get_allocated_width(body->w_fixed);
    const int content_height = gtk_widget_get_allocated_height(body->w_fixed);
    body->box_initial_x = content_width * body->drag_box->rel_x;
    body->box_initial_y = content_height * body->drag_box->rel_y;
    printf("content %d x %d box rel %f x %f initial %f,%f\n", content_width, content_height, body->drag_box->rel_x, body->drag_box->rel_y, body->box_initial_x, body->box_initial_y);

#if 0
    // move widget to top
    GtkWidget *w = dt_layout_box_widget(body->drag_box->data);
    gtk_container_remove(GTK_CONTAINER(body->w_fixed), g_object_ref(w));
    gtk_fixed_put(GTK_FIXED(body->w_fixed), w,
                  body->box_initial_x, body->box_initial_y);
    g_object_unref(w);
#endif

    // FIXME: should the cursor always be "move" when over draggable part of object?
    GdkCursor *const cursor = gdk_cursor_new_from_name(gdk_display_get_default(), "move");

    //gdk_window_set_cursor(gtk_widget_get_window(layout->drag_box->w_over), cursor);
    gdk_window_set_cursor(gtk_widget_get_window(body->w_overlay), cursor);
    g_object_unref(cursor);
  }
  else
  {
    gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_DENIED);
  }
}

static void _box_drag_gesture_update(GtkGestureDrag *gesture,
                                     gdouble offset_x, gdouble offset_y,
                                     dt_layout_body_t *body)
{
  printf("_box_drag_gesture_update\n");
  if(body->drag_box)
  {
    // FIXME: should constrain movement so that it doesn't go outside of area w/in margins
    double new_x = body->box_initial_x + offset_x;
    double new_y = body->box_initial_y + offset_y;
    gtk_fixed_move(GTK_FIXED(body->w_fixed),
                   dt_layout_box_widget(body->drag_box->data),
                   (gint)round(new_x), (gint)round(new_y));
    const int content_width = gtk_widget_get_allocated_width(body->w_fixed);
    const int content_height = gtk_widget_get_allocated_height(body->w_fixed);
    // FIXME: need to update box.screen and box.print? or don't even need to update any of this until drag is finished?
    body->drag_box->rel_x = (float) new_x / content_width;
    body->drag_box->rel_y = (float) new_y / content_height;
    // FIXME: we should trigger callback to update coords into right panel?
    printf("new %f,%f content %d x %d box rel %f x %f initial %f,%f\n", new_x, new_y, content_width, content_height, body->drag_box->rel_x, body->drag_box->rel_y, body->box_initial_x, body->box_initial_y);
  }
}

static void _box_drag_gesture_end(GtkGestureDrag *gesture,
                                  gdouble offset_x, gdouble offset_y,
                                  dt_layout_body_t *body)
{
  printf("_box_drag_gesture_end\n");
  if(body->drag_box)
  {
    //gdk_window_set_cursor(gtk_widget_get_window(layout->drag_box->w_over), NULL);
    gdk_window_set_cursor(gtk_widget_get_window(body->w_overlay), NULL);
    // FIXME: is this necessary?
    body->drag_box = NULL;
  }
}

static gboolean _event_draw_select_box(GtkWidget *widget, cairo_t *cr,
                                       dt_layout_body_t *body)
{
  GtkStyleContext *const context = gtk_widget_get_style_context(widget);
  double x = body->select_x;
  double y = body->select_y;
  double width = body->select_width;
  double height = body->select_height;

  // unlike cairo GTK render functions don't handle negative dimensions
  if(width < 0)
  {
    x = x + width;
    width = -width;
  }
  if(height < 0)
  {
    y = y + height;
    height = -height;
  }

  // FIXME: should use CSS outline property instead?
  gtk_render_background(context, cr, x, y, width, height);
  gtk_render_frame(context, cr, x, y, width, height);

  return FALSE;
}

#if 0
static void _fixed_gesture_pressed(GtkGestureMultiPress *gesture,
                                   guint n_press, gdouble x, gdouble y,
                                   gpointer user_data)
{
  GtkWidget *w = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
  printf("_fixed_gesture_pressed widget %p\n", w);
}
#endif

static void _box_gesture_pressed(GtkGestureMultiPress *gesture,
                                 guint n_press, gdouble x, gdouble y,
                                 _box_info_t *binfo)
{
  // press gesture on a box widget records which one we want to move,
  // this is followed by drag gesture in the parent widget to move
  // that box
  GtkWidget *w = gtk_event_controller_get_widget(GTK_EVENT_CONTROLLER(gesture));
  printf("_box_gesture_pressed widget %p binfo %p\n", w, binfo);
  // FIXME: hook this up better
  // FIXME: is this a case to set the parent widget data, or can we select/activate this widget, which makes us know where to look for it?
  // FIXME: hacky, figure out another way to set this, e.g. this callback could take to user data params?
  binfo->parent->drag_box = binfo;
  // FIXME: it would be nice to move the widget to the top here, but removing and replacing it in the fixed cancels the ongoing drag gesture -- is there another way? current way dragging it moves it to top but just clicking on it doesn't. Maybe there's another action that moves widget to top (control-click, etc.) and that dragging doesn't change this order
  // FIXME: should activate/focus the widget? see https://docs.gtk.org/gtk3/class.DrawingArea.html on gtk_widget_has_focus() and gtk_render_focus()
  // fIXME: focus should trigger a callback to load coords into right panel? should mouseover of widget do that?
  // FIXME: this is where we should trigger callback to load coords into right panel? or do we select the widget, which then leads to loading coordinates?
}

static gboolean _box_drag_motion_received(GtkWidget *widget,
                                          GdkDragContext *dc,
                                          const gint x,
                                          const gint y,
                                          const guint time,
                                          gpointer user_data)
{
  // FIXME: can use gtk_widget_set_state_flags()?
  dt_gui_add_class(widget, "print-drag-target");

  printf("_box_drag_motion_received %p %d,%d\n", widget, x, y);

  return TRUE;
}

static void _box_drag_motion_leave(GtkWidget *widget,
                                   GdkDragContext *dc,
                                   const guint time,
                                   gpointer user_data)
{
  // FIXME: should highlight the widget
  printf("_box_drag_motion_leave %p\n", widget);
  dt_gui_remove_class(widget, "print-drag-target");
#if 0
  // FIXME: this cancels the drag, instead get source window from drag context and catch "drag-failed" or catch cancel from the DragContext
  gtk_widget_show(layout->w_hw_margins);
#endif
}

static void _box_drag_and_drop_received(GtkWidget *widget,
                                        GdkDragContext *context,
                                        gint x,
                                        gint y,
                                        GtkSelectionData *selection_data,
                                        guint target_type,
                                        guint time,
                                        _box_info_t *binfo)
{
  printf("_box_drop_received widget %p binfo %p %d,%d\n", widget, binfo, x, y);

  // FIXME: can we figure this out from drag-and-drop data instead?
  const dt_imgid_t imgid = binfo->parent->last_activated_img;
  // FIXME: if we have to do this, do this at page layout level
  //gtk_widget_show(layout->w_hw_margins);
  if(dt_is_valid_imgid(imgid))
  {
    dt_layout_box_set_imgid(binfo->data, imgid);
  }
}

static void _new_box_drag_gesture_begin(GtkGestureDrag *gesture,
                                    gdouble offset_x, gdouble offset_y,
                                    dt_layout_body_t *body)
{
  printf("_new_box_drag_gesture_begin %f,%f\n", offset_x, offset_y);

  // FIXME: only claim if within margins, or make this eventbox be sized to fit in margins
  gtk_gesture_set_state(GTK_GESTURE(gesture), GTK_EVENT_SEQUENCE_CLAIMED);

  // FIXME: snap to grid
  body->select_x = offset_x;
  body->select_y = offset_y;
  body->select_width = body->select_height = 0.0;
  gtk_widget_show(body->w_select_box);

  GdkCursor *const cursor = gdk_cursor_new_from_name(gdk_display_get_default(), "crosshair");
  gdk_window_set_cursor(gtk_widget_get_window(body->w_overlay), cursor);
  g_object_unref(cursor);
}

static void _new_box_drag_gesture_update(GtkGestureDrag *gesture,
                                         gdouble offset_x, gdouble offset_y,
                                         dt_layout_body_t *body)
{
  printf("_new_box_drag_gesture_update %f,%f\n", offset_x, offset_y);
  // FIXME: if go outside page margins, clamp offset
  // FIXME: snap to grid
  body->select_width = offset_x;
  body->select_height = offset_y;
  gtk_widget_queue_draw(body->w_select_box);
}

// FIXME: rename to _new_box_create?
static void _new_box_drag_gesture_end(GtkGestureDrag *gesture,
                                      gdouble offset_x, gdouble offset_y,
                                      dt_layout_body_t *body)
{
  printf("_new_box_drag_gesture_end %f,%f\n", offset_x, offset_y);
  gdk_window_set_cursor(gtk_widget_get_window(body->w_overlay), NULL);
  gtk_widget_hide(body->w_select_box);

  double x = body->select_x;
  double y = body->select_y;
  double width = body->select_width;
  double height = body->select_height;

  // unlike cairo GTK render functions don't handle negative dimensions
  if(width < 0)
  {
    x = x + width;
    width = -width;
  }
  if(height < 0)
  {
    y = y + height;
    height = -height;
  }

  // FIXME: snap to grid

  if(width < DT_PIXEL_APPLY_DPI(MIN_BOX_SIZE) ||
     height < DT_PIXEL_APPLY_DPI(MIN_BOX_SIZE))
    return;

  _box_info_t *binfo = calloc(1, sizeof(_box_info_t));
  binfo->data = dt_layout_box_new(x, y, width, height, NO_IMGID);
  binfo->parent = body;

  // FIXME: these measurements don't work when the aspect ratio changes -- fix that!
  const int content_width = gtk_widget_get_allocated_width(body->w_fixed);
  const int content_height = gtk_widget_get_allocated_height(body->w_fixed);
  // FIXME: is this the same as w_fixed allocation? if so just use that
  binfo->rel_x = x / content_width;
  binfo->rel_y = y / content_height;
  binfo->rel_width = width / content_width;
  binfo->rel_height = height / content_height;

  GtkWidget *w = dt_layout_box_widget(binfo->data);
  printf("set up box at %f,%f %fx%f binfo %p data %p widget %p\n", x, y, width, height, binfo, binfo->data, w);

  // we handle draggingg in the printable area widget, but detecting a
  // click in the layout box helps us to know when a layout widget is
  // being dragged
  GtkGesture *g = gtk_gesture_multi_press_new(w);
  // FIXME: GTK complains when we do this -- do we need to do this?
  //gtk_gesture_group(layout->drag_gesture, g);
  g_signal_connect(g, "pressed", G_CALLBACK(_box_gesture_pressed), binfo);
  gtk_drag_dest_set(w, GTK_DEST_DEFAULT_ALL,
                    // FIXME: should be target_list_internal?
                    target_list_all, n_targets_all, GDK_ACTION_MOVE);
  g_signal_connect(w, "drag-motion", G_CALLBACK(_box_drag_motion_received), NULL);
  g_signal_connect(w, "drag-leave", G_CALLBACK(_box_drag_motion_leave), NULL);
  g_signal_connect(w, "drag-data-received",
                   G_CALLBACK(_box_drag_and_drop_received), binfo);
  // FIXME: when drag fails need to show hw margins
  
  body->boxes = g_slist_append(body->boxes, binfo);

  gtk_fixed_put(GTK_FIXED(body->w_fixed), w, x, y);
}

#if 0
static gboolean _mouse_moved_fixed(GtkWidget *w, GdkEventMotion *event, gpointer user_data)
{
  printf("mouse moved thumbnail\n");
  return FALSE;
}
#endif

// FIXME: same function in layout_page -- if keep both then make a GtkDrawingAreaPassthrough sub-class widget
static void _event_realize_pass_through(GtkWidget *widget, gpointer user_data)
{
  // this drawing area blocks events going to the layout, and as it
  // has a GDK window, we must wait for that window to be realized
  // then set that window to pass-through
  // NOTE: in GTK 4 we can move to set_can_target()
  GdkWindow *window = gtk_widget_get_window(widget);
  gdk_window_set_pass_through(window, TRUE);
}

static void _thumbtable_activate_callback(gpointer instance,
                                          const dt_imgid_t imgid,
                                          dt_layout_body_t *body)
{
  // FIXME: should load image if there's a single layout box
  if(dt_is_valid_imgid(imgid))
  {
    // FIXME: once get data from drag-and-drop, won't need to save this
    body->last_activated_img = imgid;
    printf("last selected is %d\n", imgid);
  }
}

GtkWidget *dt_layout_body_get_widget(dt_layout_body_t *body)
{
  return body->w_overlay;
}

dt_layout_body_t *dt_layout_body_new(DTPageDsc *page_dsc)
{
  printf("dt_layout_body_new called\n");
  dt_layout_body_t *body = (dt_layout_body_t *)calloc(1, sizeof(dt_layout_body_t));

  body->page_dsc = g_object_ref(page_dsc);
  
  body->w_fixed = gtk_fixed_new();
  gtk_widget_set_name(body->w_fixed, "print-fixed");
  dt_gui_add_class(body->w_fixed, "dt_overlays_none");
  dt_gui_add_class(body->w_fixed, "dt_fullview");
  dt_gui_add_class(body->w_fixed, "dt_printlayout");

  // to catch events when dragging in empty paper area
  GtkWidget *event_box = gtk_event_box_new();
  gtk_widget_set_name(event_box, "print-eb-for-fixed");
  gtk_container_add(GTK_CONTAINER(event_box), body->w_fixed);

  body->w_select_box = gtk_drawing_area_new();
  gtk_widget_set_name(body->w_select_box, "print-select-box");

  body->w_overlay = gtk_overlay_new();
  gtk_widget_set_name(body->w_overlay, "print-over-paper-and-fixed");
  gtk_container_add(GTK_CONTAINER(body->w_overlay), event_box);
  gtk_overlay_add_overlay(GTK_OVERLAY(body->w_overlay), body->w_select_box);
  gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(body->w_overlay),
                                       body->w_select_box, TRUE);
  g_signal_connect(G_OBJECT(body->w_select_box), "realize",
                   G_CALLBACK(_event_realize_pass_through), NULL);

  g_signal_connect(G_OBJECT(body->w_fixed), "size-allocate",
                   G_CALLBACK(_event_size_allocate_fixed), body);
  g_signal_connect(G_OBJECT(body->w_select_box), "draw",
                   G_CALLBACK(_event_draw_select_box), body);

  // handle dragging to make a new layout box
#if 0
  GtkGesture *g = gtk_gesture_multi_press_new(event_box);
  //gtk_gesture_group(body->drag_gesture, g);
  g_signal_connect(g, "pressed",
                   G_CALLBACK(_fixed_gesture_pressed), body);
#endif
  GtkGesture *g_event = gtk_gesture_drag_new(event_box);
  g_signal_connect(g_event, "drag-begin",
                   G_CALLBACK(_new_box_drag_gesture_begin), body);
  g_signal_connect(g_event, "drag-update",
                   G_CALLBACK(_new_box_drag_gesture_update), body);
  g_signal_connect(g_event, "drag-end",
                   G_CALLBACK(_new_box_drag_gesture_end), body);
  
  // handle dragging via the fixed container, which gives us a frame
  // of reference (rather than dragging a widget relative to itself),
  // but catch the button press which inititates the drag per-widget,
  // to know which widget is being dragged
  body->drag_box = NULL;
  // this only catches events which are in layout boxes, as w_fixed
  // doesn't have its own window, so it's convenient for handling
  // dragging of existing layouts only
  // FIXME: but we could in the drag-begin function differentiate if this is a new or existing box and reject/accept gesture depending
  GtkGesture *g_fixed = gtk_gesture_drag_new(body->w_fixed);
  g_signal_connect(g_fixed, "drag-begin", G_CALLBACK(_box_drag_gesture_begin), body);
  g_signal_connect(g_fixed, "drag-update", G_CALLBACK(_box_drag_gesture_update), body);
  g_signal_connect(g_fixed, "drag-end", G_CALLBACK(_box_drag_gesture_end), body);

  // user activated a new image via the filmstrip or user entered view
  // mode which activates an image
  DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals,
                                  DT_SIGNAL_VIEWMANAGER_THUMBTABLE_ACTIVATE,
                                  G_CALLBACK(_thumbtable_activate_callback),
                                  body);

  gtk_widget_show(body->w_fixed);
  gtk_widget_show(event_box);
  gtk_widget_show(body->w_overlay);

  printf("dt_layout_body_new created body %p\n", body);
  return body;
}

void dt_layout_body_destroy(dt_layout_body_t *const body)
{
  // this will be called from dt_layout_page_destroy which, by
  // destroying its widget which contains this one, will take care of
  // unref'ing our widgets and signal handlers

  // FIXME: how to nicely deallocate boxes
  //g_slist_free_full(layout->boxes, free);
  while(body->boxes)
  {
    _box_info_t *binfo = body->boxes->data;
    body->boxes = g_slist_remove(body->boxes, binfo);
    dt_layout_box_destroy(binfo->data);
    free(binfo);
  }
  
  g_object_unref(body->page_dsc);
  free(body);
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
