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
#include "dtgtk/layout_page.h"
#include "dtgtk/thumbtable.h"
#include "gui/accelerators.h"
#include "gui/gtk.h"
#include "gui/presets.h"
#include "print/page_dsc.h"
#include "views/view.h"
#include "views/view_api.h"

#include <gdk/gdkkeysyms.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

DT_MODULE(1)

// planning for MVC
//
// we need a model, maybe here or in common/printing.c which contains
//   printer name
//   paper name and dimensions
//   media type
//   hardware margins (yes?)
//   paper orientation
//   profile/intent for conversion to intermediate colorspace, style, and style mode
//   profile/intent for conversion to printe colorspace + black point compensation flag
//   a list of image boxes
//     each with origin, dimensions and alignment (horizontal/vertical start/middle/end)
//   whether in borderless model
//   whether to snap to grid, display grid, and grid intervals   (does this pertain to model?)
//   display units for dimensions (could be in a separate model with grid info?)
// all measurements can be in mm
//
// The catch is the print_settings.c needs to be able to save this
// model as a preset.  So do we keep a canonical gobject model in
// printing.c with signals, and print_settings.c updates itself or
// updates the model as it changes?
//
// Print settings just needs to be apprised of all the current
// settigngs to save them in get_params, it coule be bound to
// them. And it just sets GUI in set_params, which will set model and
// by bindings/signals can set print layout properly.
//
// we have two views
//
// 1.) print settings panel
//     which views (from the model)
//       printer, paper size, media type
//       paper orientation
//       profile/intent for conversion to intermediate colorspace
//       profile/intent for conversion to printe colorspace + black point compensation flag
//       current image box origin, dimensions and alignment (horizontal/vertical start/middle/end)
//       current image box's image dimensions (this is updated in model by center view?)
//       current image's scale factor and dpi (is this calculated here?)
//       grid intervals
//       whether in borderless mode
//     which allows for controlling (which alters the model)
//       printer, paper size, media type
//       paper orientation
//       profile/intent for conversion to intermediate colorspace, style, and style mode
//       profile/intent for conversion to printe colorspace + black point compensation flag
//       current image box origin, dimensions and alignment (horizontal/vertical start/middle/end)
//       making a new image box (is this necessary, or just drag in center view?)
//       deleting current image box
//       deleting all image boxes
//       grid intervals
//       whether to snap to grid
//       whether to display grid
//
// 2.) center view
//     which views
//        paper size (as a ratio)
//        hardware margins
//        paper orientation
//        all image boxes origin/dimensions/alignment
//        all images within boxes origin/dimensions
//        knockouts for layout margins/dimensions when mouseover an image
//        knockouts for image margins/dimensions would also be nice
//        grid
//      which controls
//        paper orientation (when create a single image or drag image into a single layout box -- does this happen in views/print.c?)
//        current image box origin/dimensions
//        imgid in box (when drag-and-drop image onto a box)
//
// make gobjects to hold the model
// TODO: can use g_value_register_transform_func() to convert between mm/cm/in units? or is that more for type casting?
// TODO: actually looks like the right way is to use a binding https://docs.gtk.org/gobject/class.Binding.html
// TODO: set up properties via g_param_spec_double() or g_param_spec_int() or g_param_spec_pointer()
// TODO: use g_value_get_int() and ilk to work with gobject via
// TODO: how do use GParamSpec to work on GObject properties? via g_object_interface_install_property() of derived class? via GTypeInfo? see https://developer-old.gnome.org/gobject/unstable/gobject-properties.html for how to set up a gobject subclass with properties
// TODO: can use https://docs.gtk.org/gobject/method.Object.notify_by_pspec.html to notify within classes
// TODO: use view's init() function to set up GObjects, then the gui_init() funcs set up callbacks/bindings
// TODO: start with landscape as the first trial one
// TODO: think about the different common files which handle printing, and perhaps make a src/pront directory to include s/cups_print/cups/, s/printprof/profile/, and printing.h (which may be obviated by layout_page) and any gobject clsses to encapsulate dt_print_info_t and dt_images_box
// TODO: do we want to register a new class with properties for storage or a new type https://docs.gtk.org/gobject/func.type_register_static.html
// TODO: can the instance init/constructor function do some helpful work setting up relationships? (see GObject concepts), maybe class_init/instance_init
// TODO: "One of GObject’s nice features is its generic get/set mechanism for object properties. When an object is instantiated, the object’s class_init handler should be used to register the object’s properties with g_object_class_install_properties()." (from GObject concepts)
// TODO: use signal detail parameter on notify to figure out which properties are changed
// TODO: if there's only one instance of layout/print data, should attach an instance signal connection or connect it via the class?

typedef struct dt_print_t
{
  // pinfo and imgs are pointers to the actual data in
  // dt_lib_print_settings_t and are set by _view_print_settings()
  // aka dt_view_print_settings()
  //
  // model: printer/page/paper/medium
  // FIXME: this is used for now obsolete routines to draw the page here, and for sending its components to layout_page creation -- but maybe layout_page can just catch when its components are set, and otherwise be created empty/unsized? then don't have to store this here at all!
  dt_print_info_t *pinfo;
  // model: image layout boxes and full page dimensions
  // FIXME: this overlays a bit with pinfo
  dt_images_box *imgs;
  dt_imgid_t last_selected;

  // view: page size/ratio, image boxes, images
  // control: image box position/size
  dt_layout_page_t *page;

  DTPageDsc *page_dsc;
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

// FIXME: do not need to call this if thumbnail is handling mipmap updates?
static void _print_mipmaps_updated_signal_callback(gpointer instance,
                                                   dt_imgid_t imgid,
                                                   gpointer user_data)
{
  printf("_print_mipmaps_updated_signal_callback()\n");
  dt_control_queue_redraw_center();
}

static void _film_strip_activated(const dt_imgid_t imgid, void *data)
{
  const dt_view_t *self = (dt_view_t *)data;
  dt_print_t *prt = (dt_print_t *)self->data;

  prt->last_selected = imgid;

  // only select from filmstrip if there is a single image displayed, otherwise
  // we will drag and drop into different areas.

  if(prt->imgs->count != 1) return;

  // if the previous shown image is selected and the selection is unique
  // then we change the selected image to the new one
  if(dt_is_valid_imgid(prt->imgs->box[0].imgid))
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
      if(sqlite3_column_int(stmt, 0) == prt->imgs->box[0].imgid
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

  prt->imgs->box[0].imgid = imgid;

  dt_thumbtable_set_offset_image(dt_ui_thumbtable(darktable.gui->ui), imgid, TRUE);

  // update the active images list
  g_slist_free(darktable.view_manager->active_images);
  darktable.view_manager->active_images = g_slist_prepend(NULL, GINT_TO_POINTER(imgid));
  DT_DEBUG_CONTROL_SIGNAL_RAISE(darktable.signals, DT_SIGNAL_ACTIVE_IMAGES_CHANGE);

  // force redraw
  dt_control_queue_redraw();
}

static void _view_print_filmstrip_activate_callback(gpointer instance,
                                                    dt_imgid_t imgid,
                                                    gpointer user_data)
{
  if(dt_is_valid_imgid(imgid)) _film_strip_activated(imgid, user_data);
}

#if 0
static void _view_new_layout_box(const dt_view_t *view,
                                 const dt_image_box *const box,
                                 const dt_imgid_t imgid)
{
  printf("_view_new_layout_box %fx%f + %fx%f imgid %d\n", box->screen.x, box->screen.y, box->screen.width, box->screen.height, imgid);
  dt_print_t *prt = (dt_print_t *)view->data;
  // FIXME: should this call into print_settings to update dt_images_box there? or should there be some sort of callback/binding
  if(prt->print_layout)
    dt_layout_box_new(box, imgid);
  else
    dt_print(DT_DEBUG_ALWAYS, "print.c: _view_new_layout_box called before layout created\n");    
}
#endif

#if 0
// FIXME: how to implement this? can the UI know when to destroy the active layout box and handle it directly
void _view_destroy_layout_box(const dt_view_t *view, dt_thumbnail_t *const box)
{
  //dt_print_t *prt = (dt_print_t *)view->data;
  dt_thumbnail_destroy(box);
}
#endif

// FIXME: this is a catchall for any print settings updated by libs/print_settings.c, instead have specific callbacks for specific updates and store more state here, and respond to more events via widgets here
// FIXME: at the least there could be one call for pinfo changed (e.g. orientation) and one for imgs changed (e.g. new image)
static void _view_print_settings(const dt_view_t *view,
                                 dt_print_info_t *pinfo,
                                 dt_images_box *imgs)
{
  dt_print_t *prt = (dt_print_t *)view->data;

  prt->pinfo = pinfo;
  prt->imgs = imgs;
  dt_control_queue_redraw();
}

#if 0
static void _drag_and_drop_received(GtkWidget *widget,
                                    GdkDragContext *context,
                                    gint x,
                                    gint y,
                                    GtkSelectionData *selection_data,
                                    guint target_type,
                                    guint time,
                                    gpointer data)
{
  const dt_view_t *self = (dt_view_t *)data;
  dt_print_t *prt = (dt_print_t *)self->data;

  const int bidx = dt_printing_get_image_box(prt->imgs, x, y);

  if(bidx != -1)
    dt_printing_setup_image(prt->imgs, bidx, prt->last_selected,
                            100, 100, ALIGNMENT_CENTER);

  prt->imgs->motion_over = -1;
  dt_control_queue_redraw_center();
}

static gboolean _drag_motion_received(GtkWidget *widget,
                                      GdkDragContext *dc,
                                      const gint x,
                                      const gint y,
                                      const guint time,
                                      gpointer data)
{
  const dt_view_t *self = (dt_view_t *)data;
  dt_print_t *prt = (dt_print_t *)self->data;

  const int bidx = dt_printing_get_image_box(prt->imgs, x, y);
  prt->imgs->motion_over = bidx;

  if(bidx != -1) dt_control_queue_redraw_center();

  return TRUE;
}
#endif

void
init(dt_view_t *self)
{
  self->data = calloc(1, sizeof(dt_print_t));
  dt_print_t *prt = (dt_print_t *)self->data;

  // FIXME: init this when enter view? it would make disconnecting signal handlers from this easier
  prt->page_dsc = dt_page_dsc_new();
  // FIXME: is this right way to "copy" a pinter?
  darktable.view_manager->proxy.print.page_dsc = g_object_ref(prt->page_dsc);

  /* initialize CB to get the print settings from corresponding lib module */
  darktable.view_manager->proxy.print.view = self;
  //darktable.view_manager->proxy.print.new_layout_box = _view_new_layout_box;
  //darktable.view_manager->proxy.print.destroy_layout_box = _view_destroy_layout_box;
  darktable.view_manager->proxy.print.print_settings = _view_print_settings;
}

void cleanup(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t *)self->data;

  // FIXME: why doesn't this work?
  //g_clear_object(darktable.view_manager->proxy.print.page_dsc);
  //g_clear_object(prt->page_dsc);
  g_object_unref(darktable.view_manager->proxy.print.page_dsc);
  g_object_unref(prt->page_dsc);

#if 0
  g_clear_object(&prt->drag_gesture);
  g_clear_object(&prt->multipress_gesture);
#endif

  free(prt);
}

static void _expose_print_page(dt_view_t *self,
                               cairo_t *cr,
                               const int32_t width,
                               const int32_t height,
                               const int32_t pointerx,
                               const int32_t pointery)
{
  dt_print_t *prt = (dt_print_t *)self->data;

  if(prt->pinfo == NULL)
    return;

  float px=.0f, py=.0f, pwidth=.0f, pheight=.0f;
  float ax=.0f, ay=.0f, awidth=.0f, aheight=.0f;

  gboolean borderless = FALSE;

  dt_get_print_layout(prt->pinfo, width, height,
                      &px, &py, &pwidth, &pheight,
                      &ax, &ay, &awidth, &aheight, &borderless);

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

#if 0
  // x page -> x display
  // (x / pg_width) * p_width + p_x
  cairo_set_source_rgb (cr, 0.9, 0.9, 0.9);
  cairo_rectangle (cr, px, py, pwidth, pheight);
  cairo_fill (cr);
  //printf("print expose paper area %f x %f + %f x %f\n", px, py, pwidth, pheight);
#endif

  // record the screen page dimension. this will be used to compute the actual
  // layout of the areas placed over the page.
  // FIXME: This is a bit convoluted: in view expose we figure out page dimensions in pixels, then use that here to take relative layout box dimensions and calculate their screen dimensions, which then are used in print_settings gui_post_expose() to draw them. Is there a better way?
  dt_printing_setup_display(prt->imgs,
                            px, py, pwidth, pheight,
                            ax, ay, awidth, aheight,
                            // FIXME: this is calculated then ignored
                            borderless);

if(0) {
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
  // FIXME: negative margins are no longer allowed by the print settings UI!

  cairo_rectangle (cr, np1x, np1y, np2x-np1x, np2y-np1y);
  cairo_clip(cr);

  cairo_set_source_rgb (cr, 0.77, 0.77, 0.77);
  //printf("print expose printable area %f x %f + %f x %f\n", ax, ay, awidth, aheight);
  cairo_rectangle (cr, ax, ay, awidth, aheight);
  cairo_fill (cr);
}
}

void expose(dt_view_t *self,
            cairo_t *cri,
            int32_t width_i,
            int32_t height_i,
            int32_t pointerx,
            int32_t pointery)
{
#if 0
  // clear the current surface
  dt_gui_gtk_set_source_rgb(cri, DT_GUI_COLOR_PRINT_BG);
  cairo_paint(cri);
#endif

  // print page & borders only. Images are displayed in
  // in overlaid widgets

  // FIXME: this should be a layout widget which draws with these dimensions
  _expose_print_page(self, cri, width_i, height_i, pointerx, pointery);
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

  // we need to setup the selected image
  prt->imgs->imgid_to_load = imgid;

  return FALSE;
}

void enter(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  /* scroll filmstrip to the first selected image */
  if(prt->imgs->imgid_to_load >= 0)
  {
    // change active image
    dt_thumbtable_set_offset_image(dt_ui_thumbtable(darktable.gui->ui),
                                   prt->imgs->box[0].imgid, TRUE);
    dt_view_active_images_reset(FALSE);
    dt_view_active_images_add(prt->imgs->imgid_to_load, TRUE);
  }

  DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals, DT_SIGNAL_DEVELOP_MIPMAP_UPDATED,
                            G_CALLBACK(_print_mipmaps_updated_signal_callback),
                            (gpointer)self);

  DT_DEBUG_CONTROL_SIGNAL_CONNECT(darktable.signals, DT_SIGNAL_VIEWMANAGER_THUMBTABLE_ACTIVATE,
                            G_CALLBACK(_view_print_filmstrip_activate_callback), self);

  prt->page = dt_layout_page_new(prt->page_dsc);
  gtk_overlay_add_overlay(GTK_OVERLAY(dt_ui_center_base(darktable.gui->ui)),
                          dt_layout_page_get_widget(prt->page));
  //gtk_overlay_reorder_overlay(GTK_OVERLAY(dt_ui_center_base(darktable.gui->ui)), prt->print_layout->w_main, -1);

#if 0
  // ensure the message widgets stay on top
  // FIXME: this is probably necessary if go straight to develop view then print view -- test
  gtk_overlay_reorder_overlay(ocda, gtk_widget_get_parent(dt_ui_log_msg(darktable.gui->ui)), -1);
  gtk_overlay_reorder_overlay(ocda, gtk_widget_get_parent(dt_ui_toast_msg(darktable.gui->ui)), -1);
#endif

#if 0
  GtkWidget *widget = dt_ui_center(darktable.gui->ui);
  gtk_widget_grab_focus(widget);

  // FIXME: this should be done in layout_page for each image box created -- will then always be routed the right widget!
  gtk_drag_dest_set(widget, GTK_DEST_DEFAULT_ALL,
                    target_list_all, n_targets_all, GDK_ACTION_MOVE);
  g_signal_connect(widget, "drag-data-received", G_CALLBACK(_drag_and_drop_received), self);
  g_signal_connect(widget, "drag-motion", G_CALLBACK(_drag_motion_received), self);
#endif

  dt_control_set_mouse_over_id(prt->imgs->imgid_to_load);
}

void leave(dt_view_t *self)
{
  dt_print_t *prt = (dt_print_t*)self->data;

  /* disconnect from mipmap updated signal */
  DT_DEBUG_CONTROL_SIGNAL_DISCONNECT(darktable.signals,
                                     G_CALLBACK(_print_mipmaps_updated_signal_callback),
                                     (gpointer)self);

  /* disconnect from filmstrip image activate */
  DT_DEBUG_CONTROL_SIGNAL_DISCONNECT(darktable.signals,
                                     G_CALLBACK(_view_print_filmstrip_activate_callback),
                                     (gpointer)self);

  // FIXME: should hide prt->print_layout->w_main before destroying it?
  // this should also remove the print layout widget from the center view overlay
  dt_layout_page_destroy(prt->page);

  dt_printing_clear_boxes(prt->imgs);
//  g_signal_disconnect(widget, "drag-data-received", G_CALLBACK(_drag_and_drop_received));
//  g_signal_disconnect(widget, "drag-motion", G_CALLBACK(_drag_motion_received));
}

// clang-format off
// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.py
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-spaces modified;
// clang-format on
