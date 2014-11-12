// This file is part of darktable
// Copyright (c) 2010 Tobias Ellinghaus <houz@gmx.de>.

// darktable is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// darktable is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with darktable.  If not, see <http://www.gnu.org/licenses/>.

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "pwstorage.h"

#ifdef HAVE_LIBSECRET
#include "backend_libsecret.h"
#endif

#ifdef HAVE_GKEYRING
#include "backend_gkeyring.h"
#endif

#ifdef HAVE_KWALLET
#include "backend_kwallet.h"
#endif

#include "control/conf.h"

#include <glib.h>
#include <string.h>

/** Initializes a new pwstorage context. */
const dt_pwstorage_t* dt_pwstorage_new()
{
  /* add password storage capabilities */
#ifdef HAVE_LIBSECRET
  dt_capabilities_add("libsecret");
#endif
#ifdef HAVE_KWALLET
  dt_capabilities_add("kwallet");
#endif
#ifdef HAVE_GKEYRING
  dt_capabilities_add("gnome-keyring");
#endif

  dt_pwstorage_t *pwstorage = g_malloc(sizeof(dt_pwstorage_t));
  dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] Creating new context %p\n", pwstorage);

  if(pwstorage == NULL)
    return NULL;

  gchar* _backend_str = dt_conf_get_string( "plugins/pwstorage/pwstorage_backend" );
  gint _backend = PW_STORAGE_BACKEND_NONE;

  if (strcmp(_backend_str, "auto") == 0)
  {
    const gchar *desktop = getenv("XDG_CURRENT_DESKTOP");
    if (g_strcmp0(desktop, "KDE") == 0)
      _backend = PW_STORAGE_BACKEND_KWALLET;
    else if (g_strcmp0(desktop, "GNOME") == 0)
      _backend = PW_STORAGE_BACKEND_GNOME_KEYRING;
    else if (g_strcmp0(desktop, "Unity") == 0)
      _backend = PW_STORAGE_BACKEND_GNOME_KEYRING;
    else if (g_strcmp0(desktop, "XFCE") == 0)
      _backend = PW_STORAGE_BACKEND_GNOME_KEYRING;

    dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] autodetected storage backend.\n");
  }
  else if(strcmp(_backend_str, "none") == 0)
    _backend = PW_STORAGE_BACKEND_NONE;
#ifdef HAVE_LIBSECRET
  else if(strcmp(_backend_str, "libsecret") == 0)
    _backend = PW_STORAGE_BACKEND_LIBSECRET;
#endif
#ifdef HAVE_KWALLET
  else if(strcmp(_backend_str, "kwallet") == 0)
    _backend = PW_STORAGE_BACKEND_KWALLET;
#endif
#ifdef HAVE_GKEYRING
  else if(strcmp(_backend_str, "gnome keyring") == 0)
    _backend = PW_STORAGE_BACKEND_GNOME_KEYRING;
#endif

  g_free(_backend_str);

  switch(_backend)
  {
    default:
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] unknown storage backend. Using none.\n");
    case PW_STORAGE_BACKEND_NONE:
      pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
      pwstorage->backend_context = NULL;
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] no storage backend. not storing username/password. please change in preferences, core tab.\n");
      break;
    case PW_STORAGE_BACKEND_LIBSECRET:
#ifdef HAVE_LIBSECRET
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] using libsecret backend for username/password storage");
      pwstorage->backend_context = (void*)dt_pwstorage_libsecret_new();
      if (pwstorage->backend_context == NULL)
      {
        dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] error starting libsecret. using no storage backend.\n");
        pwstorage->backend_context = NULL;
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
      } else {
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_LIBSECRET;
      }
      break;
#else
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] libsecret backend not available. using no storage backend.\n");
      pwstorage->backend_context = NULL;
      pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
#endif
    case PW_STORAGE_BACKEND_KWALLET:
#ifdef HAVE_KWALLET
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] using kwallet backend for username/password storage");
      pwstorage->backend_context = (void*)dt_pwstorage_kwallet_new();
      if(pwstorage->backend_context == NULL)
      {
        dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] error starting kwallet. using no storage backend.\n");
        pwstorage->backend_context = NULL;
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
      }
      else
      {
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_KWALLET;
      }
      dt_print(DT_DEBUG_PWSTORAGE,"  done.\n");
      break;
#else
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] kwallet backend not available. using no storage backend.\n");
      pwstorage->backend_context = NULL;
      pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
#endif
    case PW_STORAGE_BACKEND_GNOME_KEYRING:
#ifdef HAVE_GKEYRING
      dt_print (DT_DEBUG_PWSTORAGE,"[pwstorage_new] using gnome keyring backend for usersname/password storage.\n");
      pwstorage->backend_context = (void*)dt_pwstorage_gkeyring_new ();
      if (pwstorage->backend_context == NULL)
      {
        dt_print (DT_DEBUG_PWSTORAGE,"[pwstorage_new] error starting gnome keyring. using no storage backend.\n");
        pwstorage->backend_context = NULL;
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
      }
      else
        pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_GNOME_KEYRING;
#else
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] gnome keyring storage not available. using no storage backend.\n");
      pwstorage->backend_context = NULL;
      pwstorage->pw_storage_backend = PW_STORAGE_BACKEND_NONE;
#endif
      break;
  }

  switch(pwstorage->pw_storage_backend)
  {
    case PW_STORAGE_BACKEND_NONE:
      dt_conf_set_string( "plugins/pwstorage/pwstorage_backend", "none" );
      break;
    case PW_STORAGE_BACKEND_LIBSECRET:
      dt_conf_set_string( "plugins/pwstorage/pwstorage_backend", "libsecret" );
      break;
    case PW_STORAGE_BACKEND_KWALLET:
      dt_conf_set_string( "plugins/pwstorage/pwstorage_backend", "kwallet" );
      break;
    case PW_STORAGE_BACKEND_GNOME_KEYRING:
      dt_conf_set_string( "plugins/pwstorage/pwstorage_backend", "gnome keyring" );
      break;
  }

  return pwstorage;
}

/** Cleanup and destroy pwstorage context. \remarks After this point pointer at pwstorage is invalid. */
void dt_pwstorage_destroy(const dt_pwstorage_t *pwstorage)
{
  dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_new] Destroying context %p\n", pwstorage);
  switch(darktable.pwstorage->pw_storage_backend)
  {
    case PW_STORAGE_BACKEND_NONE:
      // nothing to be done
      break;
    case PW_STORAGE_BACKEND_LIBSECRET:
#ifdef HAVE_LIBSECRET
      dt_pwstorage_libsecret_destroy(pwstorage->backend_context);
#endif
      break;
    case PW_STORAGE_BACKEND_KWALLET:
#ifdef HAVE_KWALLET
      dt_pwstorage_kwallet_destroy(pwstorage->backend_context);
#endif
      break;
    case PW_STORAGE_BACKEND_GNOME_KEYRING:
#ifdef HAVE_GKEYRING
      g_free(pwstorage->backend_context);
#endif
      break;
  }
}

/** Store (key,value) pairs. */
gboolean dt_pwstorage_set(const gchar* slot, GHashTable* table)
{
  switch(darktable.pwstorage->pw_storage_backend)
  {
    case PW_STORAGE_BACKEND_NONE:
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_set] no backend. not storing anything.\n");
      break;
    case PW_STORAGE_BACKEND_LIBSECRET:
#if HAVE_LIBSECRET
      return dt_pwstorage_libsecret_set((backend_libsecret_context_t*)
          darktable.pwstorage->backend_context, slot, table);
#endif
      break;
    case PW_STORAGE_BACKEND_KWALLET:
#ifdef HAVE_KWALLET
      return dt_pwstorage_kwallet_set((backend_kwallet_context_t*)darktable.pwstorage->backend_context, slot, table);
#endif
      break;
    case PW_STORAGE_BACKEND_GNOME_KEYRING:
#ifdef HAVE_GKEYRING
      return dt_pwstorage_gkeyring_set(slot, table);
#endif
      break;
  }
  return FALSE;
}

/** Load (key,value) pairs. */
GHashTable* dt_pwstorage_get(const gchar* slot)
{
  switch(darktable.pwstorage->pw_storage_backend)
  {
    case PW_STORAGE_BACKEND_NONE:
      dt_print(DT_DEBUG_PWSTORAGE,"[pwstorage_get] no backend. not reading anything.\n");
      break;
    case PW_STORAGE_BACKEND_LIBSECRET:
#if HAVE_LIBSECRET
      return dt_pwstorage_libsecret_get((backend_libsecret_context_t*)
          darktable.pwstorage->backend_context, slot);
#endif
      break;
    case PW_STORAGE_BACKEND_KWALLET:
#ifdef HAVE_KWALLET
      return dt_pwstorage_kwallet_get((backend_kwallet_context_t*)darktable.pwstorage->backend_context, slot);
#endif
      break;
    case PW_STORAGE_BACKEND_GNOME_KEYRING:
#ifdef HAVE_GKEYRING
      return dt_pwstorage_gkeyring_get(slot);
#endif
      break;
  }

  return g_hash_table_new(g_str_hash, g_str_equal);
}

// modelines: These editor modelines have been set for all relevant files by tools/update_modelines.sh
// vim: shiftwidth=2 expandtab tabstop=2 cindent
// kate: tab-indents: off; indent-width 2; replace-tabs on; indent-mode cstyle; remove-trailing-space on;
