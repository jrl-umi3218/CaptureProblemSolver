#pragma once
/* Copyright 2018 CNRS-AIST JRL, CNRS-UM LIRMM
 *
 * This file is part of CPS.
 *
 * CPS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CPS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with CPS.  If not, see <http://www.gnu.org/licenses/>.
 */

# if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define CPS_DLLIMPORT __declspec(dllimport)
#  define CPS_DLLEXPORT __declspec(dllexport)
#  define CPS_DLLLOCAL
# else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#   define CPS_DLLIMPORT __attribute__ ((visibility("default")))
#   define CPS_DLLEXPORT __attribute__ ((visibility("default")))
#   define CPS_DLLLOCAL  __attribute__ ((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#   define CPS_DLLIMPORT
#   define CPS_DLLEXPORT
#   define CPS_DLLLOCAL
#  endif // __GNUC__ >= 4
# endif // defined _WIN32 || defined __CYGWIN__

# ifdef CPS_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define CPS_DLLAPI
#  define CPS_LOCAL
# else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef CPS_EXPORTS
#   define CPS_DLLAPI CPS_DLLEXPORT
#  else
#   define CPS_DLLAPI CPS_DLLIMPORT
#  endif // CPS_EXPORTS
#  define CPS_LOCAL CPS_DLLLOCAL
# endif // CPS_STATIC
