#pragma once

# if defined _WIN32 || defined __CYGWIN__
// On Microsoft Windows, use dllimport and dllexport to tag symbols.
#  define BMS_DLLIMPORT __declspec(dllimport)
#  define BMS_DLLEXPORT __declspec(dllexport)
#  define BMS_DLLLOCAL
# else
// On Linux, for GCC >= 4, tag symbols using GCC extension.
#  if __GNUC__ >= 4
#   define BMS_DLLIMPORT __attribute__ ((visibility("default")))
#   define BMS_DLLEXPORT __attribute__ ((visibility("default")))
#   define BMS_DLLLOCAL  __attribute__ ((visibility("hidden")))
#  else
// Otherwise (GCC < 4 or another compiler is used), export everything.
#   define BMS_DLLIMPORT
#   define BMS_DLLEXPORT
#   define BMS_DLLLOCAL
#  endif // __GNUC__ >= 4
# endif // defined _WIN32 || defined __CYGWIN__

# ifdef BMS_STATIC
// If one is using the library statically, get rid of
// extra information.
#  define BMS_DLLAPI
#  define BMS_LOCAL
# else
// Depending on whether one is building or using the
// library define DLLAPI to import or export.
#  ifdef BMS_EXPORTS
#   define BMS_DLLAPI BMS_DLLEXPORT
#  else
#   define BMS_DLLAPI BMS_DLLIMPORT
#  endif // BMS_EXPORTS
#  define BMS_LOCAL BMS_DLLLOCAL
# endif // BMS_STATIC
