//
// MATLAB Compiler: 4.11 (R2009b)
// Date: Fri Mar 25 09:46:57 2011
// Arguments: "-B" "macro_default" "-W" "cpplib:libplot_pop_2d" "-T" "link:lib"
// "plot_mop_pop_2d.m" 
//

#ifndef __libplot_pop_2d_h
#define __libplot_pop_2d_h 1

#if defined(__cplusplus) && !defined(mclmcrrt_h) && defined(__linux__)
#  pragma implementation "mclmcrrt.h"
#endif
#include "mclmcrrt.h"
#include "mclcppclass.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_libplot_pop_2d
#define PUBLIC_libplot_pop_2d_C_API __global
#else
#define PUBLIC_libplot_pop_2d_C_API /* No import statement needed. */
#endif

#define LIB_libplot_pop_2d_C_API PUBLIC_libplot_pop_2d_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_libplot_pop_2d
#define PUBLIC_libplot_pop_2d_C_API __declspec(dllexport)
#else
#define PUBLIC_libplot_pop_2d_C_API __declspec(dllimport)
#endif

#define LIB_libplot_pop_2d_C_API PUBLIC_libplot_pop_2d_C_API


#else

#define LIB_libplot_pop_2d_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_libplot_pop_2d_C_API 
#define LIB_libplot_pop_2d_C_API /* No special import/export declaration */
#endif

extern LIB_libplot_pop_2d_C_API 
bool MW_CALL_CONV libplot_pop_2dInitializeWithHandlers(
       mclOutputHandlerFcn error_handler, 
       mclOutputHandlerFcn print_handler);

extern LIB_libplot_pop_2d_C_API 
bool MW_CALL_CONV libplot_pop_2dInitialize(void);

extern LIB_libplot_pop_2d_C_API 
void MW_CALL_CONV libplot_pop_2dTerminate(void);



extern LIB_libplot_pop_2d_C_API 
void MW_CALL_CONV libplot_pop_2dPrintStackTrace(void);

extern LIB_libplot_pop_2d_C_API 
bool MW_CALL_CONV mlxPlot_mop_pop_2d(int nlhs, mxArray *plhs[], int nrhs, mxArray 
                                     *prhs[]);

extern LIB_libplot_pop_2d_C_API 
long MW_CALL_CONV libplot_pop_2dGetMcrID();


#ifdef __cplusplus
}
#endif

#ifdef __cplusplus

/* On Windows, use __declspec to control the exported API */
#if defined(_MSC_VER) || defined(__BORLANDC__)

#ifdef EXPORTING_libplot_pop_2d
#define PUBLIC_libplot_pop_2d_CPP_API __declspec(dllexport)
#else
#define PUBLIC_libplot_pop_2d_CPP_API __declspec(dllimport)
#endif

#define LIB_libplot_pop_2d_CPP_API PUBLIC_libplot_pop_2d_CPP_API

#else

#if !defined(LIB_libplot_pop_2d_CPP_API)
#if defined(LIB_libplot_pop_2d_C_API)
#define LIB_libplot_pop_2d_CPP_API LIB_libplot_pop_2d_C_API
#else
#define LIB_libplot_pop_2d_CPP_API /* empty! */ 
#endif
#endif

#endif

extern LIB_libplot_pop_2d_CPP_API void MW_CALL_CONV plot_mop_pop_2d(const mwArray& pop_y1, const mwArray& pop_y2, const mwArray& ext_y1, const mwArray& ext_y2);

#endif
#endif
