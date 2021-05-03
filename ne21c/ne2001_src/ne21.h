#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <stdlib.h>
#include <map>
#include <string>

#ifndef ne21_h__
#define ne21_h__
#ifdef __cplusplus
 extern "C" {
#endif
#include <stdio.h>
#include "f2c.h"

int  MAIN__( ) {  return 0; }

int ndir;
real r1, r2, dm, dist, sm, smtau, smtheta, smiso;
static char limit[1];

extern /* Subroutine */ int dmdsm_(real *, real *, integer *, real *,
        real *, char *, real *, real *, real *, real *, ftnlen);

 
#endif  // ne21_h__
#ifdef __cplusplus
}
#endif /* extern "C" */
