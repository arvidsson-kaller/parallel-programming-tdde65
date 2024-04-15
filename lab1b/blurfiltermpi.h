/*
  File: blurfilter.h
  Declaration of pixel structure and blurfilter function.
 */

#ifndef _BLURFILTERPMPI_H_
#define _BLURFILTERPMPI_H_

#include "pixel.h"
#include "mpi.h"

void blurfilterMPI(const int xsize, const int ysize, pixel* src, const int radius, const double *w);

#endif
