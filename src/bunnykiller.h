/* 
 * Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _GLOBALS_H_
#define _GLOBALS_H_

// Select world dimensions. I'm quite certain this is currently broken
#ifndef DIM
#define DIM 3
#endif

// Select spectrum type (Defaults to RGB = 3)
#ifndef SPECTRUM_CHANNELS
#define SPECTRUM_CHANNELS 3
#endif

// Select integrator (defaults to BDPT)
//#define _USE_PM_

// Use fast version (Y/N)?
#define _USE_EMBREE_

// Select extra effects to render
//#define _FLUORESCENCE_
//#define _POLARIZATION_

// Select floating point precision
#ifndef _DOUBLE_PRECISSION_
using Real = float;
#else
using Real = double;
#endif

// Max path length to preallocate in bidirectional methods
#ifndef MAX_BPT_VERTICES
#define MAX_BPT_VERTICES 64
#endif

// Undefine min/max macros (if present)
#ifdef min
#undef min
#endif

#ifdef max
#undef max
#endif

// We need M_PI
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#ifndef MAX_NB_OBJECTS_LEAF
#define MAX_NB_OBJECTS_LEAF 6
#endif

#ifndef _SIGMA_VISIBILITY_
#define _SIGMA_VISIBILITY_ 1e-5
#endif


#endif // _GLOBALS_H_
