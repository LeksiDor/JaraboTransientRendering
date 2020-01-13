/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _COLOR_H_
#define _COLOR_H_

#include "LinearAlgebra/Vector3.h"

namespace Color
{
	//CIE 1931 standard colorimetric observer
//	static int MIN_VISIBLE = 380;
//	static int MAX_VISIBLE = 780;
//	static int SAMPLES_SPECTRUM = 80;
//	static int DELTA_SAMPLES_SPECTRUM = 5;

	Real WL2XYZ[81][3] = {{0.001368,0.000039, 0.006450},
						{0.002236,0.000064, 0.010550},
						{0.004243,0.000120, 0.020050},
						{0.007650,0.000217, 0.036210},
						{0.014310,0.000396, 0.067850},
						{0.023190,0.000640, 0.110200},
						{0.043510,0.001210, 0.207400},
						{0.077630,0.002180, 0.371300},
						{0.134380,0.004000, 0.645600},
						{0.214770,0.007300, 1.039050},
						{0.283900,0.011600, 1.385600},
						{0.328500,0.016840, 1.622960},
						{0.348280,0.023000, 1.747060},
						{0.348060,0.029800, 1.782600},
						{0.336200,0.038000, 1.772110},
						{0.318700,0.048000, 1.744100},
						{0.290800,0.060000, 1.669200},
						{0.251100,0.073900, 1.528100},
						{0.195360,0.090980, 1.287640},
						{0.142100,0.112600, 1.041900},
						{0.095640,0.139020, 0.812950},
						{0.057950,0.169300, 0.616200},
						{0.032010,0.208020, 0.465180},
						{0.014700,0.258600, 0.353300},
						{0.004900,0.323000, 0.272000},
						{0.002400,0.407300, 0.212300},
						{0.009300,0.503000, 0.158200},
						{0.029100,0.608200, 0.111700},
						{0.063270,0.710000, 0.078250},
						{0.109600,0.793200, 0.057250},
						{0.165500,0.862000, 0.042160},
						{0.225750,0.914850, 0.029840},
						{0.290400,0.954000, 0.020300},
						{0.359700,0.980300, 0.013400},
						{0.433450,0.994950, 0.008750},
						{0.512050,1.000000, 0.005750},
						{0.594500,0.995000, 0.003900},
						{0.678400,0.978600, 0.002750},
						{0.762100,0.952000, 0.002100},
						{0.842500,0.915400, 0.001800},
						{0.916300,0.870000, 0.001650},
						{0.978600,0.816300, 0.001400},
						{1.026300,0.757000, 0.001100},
						{1.056700,0.694900, 0.001000},
						{1.062200,0.631000, 0.000800},
						{1.045600,0.566800, 0.000600},
						{1.002600,0.503000, 0.000340},
						{0.938400,0.441200, 0.000240},
						{0.854450,0.381000, 0.000190},
						{0.751400,0.321000, 0.000100},
						{0.642400,0.265000, 0.000050},
						{0.541900,0.217000, 0.000030},
						{0.447900,0.175000, 0.000020},
						{0.360800,0.138200, 0.000010},
						{0.283500,0.107000, 0.000000},
						{0.218700,0.081600, 0.000000},
						{0.164900,0.061000, 0.000000},
						{0.121200,0.044580, 0.000000},
						{0.087400,0.032000, 0.000000},
						{0.063600,0.023200, 0.000000},
						{0.046770,0.017000, 0.000000},
						{0.032900,0.011920, 0.000000},
						{0.022700,0.008210, 0.000000},
						{0.015840,0.005723, 0.000000},
						{0.011359,0.004102, 0.000000},
						{0.008111,0.002929, 0.000000},
						{0.005790,0.002091, 0.000000},
						{0.004109,0.001484, 0.000000},
						{0.002899,0.001047, 0.000000},
						{0.002049,0.000740, 0.000000},
						{0.001440,0.000520, 0.000000},
						{0.001000,0.000361, 0.000000},
						{0.000690,0.000249, 0.000000},
						{0.000476,0.000172, 0.000000},
						{0.000332,0.000120, 0.000000},
						{0.000235,0.000085, 0.000000},
						{0.000166,0.000060, 0.000000},
						{0.000117,0.000042, 0.000000},
						{0.000083,0.000030, 0.000000},
						{0.000059,0.000021, 0.000000},
						{0.000042,0.000015, 0.000000}};

	// The following functions transforms from sRGB to XYZ, and viceversa,
	// using D50 reference white.
	// Transformation matrices are obtained from:
	// http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
	Vector3 rgb2xyz(const Vector3 & rgb)
	{
		return Vector3(	0.4360747*rgb[0]+0.3850649*rgb[1]+0.1430804*rgb[2],
						0.2225045*rgb[0]+0.7168786*rgb[1]+0.0606169*rgb[2],
						0.0139322*rgb[0]+0.0971045*rgb[1]+0.7141733*rgb[2]);
	}

	Vector3 xyz2rgb(const Vector3 & xyz)
	{
		return Vector3(	+3.1338561*xyz[0]-1.6168667*xyz[1]-0.4906146*xyz[2],
						-0.9787684*xyz[0]+1.9161415*xyz[1]+0.0334540*xyz[2],
						+0.0719453*xyz[0]-0.2289914*xyz[1]+1.4052427*xyz[2]);
	}
};//Color

#endif //_COLOR_H
