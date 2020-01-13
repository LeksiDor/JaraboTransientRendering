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

#ifndef _LIGHT_SAMPLE_H_
#define _LIGHT_SAMPLE_H_

#include "bunnykiller.h"

#include "LinearAlgebra/VectorN.h"
#include "Color/Spectrum.h"

/** The LightSample class contains very basic information about a light source
	that has been hit by a shadow ray - i.e. is not occluded. Basically just the 
	direction to the light source, the distance and the amount of light. 
	Templatized to support multiple dimensions.								*/
template<unsigned D, class Radiance>
struct LightSample 
{
	VectorN<D> pos;
	VectorN<D> dir;
	Real dist;
	Radiance irradiance;
	Real instant;
}; // LightSample


#endif //_LIGHT_SAMPLE_H_
