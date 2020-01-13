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

#ifndef _CAMERA_H_
#define _CAMERA_H_

#include "bunnykiller.h"

#include <cstdio>

#include "Color/PolarizationFrame.h"
#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Ray.h"

/** Base virtual class from which inherit all types of camera.*/
template<unsigned D>
class Camera
{
public:
	Camera() {}

	virtual ~Camera() {}

	/** Get direction of viewing ray from image coords.	*/
	virtual const VectorN<D> get_ray_direction(const VectorN<D-1>& p) const = 0;
	
	/** Calculates the ray from the camera*/
	virtual Ray<D> get_ray(const VectorN<D-1>& p) const = 0;
	virtual Ray<D> get_ray(const VectorN<D-1>& p, Real &pdf) const
	{
		pdf = 1.; 
		return get_ray(p);
	}

	virtual Ray<D> get_ray(const VectorN<D-1>& p, Real &pdf, Real &time_delay) const
	{
		pdf = 1.; 
		time_delay = 0.;
		return get_ray(p);
	}

	virtual VectorN<D> get_point(const VectorN<2>&) const
	{
		return VectorN<D>(0.);
	};

	/** Get the local frame */
	virtual const PolarizationFrame<D> get_frame() const = 0;

	virtual void print(FILE *stream) const {}
}; //Camera

typedef Camera<3> Camera3D;
typedef Camera<2> Camera2D;

#endif //_CAMERA_H_
