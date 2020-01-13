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

#ifndef _AABB_H_
#define _AABB_H_

#include "bunnykiller.h"

#include <limits>

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Ray.h"

/** Axis Aligned Bounding Box.
	Templatized to support multiple dimensions. */
template<unsigned D>
class AABB
{
public:
	VectorN<D> _min, _max;
public:
	AABB<D>() :
		_min(std::numeric_limits<Real>::infinity()),
		_max(-std::numeric_limits<Real>::infinity())
	{}

	AABB<D>(const VectorN<D> &_minpoint, const VectorN<D> &_maxpoint) :
		_min(_minpoint),
		_max(_maxpoint)
	{}

	AABB<D>(const AABB<D> &bb):
		_min(bb._min),
		_max(bb._max)
	{}

	void add(const VectorN<D> &p)
	{
		_min = minimum(p, _min);
		_max = maximum(p, _max);
	}

	void add(const AABB<D> &bb)
	{
		_min = minimum(bb._min, _min);
		_max = maximum(bb._max, _max);
	}

	bool intersect(const Ray<D>& r, Intersection<D>& it, float &max_distance) const
	{
		return intersect(r, max_distance);
	}

	bool intersect(const Ray<D>& r, Real &t)const
	{
		if (_min[0] > _max[0])
			return false;

		VectorN<D> origin = r.get_origin();
		Real t0 = -std::numeric_limits<Real>::infinity(),
			 t1 = std::numeric_limits<Real>::infinity();

		/*
		* Ray-box intersection from PBRT
		*/
		for (size_t i = 0; i < D; ++i) {
			// Update interval for _i_th bounding box slab
			Real invRayDir = 1.f / r.get_direction()[i];
			Real tNear = (_min[i] - origin[i]) * invRayDir;
			Real tFar  = (_max[i] - origin[i]) * invRayDir;

			// Update parametric interval from slab intersection $t$s
			if (tNear > tFar)
				std::swap(tNear, tFar);
			t0 = tNear > t0 ? tNear : t0;
			t1 = tFar  < t1 ? tFar  : t1;
			if (t0 > t1)
				return false;
		}
		t = t0;
		return true;
	}
	
	bool intersect(const VectorN<D>& p)const
	{
		if (_min[0] > _max[0])
			return false;
		/*
		 * Point-box "intersection"
		 */
		for (size_t i = 0; i < D; ++i) {
			if (p[i] < _min[i] || p[i] > _max[i])
				return false;
		}
		return true;
	}

	const AABB<D> &operator=(const AABB<D>& bb)
	{
		_min = bb._min;
		_max = bb._max;
		
		return *this;
	}

	inline int get_max_axis() const
	{
		VectorN<D> size = _max-_min;
		Real max_dim = size[0];
		int dim = 0;
		
		for (size_t i=0; i<D; ++i) {
			if (max_dim < size[i]) {
				dim = i;
				max_dim = size[i];
			}
		}

		return dim;
	}

	inline VectorN<D> get_center() const
	{
		if (_min[0] > _max[0])
			return VectorN<D>(0.);

		VectorN<D> size = _max+_min;

		return size*.5f;
	}

	inline VectorN<D> dimensions() const
	{
		if (_min[0] > _max[0])
			return VectorN<D>(std::numeric_limits<Real>::infinity());

		return _max-_min;
	}
	
	inline Real get_max_dimension() const
	{
		int max_axis = get_max_axis();
		return _max[max_axis]-_min[max_axis];
	}

	Real get_intersection_cost()const
	{
		return 1.;
	}
	
	AABB<D> get_bounding_box() const
	{
		return *this;
	}

	Real get_diagonal() const
	{
		return dimensions().length();
	}

	Real area() const;
}; // AABB

// Specializations
/**	To be done... Oh well... */
template<unsigned D>
Real AABB<D>::area() const
{ 
	return 0;
}

/**	This function does not compute the actual area
	of a AABB in 2D, but its perimeter. This is because
	this function is mainly thought to serve as the base of
	the SAH used for building the acceleration structures*/
template<>
Real AABB<2>::area()const
{
	if (_min[0] > _max[0])
		return std::numeric_limits<Real>::infinity();

	VectorN<2> dim = dimensions();
	return 2*(dim[0] + dim[1]);
}

template<>
Real AABB<3>::area()const
{
	if (_min[0] > _max[0])
		return std::numeric_limits<Real>::infinity();

	VectorN<3> dim = dimensions();
	return 2*(dim[0]*dim[1] + dim[1]*dim[2] + dim[0]*dim[2]);
}

typedef AABB<3> AABB3D;
typedef AABB<2> AABB2D;

#endif //_AABB_H_
