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

#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "bunnykiller.h"
#include "LinearAlgebra/VectorN.h"
#include "LinearAlgebra/Vector2.h"
#include "RayTracing/Ray.h"
#include "Geometry/Object.h"

/** The Intersection class is used to store the data of the intersection.
	Templatized to support multiple dimensions.  */
template<unsigned D>
class Intersection
{		
	/// has the ray hit something ?!
	bool hit; 

	//The ray traced when found the intersection
	Ray<D> ray;

	/// the object that has been hit
	// const Object<D>* obj;           
	const Material<D>* mat;

	//The position at the surface hit
	VectorN<D> position;

	//The normal at the surface hit
	VectorN<D> normal;
	VectorN<D> tu, tv;

	//The uv-coordinates at the surface hit
	Vector2 uv;
public:
	Intersection() :
		hit(false),
		ray(),
		mat(0)
	{}

	Intersection(const Ray<D> &r, const Object<D>* _obj,
			const VectorN<D> &_n, const Vector2 &_uv) :
		hit(true),
		ray(r),
		mat(_obj ? _obj->material() : nullptr),
		position(ray.get_position()),
		normal(_n),
		uv(_uv)
	{}

	Intersection(const Ray<D> &r, const Material<D>* _mat,
			const VectorN<D> &_n, const Vector2 &_uv) :
		hit(true),
		ray(r),
		mat(_mat),
		position(ray.get_position()),
		normal(_n),
		uv(_uv)
	{}

	/// Get intersection position.
	inline const VectorN<D> &get_position() const
	{
		return position;
	}
	
	/// Get intersection normal.
	inline const VectorN<D> &get_normal() const
	{
		return normal;
	}

	inline const VectorN<D> &get_tangent_x() const
	{
		return tu;
	}

	inline const VectorN<D> &get_tangent_y() const
	{
		return tv;
	}

	/// Get intersection uv coordinates.
	inline const Vector2 &get_uv() const
	{
		return uv;
	}

	/// Get ray
	inline const Ray<D>& get_ray() const
	{
		return ray;
	}

	/// Returns true if ray has hit an object.
	inline bool did_hit() const
	{
		return hit;
	}

	/// If ray has hit, return intersected object.
	//inline const Object<D>* intersected() const {return obj;}

	void set(const Ray<D> &r, const Object<D>* _obj, const VectorN<D> &_n, const Vector2 &_uv)
	{
		hit = r.did_hit();
		ray = r;
		mat = _obj ? _obj->material() : 0;
		position = r.get_position();
		normal = _n;
		uv = _uv;
	}

	void set(const Ray<D> &r, const Material<D>* _mat, const VectorN<D> &_n, const Vector2 &_uv)
	{
		hit = r.did_hit();
		ray = r;
		mat = _mat;
		position = r.get_position();
		normal = _n;
		uv = _uv;
	}

	/// If ray has hit, return intersected material
	const Material<D>* material() const
	{
		return mat;
	}

	/// Asignation...
	const Intersection<D> &operator=(const Intersection<D> &it)
	{
		hit = it.hit;
		ray = it.ray;
		mat = it.mat;
		position = it.position;
		normal = it.normal;
		tu = it.tu;
		tv = it.tv;
		uv = it.uv;

		return *this;
	}

	void set_coordinate_system()
	{
		if (hit) {
			/* Build an orthonormal base by Hughes-Moeller method */
			VectorN<D> aux = (std::abs(normal[0]) > std::abs(normal[2])) ?
					VectorN<D>(-normal[1], normal[0], normal[2]) :
					VectorN<D>(0.0, -normal[2], normal[1]);

			tv = aux.normalized();
			tu = cross(normal, tv);
		}
	}

	void invert_coordinate_system()
	{
		normal = -normal;
		tv = -tv;
		tu = -tu;
	}
}; //Intersection

typedef Intersection<3> Intersection3D;
typedef Intersection<2> Intersection2D;

#endif //_INTERSECTION_H_
