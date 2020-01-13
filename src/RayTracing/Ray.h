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

#ifndef _RAY_H_
#define _RAY_H_

#include "bunnykiller.h"

#include <cfloat>

#include "LinearAlgebra/Vector3.h"

template<unsigned D> class Medium;

/** This is a constant used to dismiss intersections very close to previous
	intersections. Templatized to support multiple dimensions.*/
const Real SMALLEST_DIST = 1e-4;

/// This is the refraction index of air. 
const Real DEFAULT_REFRACTION_INDEX = 1.0;

/** The Ray class is used to find intersections between a ray and the scene,
	but it also stores information. For instance the Ray remembers intersection
	points, and the refraction index of material it is passing through. 
	The class is templatized, so it can be defined at any dimension.	
	[TODO] For optimization purposes, some "set" funtions would be interesting; 
	otherwise, many times we need to create a Ray, and then assign it. */
template<unsigned D>
class Ray
{		
	/// Origin of the ray
	VectorN<D> origin;

	/// The normalized direction of the ray
	VectorN<D> direction;

	/// The parameter -i.e. the distance we have traversed along the ray
	Real t;                             

	/// Level is the number of times the ray has been traced recursively.
	int level;

	/// has the ray hit something ?!
	bool hit;  

	/// Refraction index of the ray's media
	Real n;

	/// Current medium
	Medium<D> *medium;
public:
	/** Construct a ray. First argument is position. Second argument
		is the direction of the ray. The magnitude of the second argument
		is constructed as the step length. */
	Ray(const VectorN<D>& p, const VectorN<D>& d, bool normalize_d = true, int _level = 0,
			Real _n = DEFAULT_REFRACTION_INDEX, Medium<D> *_medium = nullptr):
		origin(p),
		direction(d),
		t(FLT_MAX),
		level(_level),
		hit(false),
		n(_n),
		medium(_medium)
	{
		if(normalize_d) direction.normalize();
	}

	Ray(const VectorN<D>& p, const VectorN<D>& d, Real _n, Medium<D> *_medium = nullptr) :
		origin(p),
		direction(d),
		t(FLT_MAX),
		level(0),
		hit(false),
		n(_n),
		medium(_medium)
	{
		direction.normalize();
	}

	Ray(const VectorN<D>& p, const VectorN<D>& d, int _level, Real _n = DEFAULT_REFRACTION_INDEX,
			Medium<D> *_medium = nullptr) :
		origin(p),
		direction(d),
		t(FLT_MAX),
		level(_level),
		hit(false),
		n(_n),
		medium(_medium)
	{
		direction.normalize();
	}

    Ray(Real _n = DEFAULT_REFRACTION_INDEX, Medium<D> *_medium = nullptr) :
    	t(FLT_MAX),
		level(0),
		hit(false),
		n(_n),
		medium(_medium)
    {}

    Ray(Real _t, const VectorN<D>& p, const VectorN<D>& d, Real _n = 1.0) :
		origin(p),
		direction(d),
		t(_t),
		level(0),
		hit(false),
		n(_n),
		medium(nullptr)
	{}

	Ray(const Ray<D>& r) :
		origin(r.origin),
		direction(r.direction),
		t(r.t),
		level(r.level),
		hit(r.hit),
		n(r.n),
		medium(r.medium)
	{}

	/// Get ray position.
	inline const VectorN<D> get_position() const
	{
		return origin+t*direction;
	}

	/// Get ray origin.
	inline const VectorN<D>& get_origin() const
	{
		return origin;
	}
  
	/// Get ray parameter.
	inline Real get_parameter() const
	{
		return t;
	}

	/// Get direction of ray.
	inline const VectorN<D>& get_direction() const
	{
		return direction;
	}

	inline void set_medium(Medium<D> *m)
	{
		medium = m;
	}

	inline void set_refraction_index(Real _n)
	{
		n=_n;
	}

	/** Conditional set parameter. Set parameter only if it corresponds
			to a closer hit point than stored parameter. Set the associated
			object. */
	inline bool cond_set_parameter(const Real par)
	{
		// Sets the parameter of the ray to par, if par is closer to the 
		// origin than any other previous intersection
		// If the ray hits, then the object hit is also.

		if (par < get_parameter() && par > SMALLEST_DIST) {
			hit = true;
			t = par;
			return true;
		}
		return false;
	}
	/** Non-Conditional set parameter. */
	inline void set_parameter(const Real par)
	{
		hit = true;
		t = par;
	}
	inline void increase_level()
	{
		++level;
	}

	inline bool did_hit() const
	{
		return hit;
	}

	inline int get_level() const
	{
		return level;
	}

	inline Real get_ior() const
	{
		return n;
	}

	inline Medium<D> *get_medium() const
	{
		return medium;
	}

	const Ray<D> &operator=(const Ray<D> &r)
	{
		origin = r.origin;
		direction = r.direction;
		t = r.t;
		level = r.level;
		hit = r.hit;
		n = r.n;
		medium = r.medium;
		return *this;
	}

	VectorN<D> operator+(const Real d)const
	{
		return origin + direction*d;
	}

	void print(FILE *f)const
	{
		fprintf(f, "Ray: \to [%f, %f, %f]; d [%f %f %f]\n",
				origin[0], origin[1], origin[2],
				direction[0], direction[1], direction[2]);
		if (hit)
			fprintf(f, "\tHit - %f\n", t);
	}
}; // Ray

typedef Ray<3> Ray3D;
typedef Ray<2> Ray2D;

#endif //_RAY_H_
