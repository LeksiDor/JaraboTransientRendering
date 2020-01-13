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

#ifndef _OBJECT_H_
#define _OBJECT_H_

#include <vector>

#include "bunnykiller.h"
#include "LinearAlgebra/MatrixN.h"
#include "Material/Material.h"
#include "RayTracing/Ray.h"
#include "RayTracing/AABB.h"

template<unsigned int D> class Material;
template<unsigned int D> class Intersection;
template<unsigned int D> class Ray;

/** This is the abstract ancestor of all traceable objects, both in 2D and 3D.
	Templatized to support multiple dimensions.		*/
template<unsigned int D>
class Object
{
protected:
	// An object contains a surface
	Material<D>* mat;
public:
	Object(Material<D>* _mat) :
		mat(_mat)
	{}

	virtual ~Object() {}

    /** Transform an object. Objects are represented directly in world
		coordinates. This function is used to transform the objects 
		according to an arbitrary transfer matrix. */
	virtual void transform(const MatrixN<D+1>& m) = 0;

    /** Returns the first intersection of a ray with the geometry.	*/
	bool intersect(Ray<D>& r, Intersection<D> &it) const
	{
		return intersect(r,it,std::numeric_limits<Real>::infinity());
	}
	virtual bool intersect(Ray<D>& r, Intersection<D> &it, float max_distance) const = 0;
    
	/** Returns if there is any intersection between ray and geometry.	*/ 
	bool intersect(const Ray<D>& r) const
	{
		return intersect(r, std::numeric_limits<Real>::infinity());
	}

	virtual bool intersect(const Ray<D>& r, float max_distance) const = 0;

	// Returns the material bounded to the object
	virtual const Material<D>* material() const
	{
		return mat;
	}
  
	// Gives the Bounding Box of the object
	virtual AABB<D> get_bounding_box() const = 0;

	// Gives the center of the object
	virtual VectorN<D> get_center() const = 0;
	
	// Gives primitives, if exists
	virtual void push_primitives(std::vector<Object<D>*> &objects)
	{
		objects.push_back(this);
	}

	// Gives the cost of intersecting the primitive; used for hierarchies construction
	// The values for each primitive has been computed separately using a simple program
	// that tests a ray-primitive intersection for several millions of rays. Then, they 
	// have been normalized so the one with smallest cost (AABB) has cost 1.
	virtual Real get_intersection_cost() const = 0;

	// Gets the x,y coordinates for the lambertian camera
	virtual void get_lambertiancoords(const VectorN<D-1>& c, VectorN<D>& pos, VectorN<D>& n)
	{}

}; //Object

typedef Object<3> Object3D;
typedef Object<2> Object2D;

#endif //__OBJECT_H_
