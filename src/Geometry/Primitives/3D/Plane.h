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

#ifndef _PLANE_H_
#define _PLANE_H_

#include "LinearAlgebra/Vector3.h"
#include "Geometry/Object.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include <math.h>

/** Infinite plane object. Not a solid. */
class Plane : public Object3D
{
	Vector3 pos;
	Vector3 normal;
	Real d;
public:
	Plane(const Vector3& _pos, const Vector3& _normal, Material3D* _mat):

	Object3D(_mat), pos(_pos), normal(normalize(_normal)), d(-dot(pos,normal))
	{}

	virtual void transform(const Matrix4x4& m)
	{}

	bool intersect(Ray3D& r, Intersection3D &it, float max_distance ) const
	{
		// dir must be normalized ! 
		Real dnd;
		if (fabs(dnd = dot(normal,r.get_direction())) > 1e-6)
		{		
			Real dst = -(d + dot(r.get_origin(),normal)) / dnd;
			if( dst < max_distance && r.cond_set_parameter(dst) )
			{
				it.set(r, this, normal, Vector2());
				return true;
			}
		}
		//it = Intersection();
		return false;
	}

	bool intersect(const Ray3D& r, float max_distance ) const
	{
		// dir must be normalized ! 
		Real dnd;
		if (fabs(dnd = dot(normal,r.get_direction())) > 1e-6)
		{		
			Real dst = -(d + dot(r.get_origin(),normal)) / dnd;
			if( dst < max_distance )
				return true;
		}
		//it = Intersection();
		return false;
	}

	AABB3D get_bounding_box() const
	{
		return AABB3D(pos,pos);
	}
	
	Vector3 get_center() const
	{
		return pos;
	}

	Real get_intersection_cost()const{return .5;};
}; //Plane


#endif //_PLANE_H_
