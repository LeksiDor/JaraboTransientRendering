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

#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "Geometry/Object.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"


/** The notion of a sphere should be familiar. A sphere is a solid, hence
		may transmit light. */
class Sphere : public Object3D
{
	Vector3 position;
	Real radius;
  
	Vector3 get_normal(const Vector3& surf_point) const
	{
		Vector3 norm = surf_point - position;
		norm.normalize();
		return norm;
	}

public:
	Sphere(const Vector3& pos, const Real rad, Material3D* _mat)
		: Object3D(_mat), position(pos), radius(rad) 
	{}

	
	
	void transform(const Matrix4x4&)
	{
	   //Vector3 radius_vec = m.mul_3D_point(Vector3(radius,0,0)+position);
	   //position = m.mul_3D_point(position); 
	   // The radius is scaled by the X scaling factor.
	   // Not ideal, but the best we can do without elipsoids
	   //radius_vec -= position;
	   //radius = radius_vec.length();
	}
  
	virtual bool intersect(Ray3D& r, Intersection3D &it, float max_distance) const
	{
		Vector3 sphere2ray_pos = r.get_origin() - position;

		Real A = 1.;
		Real B = 2.*dot(sphere2ray_pos, r.get_direction());
		Real C = dot(sphere2ray_pos, sphere2ray_pos) - radius*radius;

		Real Disc = B*B - 4.*A*C;
		if (Disc >= 0.) {
			Real t;

			if (Disc == 0.) {
				t = - B / 2.*A;
				if (t < max_distance && r.cond_set_parameter(t)) {
					it.set(r, this, get_normal(r.get_position()), Vector2());
					return true;
				}
			} else {
				// Calculate both intersections
				Real Disc_sqrt = sqrt(Disc);
				Real t1 = (-B + Disc_sqrt) / 2.*A;
				Real t2 = (-B - Disc_sqrt) / 2.*A;

				bool b1 = (t1 < max_distance) && r.cond_set_parameter(t1);
				bool b2 = (t2 < max_distance) && r.cond_set_parameter(t2);
				if (b1 || b2) {
					it.set(r, this, get_normal(r.get_position()), Vector2());
					return true;
				}
			}
		}
		return false;
	}

	virtual bool intersect(const Ray3D& r, float max_distance) const
	{
		Vector3 sphere2ray_pos = r.get_origin() - position;

		Real A = 1;
		Real B = 2 * dot(sphere2ray_pos, r.get_direction());
		Real C = dot(sphere2ray_pos, sphere2ray_pos)
			- radius * radius;

		Real Disc = B*B-4*A*C;
		if (Disc >= 0)
		{
			Real t;

			if (Disc == 0)
			{
				t = - B / 2 * A;
				if( t<max_distance && t > 0. )
					return true;
			}
			else
			{
				// Calculate both intersections
				Real Disc_sqrt = sqrt(Disc);
				Real t1 = (- B + Disc_sqrt) / 2 * A;
				Real t2 = (- B - Disc_sqrt) / 2 * A;

				if( (t1<max_distance && t1 > 0.) || (t2<max_distance && t2 > 0.) )
					return true;
			}
		}
		return false;
	}

	
	AABB3D get_bounding_box()const
	{
		return AABB3D(Vector3(position[0]-radius, position[1]-radius, position[2]-radius),
			Vector3(position[0]+radius, position[1]+radius, position[2]+radius));
	}

	Vector3 get_center()const
	{
		return position;
	}
	Real get_intersection_cost()const
	{
		return 5577.;
	}
}; // Sphere
#endif //__SPHERE_H__
