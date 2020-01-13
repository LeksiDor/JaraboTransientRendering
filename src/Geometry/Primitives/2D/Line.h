/* 
	Created by:			Julio Marco 
	Date of creation:	March 2013
 */

#ifndef __LINE_H
#define __LINE_H

#include "LinearAlgebra/Vector2.h"
#include "Geometry/Object.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"

/** 2D Line. */
class Line : public Object2D
{
	Vector2 v0, v1;
	Vector2 normal;
	AABB2D bb;

	/* Intersection data 
		v: v0;
		u: v1-v0; (not normalized)
		d: r.direction;
		o: r.origin;
		Intersection point: v + h*u; (0 <= h <= 1)
	*/
	Real vx, vy, ux, uy;

	void setup()
	{
		vx = v0[0];
		vy = v0[1];
		ux = v1[0] - v0[0];
		uy = v1[1] - v0[1];
		
		bb = AABB2D(Vector2(std::min(v0[0], v1[0]),std::min(v0[1], v1[1])), 
					Vector2(std::max(v0[0], v1[0]),std::max(v0[1], v1[1])));
	}

public:

	Vector2 get_normal(const Vector2& surf_point) const { return normal; }

	Line(const Vector2& _v0, const Vector2& _v1,
						Material<2>* _mat)
		: Object2D(_mat), v0(_v0), v1(_v1),
			normal(normalize(Vector2(v0[1]-v1[1], v1[0]-v0[0])))
	{
		setup();
	}

	void transform(const Matrix3x3& m)
	{
		// TODO
	}

	bool intersect(Ray2D& r, Intersection2D &it, float max_distance) const
	{
		Real dx, dy, ox, oy, h, den, num, t;
		Real dn = dot(r.get_direction(),normal); 
		
		if ( fabs(dn) < 1e-6) return false;
		
		ox = r.get_origin()[0];
		oy = r.get_origin()[1];
		dx = r.get_direction()[0];
		dy = r.get_direction()[1];
		
		den = ux*dy - uy*dx;
		if (fabs(den) < 1e-6) 
		{
			return false;
		}
		
		num = (ox-vx)*dy + (vy-oy)*dx;
		h = num/den;
		
		if (h < 0.0f || h > 1.0f) return false;
		
		if (fabs(dx) > fabs(dy))
			t = (vx+ux*h-ox)/dx;
		else
			t = (vy+uy*h-oy)/dy;
		
		if(r.cond_set_parameter(t))
		{	
			it.set(r, this, normal, Vector2(0.,0.));
			return true;
		}
		return false;
	}
	
	bool intersect(const Ray2D& r, float max_distance) const
	{
		Real dx, dy, ox, oy, h, den, num, t;
		Real dn = dot(r.get_direction(),normal); 	
		
		if ( fabs(dn) < 1e-6) return false;
		
		ox = r.get_origin()[0];
		oy = r.get_origin()[1];
		dx = r.get_direction()[0];
		dy = r.get_direction()[1];
		
		den = ux*dy - uy*dx;
		if (fabs(den) < 1e-6) return false;
		
		num = (ox-vx)*dy + (vy-oy)*dx;
		h = num/den;
		
		if (h < 0.0f || h > 1.0f) return false;
		
		if (fabs(dx) > fabs(dy))
			t = (vx+ux*h-ox)/dx;
		else
			t = (vy+uy*h-oy)/dy;
		
		if (t < SMALLEST_DIST || t > max_distance || r.get_parameter() < t) 
			return false;		
		
		return true;
	}
	
	AABB2D get_bounding_box() const
	{
		return bb;
	}
	
	Vector2 get_center() const
	{
		return 0.5f*(v1+v0);
	}
	
	Real get_intersection_cost() const
	{
		return 0;
	}
	
	virtual void get_lambertiancoords(const VectorN<1>& c, Vector2& pos, Vector2& n)
	{ //c range is [0,1]
		Vector2 u = Vector2(ux, uy);
		pos = v1 - c[0]*u;
		n = get_normal(pos);
		//printf("lambertian pos: %f %f; normal: %f %f\n", pos[0], pos[1], n[0], n[1]);
	}
};
#endif
