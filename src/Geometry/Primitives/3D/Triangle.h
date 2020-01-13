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

#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include "bunnykiller.h"
#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"
#include "Geometry/Object.h"
#include "RayTracing/AABB.h"
#include "RayTracing/Ray.h"
#include "RayTracing/Intersection.h"
#include <math.h>

/** The triangle class represents triangles. */
class Triangle : public Object3D
{
	//Base triangle data
	Vector3 v0,v1,v2;
	Vector3 n0,n1,n2;
	Vector2 c0,c1,c2;
	bool smooth;
	Vector3 normal;

	// Intersection related data
	Vector3 nf;
	unsigned char	i0,i1,i2;
	Real			v01, v02, cu1,cv1,cu2,cv2;
	bool			use_u;
	Real			iu1,iv1,idet;
	Real 			d;

	void setup()
	{
		Vector3 duv = v1; duv-=v0;
		Vector3 duw = v2; duw-=v0;		
		nf = cross(duw, duv);
		nf.normalize();

		d = - dot(v0, nf);

		i0 = nf.dominant();
		i1 = (i0==0 ? 1 : 0);
		i2 = (i0==2 ? 1 : 2);
		cu1 = duv[i1];
		cv1 = duv[i2];
		cu2 = duw[i1];
		cv2 = duw[i2];
		v01 = v0[i1];
		v02 = v0[i2];
		Real det=(cu1*cv2)-(cu2*cv1);
		

		if (fabs(det)<=0.0) throw("Degenerated Triangle...");
		
		use_u = ( fabs(cu1) > fabs(cv1) );
		idet= 1.0/det; iu1=1.0/cu1; iv1=1.0/cv1;
	}
public:

	// Need to include texture coordinates!
	Triangle(const Vector3& _v0, const Vector3& _v1, const Vector3& _v2, Material3D* _mat) :
		Object3D(_mat), v0(_v0), v1(_v1), v2(_v2),
		smooth(false), normal(normalize(cross(v1-v0,v2-v0)))
	{
		setup();
	}
	Triangle(const Vector3& _v0, const Vector3& _n0, 
			 const Vector3& _v1, const Vector3& _n1, 
			 const Vector3& _v2, const Vector3& _n2, 
			 Material3D* _mat)
		: Object3D(_mat), v0(_v0), v1(_v1), v2(_v2), 
		  n0(_n0), n1(_n1), n2(_n2), smooth(true)
	{
		setup();
	}
	Triangle(const Vector3& _v0, const Vector3& _n0, const Vector2& _c0,
			 const Vector3& _v1, const Vector3& _n1, const Vector2& _c1, 
			 const Vector3& _v2, const Vector3& _n2, const Vector2& _c2,
			 Material3D* _mat)
		: Object3D(_mat), v0(_v0), v1(_v1), v2(_v2), 
		  n0(_n0), n1(_n1), n2(_n2), c0(_c0), c1(_c1), c2(_c2), 
		  smooth(true)
	{
		setup();
	}
	void transform(const Matrix4x4&)
	{
		// TODO
	}
		
	bool intersect(Ray3D& r, Intersection3D &it, float max_distance) const
	{
		Real dn = dot(r.get_direction(),nf); 
		if ( fabs(dn) < 1e-6) return false;

		// We obtain the distance of intersection with the plane
		Real _t = -(dot( nf, r.get_origin()) + d)/dn; 

		// If the distance is less than the current closest distance, just return...
		if ( _t < SMALLEST_DIST || _t > max_distance || r.get_parameter() < _t) return false;

		// ... otherwise, test if the ray intersects with the triangle
		// i.e. between the vertices.
		Vector3 p = r.get_origin()+_t*r.get_direction();


		Real tu0 = p[i1] - v01;
		Real tv0 = p[i2] - v02;

		// Compute intersection's baricentric coordinate 'v'
		Real v = (cu1*tv0-cv1*tu0)*idet;
		// If intersected coordinate 'v' is not in the range [0,1]
		// then the ray doesn't intersect in the triangle.
		
		if ((v<0.0)||(v>1.0)) { return false; } 

		// Compute intersection's baricentric coordinate 'u'
		Real u;
		if (use_u) u = (tu0-cu2*v)*iu1;
		else u = (tv0-cv2*v)*iv1;
		// Again, if intersection is not in the triangle, return...	
		
		if ((u<0.0)||(u+v)>1.0) { return false; }
		
		// Normalize v according to u
		//v = (u>=1.0)? 0.0: v/(1.0-u);

		// Finally, it intersects, so update ray...
		// (because we already know that the intersection distance 
		// is smaller than previous one and also bigger than 
		// SMALLEST_DIST, no need for testing again...)
		r.set_parameter(_t);
		// ... and intersection...

		if( smooth )
			it.set(r, this, n0*(1.-u-v)+n1*u+n2*v, c0*(1.-u-v)+c1*u+c2*v);
		else
			it.set(r, this, normal, c0*(1.-u-v)+c1*u+c2*v);

		// ... finally return that the intersection exists! 
		return true;
	}


	bool intersect(const Ray3D& r, float max_distance ) const
	{
		Real dn = dot(r.get_direction(),nf); 
		if ( fabs(dn) < 1e-6) return false;

		// We obtain the distance of intersection with the plane
		Real _t = -(dot( nf, r.get_origin()) + d)/dn; 

		// If the distance is less than the current closest distance, just return...
		if ( _t < SMALLEST_DIST || _t > max_distance || r.get_parameter() < _t) return false;

		// ... otherwise, test if the ray intersects with the triangle
		// i.e. between the vertices.
		Vector3 p = r.get_origin()+_t*r.get_direction();


		Real tu0 = p[i1] - v01;
		Real tv0 = p[i2] - v02;

		// Compute intersection's baricentric coordinate 'v'
		Real v = (cu1*tv0-cv1*tu0)*idet;
		// If intersected coordinate 'v' is not in the range [0,1]
		// then the ray doesn't intersect in the triangle.
		
		if ((v<0.0)||(v>1.0)) { return false; } 

		// Compute intersection's baricentric coordinate 'u'
		Real u;
		if (use_u) u = (tu0-cu2*v)*iu1;
		else u = (tv0-cv2*v)*iv1;
		// Again, if intersection is not in the triangle, return...	
		
		if ((u<0.0)||(u+v)>1.0) { return false; }

		// ... finally return that the intersection exists! 
		return true;
	}


	AABB3D get_bounding_box()const
	{
		AABB3D bb(v0, v0);
		bb.add(v1); bb.add(v2);

		return bb;
	}

	Vector3 get_center() const
	{
		return v0*.3+v1*.3+v2*.3;
	}
	Real get_intersection_cost()const
	{
		return 7292; //
	}

	void tessellate(std::list<Triangle> &subtriangles, Real max_diagonal)
	{
		if (this->get_bounding_box().get_diagonal() <= max_diagonal)
		{
			subtriangles.push_back(*this);
			return;
		}

		Real l01 = (v0 - v1).length();
		Real l02 = (v0 - v2).length();
		Real l12 = (v1 - v2).length();

		Vector3 v01 = v0*.5 + v1*.5 + n0*.125*l01 - n1*.125*l01;
		Vector3 v02 = v0*.5 + v2*.5 + n0*.125*l02 - n2*.125*l02;
		Vector3 v12 = v1*.5 + v2*.5 + n1*.125*l12 - n2*.125*l12;
		Vector3 n01 = (n1 + n0) / 2;
		Vector3 n02 = (n2 + n0) / 2;
		Vector3 n12 = (n1 + n2) / 2;
		Vector2 c01 = (c1 + c0) / 2;
		Vector2 c02 = (c2 + c0) / 2;
		Vector2 c12 = (c1 + c2) / 2;

		Triangle a(v0, n0, c0, v01, n01, c01, v02, n02, c02, this->mat);
		Triangle b(v01, n01, c01, v1, n1, c1, v12, n12, c12, this->mat);
		Triangle c(v02, n02, c02, v12, n12, c12, v2, n2, c2, this->mat);
		Triangle d(v01, n01, c01, v12, n12, c12, v02, n02, c02, this->mat);

		a.tessellate(subtriangles, max_diagonal);
		b.tessellate(subtriangles, max_diagonal);
		c.tessellate(subtriangles, max_diagonal);
		d.tessellate(subtriangles, max_diagonal);

	}
}; //Triangle

#endif //_TRIANGLE_H_
