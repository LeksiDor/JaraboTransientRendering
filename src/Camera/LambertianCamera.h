/*
 *  LambertianCamera.h
 */

#ifndef __LAMBERTIANCAMERA_H
#define __LAMBERTIANCAMERA_H

#include "Utils/globals.h"
#include "LinearAlgebra/Vector2.h"
#include "LinearAlgebra/Vector3.h"
#include "Camera/Camera.h"
#include "RayTracing/Ray.h"
#include "Geometry/Object.h"


/** This class represents a 2D Lambertian camera  */
template<int D>
class LambertianCamera: public Camera<D>
{
protected:
	///Object associated to the camera
	Object<D> *ob;
	
	/// Delta distance to the object surface
	Real delta;
	bool orthographic;
	
public:
	
	/** Build the camera from viewing information.
	 */
	LambertianCamera(Object<D>* _ob, const Real &_delta, bool _orthographic = false):ob(_ob),delta(_delta), orthographic(_orthographic) { }
		
	/// Calculates the ray
	virtual const VectorN<D> get_ray_direction(const VectorN<D-1>& p) const
	{
		VectorN<D> x, n;
		VectorN<D-1> shift(0.5f);
		ob->get_lambertiancoords(p+shift, x, n);
		
		return (orthographic)?n:-n;
	}
	
	virtual Ray<D> get_ray(const VectorN<D-1>& p) const
	{
		VectorN<D> x, n;
		VectorN<D-1> shift(0.5f);
		
		ob->get_lambertiancoords(p+shift, x, n);
		
		return Ray<D>(x + delta*n, (orthographic)?n:-n);
	}
	inline void set_orthographic(bool is_orthographic)
	{
		orthographic = is_orthographic;
	}
};

typedef LambertianCamera<2> LambertianCamera2D;

#endif
