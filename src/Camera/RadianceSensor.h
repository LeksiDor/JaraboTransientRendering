/*
 *  RadianceSensor.h
 */

#ifndef __RADIANCESENSOR_H
#define __RADIANCESENSOR_H

#include "Utils/globals.h"
#include "LinearAlgebra/Vector2.h"
#include "LinearAlgebra/Vector3.h"
#include "Camera/Camera.h"
#include "RayTracing/Ray.h"
#include "Geometry/Object.h"


/** This class represents an axis aligned radiance sensor in the scene */
template<int D>
class RadianceSensor: public Camera<D>
{
protected:
	VectorN<D> center;
	Real h,w,d;
public:
	
	/** Build the camera from viewing information.
	 */
	RadianceSensor(VectorN<D> _center, Real _w, Real _h=0., Real _d=0.):center(_center), w(_w), h(_h), d(_d) {}
	/// Calculates the ray
	virtual const VectorN<D> get_ray_direction(const VectorN<D-1>& p) const
	{
		return VectorN<D>();
	}
	
	virtual Ray<D> get_ray(const VectorN<D-1>& p) const { return Ray<D>();};
	virtual VectorN<D> get_point(const VectorN<2>& p) const
	{
		//TODO: Make it generic for 2D and 3D
		return center;
	}
};
template<>
VectorN<2> RadianceSensor<2>::get_point(const Vector2& p) const
{
	return center + p*Vector2(w,h);
}
typedef RadianceSensor<2> RadianceSensor2D;

#endif

