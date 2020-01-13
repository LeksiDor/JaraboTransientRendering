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

#ifndef _LIGHTSOURCE_H_
#define _LIGHTSOURCE_H_

#include <vector>

#include "bunnykiller.h"
#include "Color/Spectrum.h"
#include "LightSource/LightSample.h"

template<unsigned int D, class Radiance> class World;
template<unsigned int D> class Intersection;

/** The abstract LightSource is the ancestor of all other light source classes.
	Templatized to support multiple dimensions.	 */
template<unsigned int D, class Radiance>
class LightSource
{
protected:
	World<D, Radiance>* world;
	Radiance intensities;
	Real time;
public:
	LightSource(World<D,Radiance>* _world, const Radiance &_intensities, Real _time = 0.)
		:world(_world), intensities(_intensities), time(_time) {}

	virtual ~LightSource() {}
		 
	/** Return the incoming direction from light to point */
	virtual VectorN<D> get_incoming_direction(const VectorN<D> &point_lighted) const = 0;

	/** Return if the light is visible from the point */
	virtual bool is_visible(const VectorN<D> &point_lighted) const = 0;
	virtual bool is_visible(const VectorN<D> &point_lighted, const VectorN<D> &direction) const
	{
		return is_visible(point_lighted);
	}

	/** Return the total light incoming the point from the light source */
	virtual Radiance get_incoming_light(const VectorN<D> &point_lighted) const = 0;

	/** Get the intensity of the light source. This is really the total emitted
		power - so one should be careful to divide by the number of samples in
		case of area lights */
	Radiance get_intensities() const {
		return intensities;
	}
	
	/** Number of samples to cast from a given point towards the light source. */
	virtual int get_no_samples() const {
		return 1;
	}

	/** Sample the light accordingly to intersection it. 
		Returns false if the light is not visible from it */
	virtual bool sample(const Intersection<D> &it, LightSample<D, Radiance> &light_sample, Real &pdf) const = 0;

	/** Sample the light source randomly, so it gives an origin to light-tracing algorithms */
	virtual void sample(LightSample<D, Radiance> &light_sample, Real &pdf) const = 0;
	
	/** Sample the light accordingly to a point in space.
		Returns false if the light is not visible from it */
	virtual bool sample(const VectorN<D> &p, LightSample<D, Radiance> &light_sample, Real &pdf) const = 0;

	virtual void print(FILE *_f_log) const {}

	virtual void setup() {}

}; //LightSource


#endif //_LIGHTSOURCE_H_


