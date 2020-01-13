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

#ifndef _AGGREGATE_H_
#define _AGGREGATE_H_

#include "Geometry/Object.h"
#include "RayTracing/AABB.h"
#include <vector>
#include <iterator>

/** Base class to model aggregates, included acceleration structures. */
template<class T, unsigned D>
class Aggregate: public Object<D>
{
protected:
	std::vector<T> primitives;
	AABB<D> bb;
	bool frozen;
public:
	Aggregate():Object<D>(0),frozen(false){}
	Aggregate(std::vector<T> &objects):Object<D>(0),frozen(false)
	{
		primitives = objects;
		for( unsigned int i=0; i<objects.size(); ++i)
			bb.add(primitives[i].get_bounding_box());
	}

	virtual ~Aggregate()
	{
		clear();
	}

	virtual void add_primitive( const T &p )
	{	
		primitives.push_back(p); 
		bb.add(p.get_bounding_box());
		frozen = false;
	}
	virtual void add_primitives( const std::vector<T> &ps )
	{
		for( unsigned int i=0; i<ps.size(); ++i)
		{
			bb.add(ps[i].get_bounding_box());
			primitives.push_back(ps[i]);
		}
		frozen = false;
	}
	virtual void freeze(){}

	virtual void clear(){ primitives.clear(); bb = AABB<D>(); }
	
	int nb_primitives() const{ return primitives.size(); }
	const T &get_primitive(int i)const{ return primitives[i]; }
	T &get_primitive(int i){ return primitives[i]; }
	virtual void transform(const MatrixN<D+1>& m){throw("No transformation here");}


	virtual bool intersect(Ray<D>& r, Intersection<D> &it, float max_distance ) const
	{
		Real t;
		if( bb.intersect(r,t) && t<r.get_parameter() && t<max_distance )
		{
			typename std::vector<T>::const_iterator pr;
			bool did_it = false;
			for( pr = primitives.begin(); pr != primitives.end(); pr++ )
				did_it = pr->intersect(r,it,max_distance) || did_it;

			return did_it;
		}
		return false;
	}
	virtual bool intersect(const Ray<D>& r, float max_distance ) const
	{
		Real t;
		if( bb.intersect(r,t) && t<r.get_parameter() && t<max_distance )
		{
			typename std::vector<T>::const_iterator pr;
			for( pr = primitives.begin(); pr != primitives.end(); pr++ )
				if( pr->intersect(r, max_distance) )
					return true;
		}
		return false;
	}

	
	// Information about the volume of the aggregate 
	AABB<D> get_bounding_box()const
	{
		return bb;
	}
	VectorN<D> get_center() const
	{
		return bb.get_center();
	}

	virtual Real get_intersection_cost()const{ return 0.; }
};

#endif //_AGGREGATE_H_
