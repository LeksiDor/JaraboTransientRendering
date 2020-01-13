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

#ifndef _OBJECT_POINTER_H_
#define _OBJECT_POINTER_H_

#include "Object.h"

/** Class that masks the pointer to an object. It is specially useful for 
	the templatized aggregates, that assume that the objects are not pointers,
	and thus are not called using the operator "->()", but ".()". 
	Templatized to support multiple dimensions.							*/
template<int D>
class ObjectPointer: public Object<D>
{
	Object<D> *o;
public:
	ObjectPointer(Object<D> *obj):Object<D>(0),o(obj){}
	
	inline void transform(const MatrixN<D+1>& m)
	{	o->transform(m);	}
  
	inline bool intersect( Ray<D>& r, Intersection<D> &it, float max_distance ) const
	{	return o->intersect(r, it, max_distance);	}
	inline bool intersect(const Ray<D>& r, float max_distance ) const
	{	return o->intersect(r, max_distance);	}


	inline const Material<D>* material() const
	{	return o->material(); }
  
	inline AABB<D> get_bounding_box()const
	{	return o->get_bounding_box();	}
	inline VectorN<D> get_center()const
	{	return o->get_center();	}
	
	inline void push_primitives(std::vector<Object<D>*> &objects)
	{	o->push_primitives(objects);	}

	inline Real get_intersection_cost()const
	{	return o->get_intersection_cost();	}
}; //ObjectPointer


typedef ObjectPointer<2> Object2DPointer;
typedef ObjectPointer<3> Object3DPointer;

#endif //_OBJECT_POINTER_H_