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

#ifndef _INSTANCE_H_
#define _INSTANCE_H_

#include "Geometry/Object.h"

/** Implementation of a class to model instances of an object. */
template<int D>
class Instance: public Object<D>
{
	MatrixN<D+1>& m_transformation;
	Object<D> *m_baseobject;
public:
	Instance(const Object<D> *obj): m_baseobject(obj), m_transformation(MatrixN<D+1>::identity()), Object<D>(obj->mat){}
	Instance(const Object<D> *obj, Material* mat): m_baseobject(obj), m_transformation(MatrixN<D+1>::identity()), Object<D>(mat){}

	virtual void transform(const MatrixN<D+1>& m)
	{
		m_transformation *= m;
	}
	
	virtual bool intersect(Ray<D>& r, Intersection<D> &it, float time = 0.) const
	{
		// Transform ray to put it in instance space
		Ray r1 = r;//*m_transformation;
		
		// Intersect transformed ray with base geometry
		bool hit = m_baseobject->intersect(r1,it);

		// And transform back the ray to get the appropiate distance

		return hit;
	}
    
	/// Gives the Bounding Box of the object
	virtual AABB<D> get_bounding_box()const = 0;

	// Gives the center of the object
	virtual Vector<D> get_center()const = 0;

	// Gives primitives, if exists
	virtual void push_primitives(std::vector<Object<D>*> &objects){objects.push_back(this);}

}; //Instance

typedef Instance<2> Instance2D;
typedef Instance<3> Instance3D;


#endif //_INSTANCE_H_