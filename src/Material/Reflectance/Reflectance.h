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

#ifndef _REFLECTANCE_H_
#define _REFLECTANCE_H_

namespace Reflectance
{
	enum Type {
		NOTAG                      = 0,
		DELTA                      = 1,
		GLOSSY                     = 2,
		DIFFUSE                    = 4,
		TRANSMISSION               = 8,
		EMISSIVE                   = 16,
		TWO_SIDED                  = 32,
		DELTA_TRANSMISSION         = DELTA  | TRANSMISSION,
		GLOSSY_TRANSMISSION        = GLOSSY | TRANSMISSION,
		DIFFUSE_TRANSMISSION       = DIFFUSE | TRANSMISSION,
		GLOSSY_N_DIFFUSE           = GLOSSY | DIFFUSE,
		DELTA_TWO_SIDED            = DELTA | TWO_SIDED,
		GLOSSY_TWO_SIDED           = GLOSSY | TWO_SIDED,
		DIFFUSE_TWO_SIDED          = DIFFUSE | TWO_SIDED,
		GLOSSY_N_DIFFUSE_TWO_SIDED = GLOSSY_N_DIFFUSE | TWO_SIDED
	};

	Type merge(const Type type1, const Type type2)
	{
		return Type(type1 | type2);
	}
}; // Reflectance

#endif //_REFLECTANCE_H_
