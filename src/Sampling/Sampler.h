/* 
 * Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom
 * the Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
 * OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _SAMPLER_H_
#define _SAMPLER_H_

#include "bunnykiller.h"

#include "Sample.h"
#include "Utils/RandomNumbers.h"

class Sampler
{
protected:
	bool single_pixel, single_scanline;
	int single_pixel_x, single_pixel_y;
	int size_x, size_y;
	int spp;
	Real inv_x, inv_y, inv_spp;
public:
	Sampler(int _size_x, int _size_y, int _spp) :
		single_pixel(false), single_scanline(false),
		single_pixel_x(-1), single_pixel_y(-1),
		size_x(_size_x), size_y(_size_y), spp(_spp),
		inv_x(1./Real(_size_x)), inv_y(1./Real(_size_y)),
		inv_spp(1./Real(_spp))
	{}
	
	virtual ~Sampler() {}

	virtual bool get_next_sample(Sample &sample) = 0;

	virtual size_t get_nb_samples() const
	{
		return size_t(size_x)*size_t(size_y)*size_t(spp);
	}

	virtual void restart() = 0;

	virtual void set_single_pixel(int x, int y)
	{
		single_pixel = true;
		single_pixel_x = x;
		single_pixel_y = y;
	}

	virtual void set_scanline(int y)
	{
		single_scanline = true;
		single_pixel_y = y;
	}
};

#endif //_SAMPLER_H_
