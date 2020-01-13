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

#ifndef _STRATIFIED_SAMPLER_H_
#define _STRATIFIED_SAMPLER_H_

#include "bunnykiller.h"

#include "Sampling/Sampler.h"



class StratifiedSampler: public Sampler
{
	int current_x, current_y;
    Real limit_x, limit_y;
	int current_subx, current_suby;
	bool jittered, end;
public:
	StratifiedSampler(int _size_x, int _size_y, int _sqrt_spp, bool _jittered = true) :
		Sampler(_size_x, _size_y, _sqrt_spp),
        current_x(0), current_y(0),
        limit_x(std::nextafter(Real(1), Real(0.))), limit_y(limit_x),
        current_subx(0), current_suby(0),
		jittered(_jittered), end(false)
	{}

	bool get_next_sample(Sample &sample)
	{
		if (end)
            return false;

        Real shift_x = (current_subx + Random::StdRNG.next_real())*inv_spp;
        Real shift_y = (current_suby + Random::StdRNG.next_real())*inv_spp;

        // Despite the RNG returning values in the [0, 1) range, sometimes
        // the result ends being (current + 1) due to rounding, so we clamp
        Real pos_x = std::min<Real>(current_x + shift_x, limit_x);
        Real pos_y = std::min<Real>(current_y + shift_y, limit_y);

        sample.position = Vector2(pos_x, pos_y);
		sample.weight = 1.;

		if (++current_subx == spp) {
			current_subx = 0;
			if (++current_suby == spp) {
				current_suby = 0;
				if (single_pixel) {
					end = true;
					return true;
				}
				
				if (++current_x == size_x) {
					current_x = 0;
					if (single_scanline) {
						end = true;
						return true;
					}

					if (++current_y == size_y) {
						end = true;
					}
				}
                limit_x = std::nextafter(Real(current_x + 1), Real(0.));
                limit_y = std::nextafter(Real(current_y + 1), Real(0.));
			}
		}

		return true;
	}

	size_t get_nb_samples() const
	{
		return single_pixel ?
				size_t(spp)*size_t(spp) :
				size_t(size_x)*size_t(size_y)*size_t(spp)*size_t(spp);
	}

	void restart()
	{
		end = false; 
		if (single_pixel) {
			current_x = single_pixel_x;
			current_y = single_pixel_y;
		} else if (single_scanline) {
			current_y = single_pixel_y;
			current_x = current_subx = current_suby = 0;
		} else {
			current_x = current_y = current_subx = current_suby = 0;
		}
        limit_x = std::nextafter(Real(current_x + 1), Real(0.));
        limit_y = std::nextafter(Real(current_y + 1), Real(0.));
	}

	virtual void set_scanline(int y)
	{
		Sampler::set_scanline(y);
		restart();
	}

	virtual void set_single_pixel(int x, int y)
	{
        Sampler::set_single_pixel(x, y);
		restart();
	}
}; // StratifiedSampler

#endif // _STRATIFIED_SAMPLER_H_
