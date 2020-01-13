/**
 * Copyright (C) 2017, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _RANDOM_NUMBERS_H_
#define _RANDOM_NUMBERS_H_

#include "bunnykiller.h"

#include <cstdint>
#include <limits>

namespace Random
{
	/**
	 * Random Number Generator based on the PCG RNG [1]
	 *
	 *	[1]	M. E. O'Neill
	 *	    PCG: A Family of Simple Fast Space-Efficient Statistically Good Algorithms for Random Number Generation
	 *	    Harvey Mudd College tech report HMC-CS-2014-0905, 2014
	 *
	 * Original PGC C minimal library by Melissa O'Neill <oneill@pcg-random.org>
	 *
	 * Licensed under the Apache License, Version 2.0 (the "License");
	 * you may not use this file except in compliance with the License.
	 * You may obtain a copy of the License at
	 *
	 *     http://www.apache.org/licenses/LICENSE-2.0
	 *
	 * Unless required by applicable law or agreed to in writing, software
	 * distributed under the License is distributed on an "AS IS" BASIS,
	 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	 * See the License for the specific language governing permissions and
	 * limitations under the License.
	 *
	 * For additional information about the PCG random number generation scheme,
	 * including its license and other licensing options, visit
	 *
	 *     http://www.pcg-random.org
	 */
	class RNG
	{
	public:
		static const uint64_t default_state  = 0x853c49e6748fea9bULL;
		static const uint64_t default_stream = 0xda3e39cb94b95bdbULL;
	public:
		constexpr RNG() :
			m_state(default_state), m_stream(default_stream)
		{}

		RNG(uint64_t state, uint64_t seq = 1)
		{
			seed(state, seq);
		}

		RNG(const RNG& rng) :
			m_state(rng.m_state), m_stream(rng.m_stream)
		{}

		/* Generates a random real number on the interval [0, 1) */
        template<class R=Real>
		inline R next_real();

		/* Generates a random float on the interval [0, 1) */
		float next_float()
		{
			union {
				uint32_t i;
				float f;
			} x;
			/* Based on the method used by the SFMT [1], generating a random number on [1, 2) and
			 * subtracting 1
			 *
			 *	[1]	M. Saito, M. Matsumoto
			 *	    A PRNG specialized in double precision floating point numbers using an affine transition
			 *      Monte Carlo and Quasi-Monte Carlo Methods, 2008
			 */
			x.i = (next_int() >> 9) | 0x3f800000UL;
			return x.f - 1.0f;
		}

		/* Generate a random double on the interval [0, 1) */
		double next_double()
		{
			union {
				uint64_t i;
				double f;
			} x;
			/* Modification of the Saito & Matsumoto method to generate doubles from 32 bits sources,
			 * adapted from Nori [1]
			 *
			 *	[1]	Wenzel Jakob
			 *	    Nori: an educational ray tracer
			 *	    https://wjakob.github.io/nori/
			 */
			x.i = (uint64_t(next_int()) << 20) | 0x3ff0000000000000ULL;
			return x.f - 1.0;
		}

		/* Generates a random 32 bits integer */
		uint32_t next_int()
		{
			uint64_t oldstate = m_state;
			m_state = oldstate * 6364136223846793005ULL + m_stream;

			uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
			uint32_t rot = oldstate >> 59u;

			return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
		}

		/* Generates a (bounded) random 32 bits integer */
		uint32_t next_int(uint32_t bound)
		{
			/* To avoid bias, we need to make the range of the RNG a multiple of
			 * bound, which we do by dropping output less than a threshold.
			 * A naive scheme to calculate the threshold would be to do
			 *
			 *     uint32_t threshold = 0x100000000ull % bound;
			 *
			 * but 64-bit div/mod is slower than 32-bit div/mod (especially on
			 * 32-bit platforms).  In essence, we do
			 *
			 *     uint32_t threshold = (0x100000000ull-bound) % bound;
			 *
			 * because this version will calculate the same modulus, but the LHS
			 * value is less than 2^32.
			 */
			uint32_t threshold = -bound % bound;

			/* Uniformity guarantees that this loop will terminate.  In practice, it
			 * should usually terminate quickly; on average (assuming all bounds are
			 * equally likely), 82.25% of the time, we can expect it to require just
			 * one iteration.  In the worst case, someone passes a bound of 2^31 + 1
			 * (i.e., 2147483649), which invalidates almost 50% of the range.  In
			 * practice, bounds are typically small and only a tiny amount of the range
			 * is eliminated.
			 */
			for (;;) {
				uint32_t r = next_int();
				if (r >= threshold)
					return r % bound;
			}

			return 0;
		}

		/* Changes the internal state of the random number generator */
		void seed(uint64_t state, uint64_t seq = default_stream)
		{
			m_state = 0U;
			m_stream = (seq << 1u) | 1u;
			next_int();
			m_state += state;
			next_int();
		}

		/* Advances or rewinds the RNG by N positions */
		void advance(int64_t N)
		{
			/* The method used here is based on Brown, "Random Number Generation
			 * with Arbitrary Stride,", Transactions of the American Nuclear
			 * Society (Nov. 1994).  The algorithm is very similar to fast
			 * exponentiation.
			 *
			 * Even though delta is an unsigned integer, we can pass a
			 * signed integer to go backwards, it just goes "the long way round".
			 */
			uint64_t delta = uint64_t(N);

			uint64_t cur_mult = 6364136223846793005ULL;
			uint64_t cur_plus = m_stream;

			uint64_t acc_mult = 1;
			uint64_t acc_plus = 0;

			while (delta > 0u) {
			   if (delta & 1u) {
				  acc_mult *= cur_mult;
				  acc_plus = acc_plus*cur_mult + cur_plus;
			   }
			   cur_plus = (cur_mult + 1u)*cur_plus;
			   cur_mult *= cur_mult;
			   delta >>= 1;
			}

			m_state =  acc_mult * m_state + acc_plus;
		}
	private:
		/* RNG state.  All values are possible. */
	    uint64_t m_state;
	    /* Which RNG sequence (stream) is selected. Must *always* be odd. */
	    uint64_t m_stream;
	}; /* RNG */

	template<>
	inline float RNG::next_real<float>()
	{
		return next_float();
	}

	template<>
	inline double RNG::next_real<double>()
	{
		return next_double();
	}

	/* Default RNG */
	RNG StdRNG = RNG();
}; /* namespace Random */

#endif /* _RANDOM_NUMBERS_H_ */
