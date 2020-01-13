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

#ifndef _PATH_MIS_H_
#define _PATH_MIS_H_

#include "bunnykiller.h"

#include "Integrator/Path.h"
#include "LinearAlgebra/VectorN.h"

//=============================================================================
// Compute path weights using multiple importance sampling (MIS).
//  The function computes the MIS weight for a bidirectional sub-path
//  X_{s,t} using the power heuristic [1], as described in Veach's thesis
//  [2]. The bidirectional path is stored in 'PathL' (the light path) and
//  'PathE' (the eye path). The power heuristic is parametrized using the
//  parameter 'beta'.
//  
//  [1] E. Veach and L. Guibas. 1995. Optimally Combining Sampling Techniques for
//      Monte Carlo Rendering. In SIGGRAPH '95.
//      http://www-graphics.stanford.edu/papers/combine/
//   [2] E. Veach. 1997. Robust Monte Carlo Methods for Light Transport
//       Simulation. PhD Thesis (Stanford University).
//       http://graphics.stanford.edu/papers/veach_thesis/

namespace PathMIS
{
	//=========================================================================
	template<unsigned D, class Radiance, class RadianceAttenuation>
	Real connection_probability(const Path<D, Radiance, RadianceAttenuation> &path0, const Path<D, Radiance, RadianceAttenuation> &path1,
						unsigned int i0, unsigned int i1 )
	{
		VectorN<D> d = path1[i1].get_vertex_position() - path0[i0].get_vertex_position();
		d.normalize();

		Real p = 0.;
		path0[i0].compute_probability(d, p);

		return p;
	}

	//=========================================================================
	template<unsigned D, class Radiance, class RadianceAttenuation>
	Real compute_weight(const Path<D, Radiance, RadianceAttenuation> &eye, const Path<D, Radiance, RadianceAttenuation> &light,
						unsigned int t, unsigned int s, const Real beta)
	{
		Real iw = 1;

		// Compute light path
		if (s > 1) {
			// Compute \p_{s-1,t+1}/\p_{s,t}
			Real p_s_1 = connection_probability(eye, light, t, s);
			p_s_1 /= light[s].get_vertex_pdf();

			iw += pow(p_s_1,beta);

			Real p_i_0 = p_s_1;

			// Compute p_{i+1,i} for all i\in[s-1,0] //FIX COMMENT
			for (size_t i=s-1; i>0; --i) {
				Real p_i = connection_probability(light, light, i+1, i);
				p_i = p_i_0 + p_i/light[i].get_vertex_pdf();

				iw = iw + pow(p_i,beta);

				p_i_0 = p_i;
			}
		}

		// Compute eye path
		if ( t > 1 ) {
			// Compute \p_{s+1,t-1}/\p_{s,t}
			Real p_t_1 = connection_probability(light, eye, s, t);
			p_t_1 /= eye[t].get_vertex_pdf();

			iw += pow(p_t_1,beta);

			Real p_i_0 = p_t_1;
			// Compute p_{i-1,i} for all i\in[s-1,0] //FIX COMMENT
			for (size_t i=t-1; i>0; --i) {
				Real p_i = connection_probability(eye, eye,i+1, i);
				p_i = p_i_0 + p_i/eye[i].get_vertex_pdf();

				iw = iw + pow(p_i,beta);

				p_i_0 = p_i;
			}
		}

		return 1./iw;
	}
}; //PathMIS

#endif //_PATH_MIS_H_
