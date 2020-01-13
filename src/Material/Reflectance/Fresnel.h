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

#ifndef _FRESNEL_H_
#define _FRESNEL_H_

#include "LinearAlgebra/VectorN.h"



namespace Fresnel
{
	// Dielectrics
	
	template<unsigned D>
	void fresnel_reflectance(const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n, Real &Rp, Real &Rs)
	{
		Real m_n = n;
		Real m_n2 = m_n*m_n;

		Real costhita = std::min<Real>(dot_abs(omega_o, normal), Real(1.));
		Real sinthita = sqrt(Real(1.) - costhita*costhita);
		Real tanthita = sinthita / costhita;

		// Check Brewster & TIR angle
//		Real cos_brewster = cos(atan(n));
		Real cos_tir = (n > Real(1.)) ? Real(-1.) : cos(asin(n));

		Real sin2 = sinthita*sinthita, cos2 = costhita*costhita, tan2 = tanthita*tanthita;

		Real a2 = m_n2 - sin2;
		Real a = sqrtf(a2);

		Rs = (a2 - Real(2.) * a * costhita + cos2) / (a2 + Real(2.) * a*costhita + cos2);
		Rp = (a2 - Real(2.) * a * sinthita*tanthita + sin2*tan2) /
			 (a2 + Real(2.) * a * sinthita*tanthita + sin2*tan2) * Rs;

		if (costhita < cos_tir) {
			Rs = Real(1.);
			Rp = Real(1.);
		}
	}

	template<unsigned D>
	inline Real fresnel_reflectance(const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n)
	{
		Real Rp, Rs;
		fresnel_reflectance(omega_o, normal, n, Rp, Rs);
		return (Rp + Rs) / Real(2.);
	}

	template<unsigned D>
	void fresnel_reflectance(const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n, Real &Rp, Real &Rs, Real &Pp, Real &Ps)
	{
		Real m_n = n, m_n2 = n*n;

		Real costhita = std::min<Real>(dot_abs(omega_o, normal), Real(1.));
		Real sinthita = sqrt(Real(1.) - costhita*costhita);
		Real tanthita = sinthita / costhita;

		// Check Brewster & TIR angle
		Real cos_brewster = cos(atan(n));
		Real cos_tir = (n > Real(1.)) ? Real(-1.) : cos(asin(n));

		Real sin2 = sinthita*sinthita, cos2 = costhita*costhita, tan2 = tanthita*tanthita;

		Real a2 = m_n2 - sin2;
		Real a = sqrt(a2);

		Rs = (a2 - Real(2.) * a * costhita + cos2) / (a2 + Real(2.) * a*costhita + cos2);
		Rp = (a2 - Real(2.) * a * sinthita*tanthita + sin2*tan2) /
			 (a2 + Real(2.) * a * sinthita*tanthita + sin2*tan2) * Rs;

		// Travelling to a denser medium
		if (cos_tir < 0.)
		{
			Ps = M_PI;
			Pp = (costhita > cos_brewster) ? 0. : M_PI;
		}
		else
		{
			if (costhita > cos_tir)
			{
				Ps = 0.;
				Pp = (costhita > cos_brewster) ? M_PI : 0.;
			}
			else // Total Internal Reflection
			{
				Rs = 1.;
				Rp = 1.;

				// From [Azzam 2004]
				Real N = 1. / m_n;
				Real inner = sqrt(N*N*sinthita*sinthita - 1.);
				Pp = 2. * atan2(N*inner, costhita);
				Ps = 2. * atan2(inner, (N*costhita));
			}
		}
	}

	// Conductors
	
	template<unsigned D>
	void fresnel_reflectance( const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n, const Real k, Real &Rp, Real &Rs )
	{
		Real m_n = n;
		Real m_k = k;

		Real m_n2 = m_n*m_n;
		Real m_k2 = m_k*m_k;
		
		Real costhita = dot_abs(omega_o,normal);
		Real sinthita = sqrt(1-costhita*costhita);
		Real tanthita = sinthita / costhita;

		Real sin2 = sinthita*sinthita, cos2 = costhita*costhita, tan2 = tanthita*tanthita;
		    
		Real sqnk = sqrt((m_n2-m_k2-sin2)*(m_n2-m_k2-sin2)+4*m_n2*m_k2);
		Real a2 = (sqnk+m_n2-m_k2-sin2)*.5; a2 = a2<0.?0:a2;
		Real b2 = (sqnk-m_n2+m_k2+sin2)*.5; b2 = b2<0.?0:b2;
		
		Real a = sqrt(a2);
//		Real b = sqrt(b2);
		    
		Rs = (a2+b2-2*a*costhita+cos2)/(a2+b2+2*a*costhita+cos2),
		Rp = (a2+b2-2*a*sinthita*tanthita+sin2*tan2) /
			 (a2+b2+2*a*sinthita*tanthita+sin2*tan2) * Rs;	
	}

	template<unsigned D>
	inline Real fresnel_reflectance( const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n, const Real k)
	{
		Real Rp, Rs;
		fresnel_reflectance(omega_o, normal, n, k, Rp, Rs);
		return (Rp+Rs)/2;
	}

	template<unsigned D>
	void fresnel_reflectance(const VectorN<D> &omega_o, const VectorN<D> &normal, const Real n, const Real k, Real &Rp, Real &Rs, Real &Pp, Real &Ps)
	{
		Real m_n = n, m_k = k, m_n2 = n*n, m_k2 = k*k;

		Real costhita = dot_abs(omega_o, normal);
		Real sinthita = sqrtf(1 - costhita*costhita);
		Real tanthita = sinthita / costhita;

		Real sin2 = sinthita*sinthita, cos2 = costhita*costhita, tan2 = tanthita*tanthita;

		Real sqnk = sqrtf((m_n2 - m_k2 - sin2)*(m_n2 - m_k2 - sin2) + 4 * m_n2*m_k2);
		Real a2 = (sqnk + m_n2 - m_k2 - sin2)*.5; a2 = a2<0. ? 0 : a2;
		Real b2 = (sqnk - m_n2 + m_k2 + sin2)*.5; b2 = b2<0. ? 0 : b2;

		Real a = sqrtf(a2),
			 b = sqrtf(b2);

		Rs = (a2 + b2 - 2 * a*costhita + cos2) / (a2 + b2 + 2 * a*costhita + cos2),
		Rp = (a2 + b2 - 2 * a*sinthita*tanthita + sin2*tan2) /
			 (a2 + b2 + 2 * a*sinthita*tanthita + sin2*tan2) * Rs;
		
		Real tPs = 2*costhita/(cos2-a2-b2),
			 tPp = 2 * b*costhita * ((m_n2 - m_k2)*b - 2 * m_n*m_k*a)/
				((m_n2 + m_k2)*(m_n2 + m_k2) * cos2 - a2 - b2);

		Pp = atan(tPp);
		Ps = atan(tPs);

		/*
		dt_2pi = wl/c;

		Tp = abs(atan(tPp))/(2*pi)*dt_2pi; % In seconds
		Ts = abs(atan(tPs))/(2*pi)*dt_2pi; % In seconds
		*/
	}
}; //Fresnel

#endif //_FRESNEL_H_
