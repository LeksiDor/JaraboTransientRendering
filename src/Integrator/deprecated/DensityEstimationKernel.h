#ifndef _DENSITY_ESTIMATION_KERNEL_H_
#define _DENSITY_ESTIMATION_KERNEL_H_

#include "Utils/globals.h"
#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Intersection.h"

namespace DensityEsimationKernels
{
	//=============================================================
	template<int D>
	class Perlin
	{
	public:
		Real operator()(const Real t, const Real bandwidth)const
		{
			Real tn = t/bandwidth;
			Real k_t = 1.-6.*pow(tn,5)+15.*pow(tn,4)-10.*pow(tn,3);
			return c(bandwidth)*k_t;
		}
		Real c(const Real bandwidth)const;
	};
	//-------------------------------------------------------------
	template<>
	Real Perlin<2>::c(const Real bandwidth)const
	{
		return 7./(2.*M_PI*bandwidth*bandwidth);
	}
	//-------------------------------------------------------------
	template<>
	Real Perlin<1>::c(const Real bandwidth)const
	{
		return 1/bandwidth;
	}

	//=============================================================
	template<int D>
	class Box
	{
	public:
		Real operator()(const Real t, const Real bandwidth)const
		{
			return c(bandwidth);
		}
		Real c(const Real bandwidth)const;
	};
	//-------------------------------------------------------------
	template<>
	Real Box<2>::c(const Real bandwidth)const
	{
		return 1./(M_PI*bandwidth*bandwidth);
	}
	//-------------------------------------------------------------
	template<>
	Real Box<1>::c(const Real bandwidth)const
	{
		return 1/(2*bandwidth);
	}

}; //DensityEsimationKernels


#endif //_DENSITY_ESTIMATION_KERNEL_H_