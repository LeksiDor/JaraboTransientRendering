#ifndef _TERMINATION_CRITERIA_H_
#define _TERMINATION_CRITERIA_H_

#include "Utils/globals.h"
template<int D> class Intersection;

namespace TerminationCriteria
{
	static int max_nb_bounces = MAX_NB_BOUNCES;
	static int max_splitting_bounces = MAX_NB_BOUNCES;

	struct recursion_depth
	{
		template<int D>
		inline bool operator()(const Intersection<D> &it, Real &pdf)const
		{
			pdf = 1.;
			if( it.get_ray().get_level() < max_nb_bounces )
				return false;

			return true;
		}
	};

	struct russian_roulette
	{
		template<int D>
		inline bool operator()(const Intersection<D> &it, Real &pdf)const
		{
			pdf = 1.;
			if( it.get_ray().get_level() < max_nb_bounces )
				return false;

			return true;
		}
	};

	struct split_and_roussian_roulette: public russian_roulette
	{
		template<int D>
		inline bool operator()(const Intersection<D> &it, Real &pdf)const
		{
			if( it.get_ray().get_level() < max_splitting_bounces ) 
			{
				pdf = 1.;
				return false;
			}
			else
				return russian_roulette::operator()(it, pdf);
		}
	};
};


#endif //_TERMINATION_CRITERIA_H_