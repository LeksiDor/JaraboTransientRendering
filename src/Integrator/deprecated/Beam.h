/*
 *  Beam.h
 *  BunnyKiller
 *
 *  Created by Julio on 10/06/13.
 *
 */

#ifndef __BEAM_H
#define __BEAM_H

#ifdef _WIN32 || _WIN64				// Conclict between WinDef.h definitions for min, max and std::min, std::max
	#define NOMINMAX
#endif

#include "Utils/Timer.h"
#include "Photon.h"
#include "LinearAlgebra/Vector2.h"
#include "LinearAlgebra/Vector3.h"
#include "RayTracing/Ray.h"
#include "Media/Medium.h"
#include "RayTracing/Intersection.h"
#include "RadianceSample.h"
#include "RayTracing/AABB.h"
#include "Geometry/Aggregate/BVH.h"
#include <math.h>


enum BeamBlurRegion { BR1D = 1, BR2D = 2, BR3D = 3};

#ifndef _MAX_DIVISIONS_BEAM_
#define _MAX_DIVISIONS_BEAM_ 50
#endif 

using namespace std;
/** 2D Line. */
template<int D, BeamBlurRegion BR>
class Beam: public Ray<D>
{
	Real radius, cross_section;
	Photon<D> *photon;
	AABB<D> bb;
	void setup();

	Real poly (Real *K, Real x, const short int order) const
	{
		// K is a 0..order vector (i.e. k is order+1 size) with the poly coefficients
		// e.g.: 10x^5 + 0.3x^4 + x^2 + 8 would be K = {8,0,1,0,0.3,10}
		
		Real y=K[0];
		for (int i=1; i <= order; i++)
			if (K[i] != 0) y += K[i]*pow(x,i);
		
		return y;
	}
	Real solve_secant_cubic (Real *K, const Real start, const Real end) const
	{
		Real x0, x1, x2 = start, x3 = end;
		Real y0,y1,y2,y3;
		do
		{
			x0 = x2;
			x1 = x3;
			y0 = poly(K,x0,3);
			y1 = poly(K,x1,3);
			x2 = x1 - y1*(x1 - x0)/(y1 - y0);
			y2 = poly(K,x2,3);
			x3 = x2 - y2*(x2 - x1)/(y2 - y1);
			y3 = poly(K,x3,3);
		} while (fabs(y3) > 1e-4);
		return x3;
	}
public:
	inline Spectrum power_at (Real _t) const
	{
		Spectrum u_t = (this->get_medium() != NULL)?this->get_medium()->get_extinction(this->get_origin()):Spectrum(0.);
		return exp(-u_t*_t)*photon->m_power;
	}
	
	inline Real time_of_flight(Real distance) const
	{
		return distance/0.000299792458;
	}
	inline Real get_time_of_flight() const
	{
		return photon->m_time;
	}
	inline Photon<D>* get_photon()const{return this->photon;}
	void scale_power (const Real scale) const
	{
		this->photon->m_power *= scale;
	}
	
	Beam(Photon<D> *_p, Real _radius, Real _t=0., bool normalize_d=true, int _level=0,
		 Real _n = DEFAULT_REFRACTION_INDEX, Medium<D> *_medium = 0)
	: Ray<D>(_p->m_position, _p->m_direction, normalize_d, _level, _n, _medium), photon(_p), radius(_radius), cross_section(M_PI*radius*radius)
	
	{
		Ray<D>::set_parameter(_t);
		setup();
	}
	
	Beam()
	: Ray<D>()
	{}
	
	VectorN<4> get_intersection_points (const Ray<D>& r) const;
	Spectrum contribution (const Ray<D>& r) const;
	Spectrum contribution (Intersection<D>* it) const;
	void contribution (const Ray<D>& r, RadianceSample& sample) const;
	void contribution (const VectorN<D>& p, RadianceSample& sample) const;

	int subdivide (std::vector<Beam<D,BR> >& beams, const int max_div = 30)
	{
		int nb = 0;
		Real bstep = 1.;//this->get_parameter()/nb;
		/*nb = (int)ceil(this->get_parameter()/bstep);
		for (int i=0; i < nb; i++)
		{
			Real p_time = photon->m_time; // TO-DO: Has to be fixed for transient (just in case we do create extra photons)
			Photon<D> *new_photon = new Photon<D>(photon->m_position + photon->m_direction*i*bstep, 
												  photon->m_direction, 
												  power_at(i*bstep), 
												  p_time, //to be fixed for transient
												  photon->m_level);
			Beam<D,BR> *b = new Beam<D,BR>(new_photon, radius, (i == nb-1?this->get_parameter() - i*bstep:bstep), true, photon->m_level,
									 DEFAULT_REFRACTION_INDEX, this->get_medium());
			beams.push_back(*b);
		}*/

		// Note: In theory, we cannot use a variable to specify the size of an static vector.
		// I've added this crap to handle this...
		Real divs = max_div;
		if( max_div > _MAX_DIVISIONS_BEAM_ )
			divs = _MAX_DIVISIONS_BEAM_;

		Real t = this->get_parameter();
		Real step = 1./(Real)divs;
		Real bsteps[_MAX_DIVISIONS_BEAM_];
		int i;
		for (i=0;i < divs; i++)
			bsteps[i] = 1. - step*i;
		
		Real log_bsteps[_MAX_DIVISIONS_BEAM_];
		for (i=0;i < divs; i++)
			log_bsteps[i] = t*(1+bsteps[i]*(log(bsteps[i])-1));
		
		for (i=0; i < divs; i++)
		{
			Real p_time = photon->m_time; // TO-DO: Has to be fixed for transient (just in case we do create extra photons)
			Photon<D> *new_photon = new Photon<D>(photon->m_position + photon->m_direction*log_bsteps[i], 
												  photon->m_direction, 
												  power_at(log_bsteps[i]), 
												  p_time, //to be fixed for transient
												  photon->m_level);
			//Beam<D,BR> *b = new Beam<D,BR>(new_photon, radius, (i==0?log_bsteps[i]:(log_bsteps[i]-log_bsteps[i-1])), true, photon->m_level,
			//							   DEFAULT_REFRACTION_INDEX, this->get_medium());
			beams.push_back(Beam<D,BR>(new_photon, radius, (i==max_div-1?(t-log_bsteps[i]):(log_bsteps[i+1]-log_bsteps[i])), true, photon->m_level,
									   DEFAULT_REFRACTION_INDEX, this->get_medium()));
		}
		
		return nb;
	}
	
	bool intersect(const Ray<D>& r, Intersection<D>& it, float max_distance) const
	{	
		VectorN<4> pts = get_intersection_points(r);
		bool hit = !((pts[0] <= 0.0 && pts[1] <= 0.0) 
					 || (pts[3] <= 0.0 && pts[2] <= 0.0)
					 || (pts[1] > r.get_parameter() && pts[0] > r.get_parameter())
					 || (pts[3] > this->get_parameter() && pts[2] > this->get_parameter()));
		if (hit) 
		{
			//it.b = this;
			//it.beam_it = pts;
		}
		return hit;
	}
	
	bool intersect(const VectorN<D>& p) const
	{	
		VectorN<D> p_origin = p-this->get_origin();
		
		if (this->get_direction().dot(p_origin.normalize()) < 0.) return false;
		
		Real dist2 = this->get_direction().cross(p-this->get_origin()).length2();
		return dist2 < radius*radius;
	}
	
	bool intersect(const Ray<D>& r, float& max_distance) const
	{	
		bool hit = get_bounding_box().intersect(r,max_distance);
		if (!hit) return false;
		VectorN<4> pts = get_intersection_points(r);
		hit = !((pts[0] <= 0.0 && pts[1] <= 0.0) 
					|| (pts[3] <= 0.0 && pts[2] <= 0.0)
					|| (pts[1] > r.get_parameter() && pts[0] > r.get_parameter())
					|| (pts[3] > this->get_parameter() && pts[2] > this->get_parameter()));
		return hit;
	}
	
	inline Real get_radius() const
	{	return radius;}
	
	Real get_intersection_cost()const
	{
		return 7000;
	}
	
	AABB<D> get_bounding_box() const
	{
		return bb;
	}
	
	VectorN<D> get_center()const
	{
		return this->get_origin()+0.5*this->get_parameter()*this->get_direction();
	}
	
	inline void set_parameter(const Real par)
	{
		Ray<D>::set_parameter(par);
		setup();
	}
};
template<>
VectorN<4> Beam<3,BR1D>::get_intersection_points(const Ray<3>& r) const
{
	Real A,B,C,tr,tb;
	VectorN<4> v;
	Vector3 ro = r.get_origin(), rw = r.get_direction(), bo = this->get_origin(), bw = this->get_direction();

	A = dot(ro-bo, bw);
	B = dot(ro-bo, rw);
	C = dot(bw, rw);
	
	// CROSS PRODUCT FORMULATION
	if (fabs(C*C - 1.) < 1e-5) return VectorN<4>(-1.); //beam and ray are parallel
	else 
	{
		Real dist2 = (bo - ro - (dot(bo - ro, bw))*bw).length2(); //distance from camera origin to beam
		if (dist2 < radius*radius) //if ro is within the kernel radius
		{	
			/*	PLANE 1: Contains {bo, bo+bw, ro}
				PLANE 2: Orthogonal to PLANE 1 and containing beam ray (therefore containing bo)
				tr: Intersection distance of 
			 */
			Vector3 wx1 = cross(bw,bo-ro).normalize();	//perpendicular to PLANE 1 defined by {bo, bo+bw, ro}
			Vector3 wx2 = cross(wx1, bw);	//perpendicular to PLANE 2 (orthogonal to PLANE 1 and containing beam ray)

			tr = dot(bo-ro, wx2)/dot(rw,wx2);
			tb = A + tr*C;
		}
		else //do as usual
		{
			tr = (B-C*A)/(C*C - 1.);
			tb = A + C*tr;
		}
	}
	
	Vector3 RtoB = bo + bw*tb - ro - rw*tr;
	Real len2 = RtoB.length2();
	if (len2 > radius*radius) return VectorN<4>(-1.);
	
	v[0] = tr; v[1] = tb; v[2] = len2;
	return v;
}

template<>
VectorN<4> Beam<3,BR2D>::get_intersection_points(const Ray<3>& r) const
/* Returns Vector4(t_r-, t_r+, t_b-, t_b+) where 
 t_r- : closest intersection point distance between ray and beam blur region
 t_r+ : furthest intersection point distance between ray and beam blur region
 t_b- : projection of t_r- over the beam ray 
 t_b+ : projection of t_r+ over the beam ray 
 WITHOUT CLIPPING!!
 */
{
	Vector3 rw, ro, bw, bo, u1, u2;
	Real k1, k2, u_A, u_B, u_C, discr;
	Real tr0, tr1, tb0, tb1;
	const Real eps = 1e-5;
	
	bo = this->get_origin();
	bw = this->get_direction();
	ro = r.get_origin();
	rw = r.get_direction();
	
	k1 = dot(rw, bw);
	if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
	{	
		return VectorN<4>(-1.); 
	}
	
	k2 = dot(ro,bw) - dot(bo, bw);
	u1 = ro - bo - k2*bw;
	u2 = rw - k1*bw;
	u_A = dot(u2, u2);
	u_B = 2.*dot(u1, u2);
	u_C = dot(u1, u1) - radius*radius;
	discr = u_B*u_B - 4.*u_A*u_C;
	
	if (discr < eps) //there's no solution to the quadratic equation or there's only one solution => no intersection
		return VectorN<4>(-1.);
	else 
	{
		Real discr_sqrt = sqrt(discr);
		tr1 = (-u_B + discr_sqrt)/(2.*u_A); 
		tr0 = (-u_B - discr_sqrt)/(2.*u_A);
		tb0 = tr0*k1 + k2;
		tb1 = tr1*k1 + k2;
		if (fabs(tr1 - tr0) < eps) 
			return VectorN<4>(-1.);
	}

	VectorN<4> v;
	v[0] = tr0; v[1] = tr1; v[2] = tb0; v[3] = tb1;
	return v;
}	


template<>
bool Beam<3,BR1D>::intersect(const Ray<3>& r, float& max_distance) const
{
	VectorN<4> pts = get_intersection_points(r);
	return pts[0] >= 0. && pts[1] >= 0.;
}

template<>
VectorN<4> Beam<2,BR2D>::get_intersection_points(const Ray<2>& r) const
/* Returns Vector4(t_r-, t_r+, t_b-, t_b+) where 
		t_r- : closest intersection point distance between ray and beam blur region
		t_r+ : furthest intersection point distance between ray and beam blur region
		t_b- : projection of t_r- over the beam ray 
		t_b+ : projection of t_r+ over the beam ray 
	WITHOUT CLIPPING!!
*/
{
	Real c_theta = dot(get_direction(), r.get_direction()), abs_c_theta;
	Real rdx, rdy, bdx, bdy; //direction
	Real rox, roy, box, boy; //origin
	Real den, tr, tr0, tr1, tb, tb0, tb1, h;

	abs_c_theta = fabs(c_theta);
	
	rdx = r.get_direction()[0]; 
	rdy = r.get_direction()[1];
	rox = r.get_origin()[0]; 
	roy = r.get_origin()[1];
	
	bdx = this->get_direction()[0];
	bdy = this->get_direction()[1];
	box = this->get_origin()[0];
	boy = this->get_origin()[1];
	
	den = bdy*rdx - bdx*rdy;
	
	if (fabs(den) < 1e-5 || (1.-abs_c_theta*abs_c_theta) < 1e-5)
		return VectorN<4>(-1.0); //no intersection
	
	//intersection between beam ray and ray
	tr = (roy*bdx - boy*bdx + box*bdy - rox*bdy)/den;
	//intersection with blur region
	h = radius/sqrtf(1.-abs_c_theta*abs_c_theta);

	tr1 = tr+h;
	tr0 = tr-h;
	
	Vector2 rtr = r.get_origin() + tr*r.get_direction();
	//perpendicular projections of tr1 and tr0 over beam ray
	tb = sqrtf(dot(get_origin()-rtr, get_origin()-rtr));
	if (fabs(get_origin()[0] + tb*get_direction()[0] - rtr[0]) > 0.0001 
		|| fabs(get_origin()[1] + tb*get_direction()[1] - rtr[1]) > 0.0001)
		tb = -tb;
	
	tb0 = tb - h*c_theta;
	tb1 = tb + h*c_theta;
	VectorN<4> v;
	v[0] = tr0; v[1] = tr1; v[2] = tb0; v[3] = tb1;
	return v;
}

template<int D, BeamBlurRegion BR>
Spectrum Beam<D, BR>::contribution (Intersection<D>* it) const
{
	Ray<D> r = it->get_ray();
	VectorN<4> it_points = it->beam_it;
	Real c_theta = dot(this->get_direction(), r.get_direction()), abs_c_theta;
	Real tr0, tr1, tb0, tb1;
	Spectrum f, aux(1.); //phase function

	
	if (BR == BR2D) //2D(1) blur region, perpendicular to beam
	{
		tr0 = it_points[0]; tr1 = it_points[1]; tb0 = it_points[2], tb1 = it_points[3];
		
		abs_c_theta = fabs(c_theta);
		
		//CLIPPING INTERSECTIONS WITHIN BEAM AND RAY BOUNDARIES)
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > r.get_parameter() && tr0 > r.get_parameter())
			|| (tb1 > this->get_parameter() && tb0 > this->get_parameter()))
			return Spectrum(0.0);
		
		/*Real dist2 = (this->get_origin() - r.get_origin() - (dot(this->get_origin() - r.get_origin() , this->get_direction()))*this->get_direction()).length2(); //distance from camera origin to beam
		 if (dist2 < radius*radius) aux= Spectrum(0., 3., 0.);*/
		
		//CHECK BOUNDS OF RAY DISTANCES
		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > r.get_parameter())
		{
			tb1 = tb1 - (tr1 - r.get_parameter())*c_theta;
			tr1 = r.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
			tb0 = 0.0f;
		}
		else if (tb0 > this->get_parameter())
		{
			tr0 = tr0 - (tb0-this->get_parameter())/c_theta;
			tb0 = this->get_parameter();
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > this->get_parameter())
		{
			tr1 = tr1 - (tb1-this->get_parameter())/c_theta;
			tb1 = this->get_parameter();
		}
		
		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
			return Spectrum(0.);
		
		//DO THE MATH
		Medium<D> *m = this->get_medium();
		Spectrum num, den, u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		if (c_theta < 0.) //tb0 > tb1
		{	
			num = exp(-u_t*(tr0-tr1)*(abs_c_theta - 1.)) - 1.;
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta - 1.);
		}
		else // tb0 < tb1
		{	
			num = 1. - exp(-u_t*(tr1-tr0)*(abs_c_theta + 1.));
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta + 1.);
		}
		
		Spectrum contr = num/den;
		
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		return f*contr;
	}
	else if (BR == BR1D) //1D blur region
	{
		Real s_theta;
		Spectrum aux(1.);
		
		tr0 = it_points[0]; 
		tb0 = it_points[1];
		if (tr0 < 0. || tb0 < 0. || tr0 > r.get_parameter() || tb0 > this->get_parameter())
			return Spectrum(0.);
		
		s_theta = sqrtf(1.-c_theta*c_theta);
		
		Medium<D> *m = this->get_medium();
		Spectrum u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		
		Real K = 1., x;
		//KERNEL
		/*x = it_points[2]/(radius*radius);
		 x = 1.-x;
		 x *= x;
		 K = (15./16.)*x;
		 */
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		Spectrum contr = exp(-u_t*tb0)*exp(-u_t*tr0)/s_theta;
		float contr_max = max(contr[0], max(contr[1], contr[2]));
		if (contr_max > 1.) contr /= contr_max;
		
		return f*contr;//(contr.avg()>1. ? Spectrum(1.) : contr);
	}
	
}


template<int D, BeamBlurRegion BR>
Spectrum Beam<D, BR>::contribution (const Ray<D>& r) const
{
	Real c_theta = dot(this->get_direction(), r.get_direction()), abs_c_theta;
	Real tr0, tr1, tb0, tb1;
	VectorN<4> it_points;
	Spectrum f, aux(1.); //phase function
	
	it_points = get_intersection_points(r);
	if (BR == BR2D) //2D(1) blur region, perpendicular to beam
	{
		tr0 = it_points[0]; tr1 = it_points[1]; tb0 = it_points[2], tb1 = it_points[3];
		
		abs_c_theta = fabs(c_theta);
		
		//CLIPPING INTERSECTIONS WITHIN BEAM AND RAY BOUNDARIES)
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > r.get_parameter() && tr0 > r.get_parameter())
			|| (tb1 > this->get_parameter() && tb0 > this->get_parameter()))
			return Spectrum(0.0);
		
		/*Real dist2 = (this->get_origin() - r.get_origin() - (dot(this->get_origin() - r.get_origin() , this->get_direction()))*this->get_direction()).length2(); //distance from camera origin to beam
		 if (dist2 < radius*radius) aux= Spectrum(0., 3., 0.);*/
		
		//CHECK BOUNDS OF RAY DISTANCES
		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > r.get_parameter())
		{
			tb1 = tb1 - (tr1 - r.get_parameter())*c_theta;
			tr1 = r.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta;
			tb0 = 0.0f;
		}
		else if (tb0 > this->get_parameter())
		{
			tr0 = tr0 - (tb0-this->get_parameter())/c_theta;
			tb0 = this->get_parameter();
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > this->get_parameter())
		{
			tr1 = tr1 - (tb1-this->get_parameter())/c_theta;
			tb1 = this->get_parameter();
		}
		
		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
			return Spectrum(0.);
		
		//DO THE MATH
		Medium<D> *m = this->get_medium();
		Spectrum num, den, u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		
		//with no ray extinction
		//num = exp(-tb0*u_t) - exp(-u_t*(tb0 + c_theta*(tr1 - tr0)));
		//den = u_t*c_theta;
		
		num = exp(-u_t*(tb0 + tr0)) - exp(-u_t*(tb0+c_theta*(-tr0+tr1)+tr1));
		den = u_t*(c_theta+1.);
		
		//printf("%f %f %f %f %f\n", tr0,tr1,tb0,tb1, fabs(tb0-tb1));
		/*if (c_theta < 0.) //tb0 > tb1
		{	
			num = exp(-u_t*(tr0-tr1)*(abs_c_theta - 1.)) - 1.;
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta - 1.);
		}
		else // tb0 < tb1
		{	
			num = 1. - exp(-u_t*(tr1-tr0)*(abs_c_theta + 1.));
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta + 1.);
		}*/
		
		Spectrum contr = num/den;
		
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		return f*contr;
	}
	else if (BR == BR1D) //1D blur region
	{
		Real s_theta;
		Spectrum aux(1.);
		
		tr0 = it_points[0]; 
		tb0 = it_points[1];
		if (tr0 < 0. || tb0 < 0. || tr0 > r.get_parameter() || tb0 > this->get_parameter())
			return Spectrum(0.);
		
		s_theta = sqrtf(1.-c_theta*c_theta);
		
		Medium<D> *m = this->get_medium();
		Spectrum u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		
		Real K = 1., x;
		//KERNEL
		/*x = it_points[2]/(radius*radius);
		 x = 1.-x;
		 x *= x;
		 K = (15./16.)*x;
		 */
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		Spectrum contr = exp(-u_t*tb0)*exp(-u_t*tr0)/s_theta;
		float contr_max = max(contr[0], max(contr[1], contr[2]));
		if (contr_max > 1.) contr /= contr_max;
		
		return f*contr;//(contr.avg()>1. ? Spectrum(1.) : contr);
	}
	
}
template<int D, BeamBlurRegion BR>
void Beam<D, BR>::contribution (const VectorN<D>& p, RadianceSample& sample) const
{
	if (BR == BR2D)
	{
		VectorN<D> p_origin = p - this->get_origin();
		Real tb = this->get_direction().dot(p_origin);
		if (tb > this->get_parameter())
		{
			sample = RadianceSample();
			return;
		}
		Medium<D> *m = this->get_medium();
		Spectrum u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		//no phase function, integrating the radiance arriving at p, not outgoing
		sample = RadianceSample(exp(-u_t*tb), get_time_of_flight() + this->time_of_flight(tb)*this->get_refraction_index(), this->get_level(), 0.);
	}
}
template<int D, BeamBlurRegion BR>
void Beam<D, BR>::contribution (const Ray<D>& r, RadianceSample& sample) const
{
	Real c_theta = dot(this->get_direction(), r.get_direction()), abs_c_theta;
	Real tr0, tr1, tb0, tb1;
	VectorN<4> it_points;
	Spectrum f, aux(1.); //phase function
	
	it_points = get_intersection_points(r);
	
	if (BR == BR2D) //2D(1) blur region, perpendicular to beam
	{
		tr0 = it_points[0]; tr1 = it_points[1]; tb0 = it_points[2], tb1 = it_points[3];
		
		abs_c_theta = fabs(c_theta);
		
		//CLIPPING INTERSECTIONS WITHIN BEAM AND RAY BOUNDARIES)
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > r.get_parameter() && tr0 > r.get_parameter())
			|| (tb1 > this->get_parameter() && tb0 > this->get_parameter()))
		{
			sample = RadianceSample();
			return;
		}
		/*Real dist2 = (this->get_origin() - r.get_origin() - (dot(this->get_origin() - r.get_origin() , this->get_direction()))*this->get_direction()).length2(); //distance from camera origin to beam
		 if (dist2 < radius*radius) aux= Spectrum(0., 3., 0.);*/
		
		//CHECK BOUNDS OF RAY DISTANCES
		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > r.get_parameter())
		{
			tb1 = tb1 - (tr1 - r.get_parameter())*c_theta;
			tr1 = r.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
			tb0 = 0.0f;
		}
		else if (tb0 > this->get_parameter())
		{
			tr0 = tr0 - (tb0-this->get_parameter())/c_theta;
			tb0 = this->get_parameter();
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > this->get_parameter())
		{
			tr1 = tr1 - (tb1-this->get_parameter())/c_theta;
			tb1 = this->get_parameter();
		}
		
		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
		{
			sample = RadianceSample();
			return;
		}
		//DO THE MATH
		Medium<D> *m = this->get_medium();
		Spectrum num, den, u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		if (c_theta < 0.) //tb0 > tb1
		{	
			num = exp(-u_t*(tr0-tr1)*(abs_c_theta - 1.)) - 1.;
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta - 1.);
		}
		else // tb0 < tb1
		{	
			num = 1. - exp(-u_t*(tr1-tr0)*(abs_c_theta + 1.));
			den = exp(u_t*(tr0 + tb0))*u_t*(abs_c_theta + 1.);
		}
		
		Spectrum contr = num/den;
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		//sample = RadianceSample(f*contr, get_time_of_flight() + this->time_of_flight(0.5*(tr0+tr1+tb0+tb1))*r.get_refraction_index(), this->get_level(), this->time_of_flight(0.5*(tr0+tr1))*r.get_refraction_index());
        
        sample = RadianceSample(f*contr, 
								get_time_of_flight() + this->time_of_flight(std::min(tr0+tb0, tr1+tb1) + 0.5*fabs(tb1+tr1-tb0-tr0))*r.get_refraction_index(), 
								this->get_level(), 
								this->time_of_flight(0.5*(tr0+tr1))*r.get_refraction_index());
        //printf("c_r + b_r = %f\t (1+cos_theta)*c_r = %f\t aux = %f\n", fabs(tr0 - tr1) + fabs(tb0 - tb1), (1+c_theta)*fabs(tr0 - tr1), tb1+tr1 - tb0 - tr0);
		sample.camera_range = fabs(tr0 - tr1);
		sample.beam_range = fabs(tb1+tr1 - tb0 - tr0);
        //sample.beam_range = fabs(tb0 - tb1);
	}
	else if (BR == BR1D) //1D blur region
	{
		Real s_theta, len2;
		Spectrum aux(1.);
		
		tr0 = it_points[0]; 
		tb0 = it_points[1];
		len2 = it_points[2];
		if (tr0 < 0. || tb0 < 0. || tr0 > r.get_parameter() || tb0 > this->get_parameter())
			return; Spectrum(0.);
		
		s_theta = sqrtf(1.-c_theta*c_theta);
		
		Medium<D> *m = this->get_medium();
		Spectrum u_t = (m)?m->get_extinction(this->get_origin()):Spectrum(0.);
		
		Real K = 1., x;
		//KERNEL
		/*x = it_points[2]/(radius*radius);
		 x = 1.-x;
		 x *= x;
		 K = (15./16.)*x;
		 */
		//PHASE FUNCTION
		f = m->f(r.get_origin()+tr0*r.get_direction(),this->get_direction(),-r.get_direction());
		
		Spectrum contr = exp(-u_t*tb0)*exp(-u_t*tr0)/s_theta;
		float contr_max = max(contr[0], max(contr[1], contr[2]));
		if (contr_max > 1.) contr /= contr_max;
		
		sample = RadianceSample(f*contr, get_time_of_flight() + this->time_of_flight(tr0 + tb0)*r.get_refraction_index(), this->get_level(), this->time_of_flight(tr0)*r.get_refraction_index());
	}
	
}

template<>
void Beam<2,BR2D>::setup()
{
	//bounding box using blur region radius
	float abs_cos = fabs(dot(get_direction(), Vector2(1., 0.))); //increment y
	float abs_sin = sqrtf(1.-abs_cos*abs_cos); //increment x
	abs_cos*=radius;
	abs_sin*=radius;
	Vector2 end = this->get_origin() + this->get_parameter()*this->get_direction();
	
	bb = AABB2D(Vector2(min(this->get_origin()[0], end[0])-abs_sin, min(this->get_origin()[1], end[1])-abs_cos),
				  Vector2(max(this->get_origin()[0], end[0])+abs_sin, max(this->get_origin()[1], end[1])+abs_cos));
}

template<>
void Beam<3,BR2D>::setup()
{
	Vector3 end = this->get_origin()+ this->get_parameter()*this->get_direction();
	bb = AABB3D(Vector3(min(this->get_origin()[0], end[0])-radius, min(this->get_origin()[1], end[1])-radius, min(this->get_origin()[2], end[2])-radius),
				Vector3(max(this->get_origin()[0], end[0])+radius, max(this->get_origin()[1], end[1])+radius, max(this->get_origin()[2], end[2])+radius));
}

template<>
void Beam<3,BR1D>::setup()
{
	Vector3 end = this->get_origin()+ this->get_parameter()*this->get_direction();
	bb = AABB3D(Vector3(min(this->get_origin()[0], end[0]) - radius, min(this->get_origin()[1], end[1]) - radius, min(this->get_origin()[2], end[2]) - radius),
				Vector3(max(this->get_origin()[0], end[0]) + radius, max(this->get_origin()[1], end[1]) + radius, max(this->get_origin()[2], end[2]) + radius));
}

typedef Beam<2,BR2D> Beam2D;
typedef Beam<3,BR2D> Beam3D;
#endif
