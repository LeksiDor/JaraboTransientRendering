/*
 * Copyright (C) 2017, Julio Marco (http://webdiis.unizar.es/~juliom/)
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

#ifndef _VOLUMETRIC_RADIANCE_ESTIMATORS_H_
#define _VOLUMETRIC_RADIANCE_ESTIMATORS_H_

#include "bunnykiller.h"

#include <cmath>
#include <cstdio>
#include <vector>

#include "Color/Spectrum.h"
#include "DensityEstimation/DensityEstimationKernel.h"
#include "Integrator/Particle.h"
#include "Integrator/RadianceSample.h"
#include "Media/Medium.h"
#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include "RayTracing/World.h"

#ifdef _USE_EMBREE_
#define EMBREE_STATIC_LIB
#include "External/Embree/include/embree2/rtcore_geometry_user.h"
#endif

namespace VRE{
	template<unsigned D, class Radiance, class RadianceAttenuation>
	class VolumetricRadianceEstimator
	{
	protected:
		using RadianceSampleR = RadianceSample<Radiance>;
		using ParticleR = Particle<D, Radiance, RadianceAttenuation>;
	protected:
		unsigned int m_nb_samples_ray;
		Real m_inv_nb_samples_ray;
		Real m_search_radius;
		Real m_time_exp;
		World<D, Radiance> *m_world;
	public:
		VolumetricRadianceEstimator() :
			m_nb_samples_ray(1), m_inv_nb_samples_ray(1.)
		{}

		VolumetricRadianceEstimator(unsigned int nb_samples) :
			m_nb_samples_ray(nb_samples), m_inv_nb_samples_ray(1./nb_samples)
		{}
		
		virtual ~VolumetricRadianceEstimator()
		{}

		// Steady
		// For rays
		virtual Radiance operator()(const Ray<D> &c) const = 0;
		virtual Radiance operator()(const Ray<D> &c, const ParticleR &p, const typename ParticleR::PIsect &isec, Real &pdf) const = 0;
		// For points
		virtual Radiance operator()(const VectorN<D> &c) const = 0;
		virtual Radiance operator()(const VectorN<D> &c, const ParticleR &p, const typename ParticleR::PIsect &isec, Real &pdf) const = 0;

		// Transient
		virtual void operator()(const Ray<D> &c, std::vector<RadianceSampleR> &samples) const = 0;
		
	public:
		Real sample_media_distance(const Ray<D> &r, Real &p) const;

		virtual Real K(const Ray<D> &c, const ParticleR &p) const = 0;

		virtual Real K(const Real r) const {
			switch (D) {
				case 1: return 1./(2.*r);
				case 2: return 1./(M_PI*r*r);
				case 3: return 3./(4.*M_PI*r*r*r);
			}
		};
		
		virtual bool intersect(const Ray<D> &c, const ParticleR &p, typename ParticleR::PIsect &it) const
		{
			return false;
		};

		virtual bool intersect(const VectorN<D> &c, const ParticleR &p, typename ParticleR::PIsect &it) const
		{
			return false;
		};
	};
#if 0
	/*******************************************************/
	// Point data - Point query 3D (only for a 3D space)
	template<class Radiance, class RadianceAttenuation>
	class PP3D: public VolumetricRadianceEstimator<3, Radiance, RadianceAttenuation>
	{
	protected:
		using VRE = VolumetricRadianceEstimator<3, Radiance, RadianceAttenuation>;
		using RadianceSampleR = typename VRE::RadianceSampleR;
		using ParticleR = typename VRE::ParticleR;
	private:
		KDTree<ParticleR, 3> m_global_vpoints_tree, m_caustic_vpoints_tree;
		typedef typename KDTree<ParticleR, 3>::Node KDTreeNode;
	public:
		PP3D() :
			VRE(1)
		{}
		PP3D(std::vector<ParticleR> particles, unsigned int nb_samples = 1) :
			VRE(nb_samples)
		{
			typename std::vector<ParticleR>::iterator it;
			unsigned int nb_particles = 0;
			VRE::m_search_radius = (*particles.begin()).m_r;

			for (it=particles.begin(); it != particles.end(); it++) {
				if (it->m_type == ParticleR::VolPoint) {
					nb_particles++;
					m_global_vpoints_tree.store(std::vector<Real>(it->m_v->get_vertex_position().m_data, 
						it->m_v->get_vertex_position().m_data+3), *it);
				}
			}
			if (nb_particles)
				m_global_vpoints_tree.balance();
		}

		virtual Radiance operator()(const Ray<3> &c) const override
		{ 
			Spectrum I(0.);
			if (!m_global_vpoints_tree.is_empty()) {
				Real pr, pp, dist;
				const Medium<3>* m = c.get_medium();

				for (unsigned int i = 0; i < VRE::m_nb_samples_ray; i++) {
					dist = m->sample_mean_free_path(c, pr);
					if (dist > c.get_parameter())
						continue;
					
					Ray<3> b(c.get_origin() + dist*c.get_direction(), c.get_direction());
					std::vector<Real> ori {b.get_origin()[0], b.get_origin()[1], b.get_origin()[2]};

					std::vector<const typename KDTree<ParticleR, 3>::Node*> photons;
					m_global_vpoints_tree.find(ori, VRE::m_search_radius, &photons);
		
					typename ParticleR::PIsect isec;
					for (const KDTreeNode* ph: photons) {
						I += (*this)(b, ph->data(), isec, pp);
					}
				}
				I /= Real(VRE::m_nb_samples_ray);
			}
			return I;
		}

		virtual Radiance operator()(const Ray<3> &b, const Particle<3> &p, const typename Particle<3>::PIsect &isec, Real &pdf) const override
		{ 
			Spectrum f, u_t, u_s;
			const Medium<3>* m = dynamic_cast<const Path<3>::VolumeVertex*>(p.m_v)->m_medium;
			u_t = m->get_extinction(p.m_v->m_position);
			f = m->f(b.get_origin(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*K(b,p)*f/u_t;
		}

		Real K(const Ray<3> &b, const Particle<3> &p) const
		{
			return 3./(4.*M_PI*p.m_r*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<3> &b,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
		{
			it.did_hit = p.dist(b.get_origin()) < p.m_r;
			return it.did_hit;
		}
	};

	/*******************************************************/
	// Point data - Point query 2D (only for a 2D space)
	template<int D>
	class PP2D: public VolumetricRadianceEstimator<2>
	{
	private:
		KDTree<Particle<2>, 2> m_global_vpoints_tree, m_caustic_vpoints_tree;
	public:
		PP2D():VolumetricRadianceEstimator(1){}
		PP2D(std::vector<Particle<2> > particles, unsigned int nb_samples = 1):VolumetricRadianceEstimator(nb_samples){
			std::vector<Particle<2> >::iterator it;
			unsigned int nb_particles = 0;
			m_search_radius = (*particles.begin()).m_r;
			for (it=particles.begin(); it != particles.end(); it++)
			{
				if (it->m_type == Particle<2>::VolPoint)
				{	
					nb_particles++;
					m_global_vpoints_tree.store(std::vector<Real>(it->m_v->get_vertex_position().m_data, 
						it->m_v->get_vertex_position().m_data+2), *it);
				}
			}
			if (nb_particles)
				m_global_vpoints_tree.balance();
		}

		Spectrum operator()(const Ray<2> &c) const
		{ 
			Spectrum I(0.);
			if ( !m_global_vpoints_tree.is_empty() )
			{
				Real pr, pp, dist;
				unsigned int ns = 0;
				const Medium<2>* m = c.get_medium();
				for (unsigned int i = 0; i < m_nb_samples_ray; i++)
				{
					dist = m->sample_mean_free_path(c, pr);
					if (dist > c.get_parameter()) continue;
					ns++;
					Ray<2> b (c.get_origin()+dist*c.get_direction(), c.get_direction());
			
					std::vector<const Particle<2> *>::const_iterator it;
			
					typename Particle<2>::PIsect isec;
			
					typedef typename KDTree<Particle<2>, 2>::Node KDTreeNode;
					std::list<const KDTreeNode*> photons;
					typename std::list<const KDTreeNode*>::iterator ph;
					
					m_global_vpoints_tree.find(std::vector<Real>(b.get_origin().m_data,b.get_origin().m_data+2), m_search_radius, &photons);
					Real pa;
		
					for( ph = photons.begin(); ph != photons.end(); ph++ )
						I += (*this)(b,((*ph)->data()), isec, pp);
				}
				I/=m_nb_samples_ray;
			}
			return I;
		}

		Spectrum operator()(const Ray<2> &b, const Particle<2> &p, const typename Particle<2>::PIsect &isec, Real &pdf) const
		{ 
			Spectrum f, u_t, u_s;
			const Medium<2>* m = dynamic_cast<const Path<2>::VolumeVertex*>(p.m_v)->m_medium;
			u_t = m->get_extinction(p.m_v->m_position);
			f = m->f(b.get_origin(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*K(b,p)*f/u_t;
		}

		Spectrum operator()(const VectorN<2> &b) const
		{ 
			Spectrum I(0.);
			if ( !m_global_vpoints_tree.is_empty() )
			{
				Real pr, pp, dist;
				std::vector<const Particle<2> *>::const_iterator it;
			
				typedef typename KDTree<Particle<2>, 2>::Node KDTreeNode;
				std::list<const KDTreeNode*> photons;
				typename std::list<const KDTreeNode*>::iterator ph;
					
				m_global_vpoints_tree.find(std::vector<Real>(b.m_data,b.m_data+2), m_search_radius, &photons);
				Real pa;
			
				for ( ph = photons.begin(); ph != photons.end(); ph++ )
				{
					Spectrum f, u_t, u_s;
					const Medium<2>* m = dynamic_cast<const Path<2>::VolumeVertex*>((*ph)->data()->m_v)->m_medium;
					u_t = m->get_extinction((*ph)->data()->m_v->m_position);
					u_s = m->get_scattering((*ph)->data()->m_v->m_position);
					I += (*ph)->data()->m_power*K(b,(*ph)->data())*u_s/u_t; 
				}
			}
			return I;
		}

		Real K(const VectorN<2> &b, const Particle<2> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}

		Real K(const Ray<2> &b, const Particle<2> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<2> &b,  const Particle<2> &p, typename Particle<2>::PIsect &it) const
		{
			it.did_hit = p.dist(b.get_origin()) < p.m_r;
			return it.did_hit;
		}
	};

	/*******************************************************/
	// Point data - Long Beam query 3D
	template<int D>
	class PLB3D: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_vpoints_bvh, *m_caustic_vpoints_bvh;
	public:
		PLB3D():VolumetricRadianceEstimator(1)
		{}

		PLB3D(std::vector<Particle<D> > particles, unsigned int nb_samples) :
			VolumetricRadianceEstimator(nb_samples)
		{
			m_global_vpoints_bvh = new BVH<Particle<D>, D>(particles);
		}

		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<const Particle<D>*>::const_iterator it;
			typename Particle<D>::PIsect isec; 
			Real pp;
			
			m_global_vpoints_bvh->intersect_all(c, c.get_parameter(), found);
			
			for (it = found.begin(); it != found.end(); it++) {
				if (intersect(c, **it, isec)) {
					I += (*this)(c, **it, isec, pp);
				}
			}

			return I;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			Spectrum f, u_t, u_s, Tr;
			const Medium<D>* m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			typename Particle<D>::PIsect isec2;
			
			u_t = m->get_extinction(p.m_v->m_position);
			Tr =  m->get_transmittance(isec.tr0) -  m->get_transmittance(isec.tr1);
			Tr /= u_t;
			
			f = m->f(b.get_origin(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*Tr*K(b,p)*f;
		}

		Real K(const Ray<D> &b, const Particle<D> &p) const
		{
			return 3./(4.*M_PI*p.m_r*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<D> &b,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};
	
	template<>
	bool PLB3D<3>::intersect(const Ray<3> &b,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 sphere2ray_pos = b.get_origin() - p.m_v->m_position;

		Real A = 1;
		Real B = 2 * dot(sphere2ray_pos, b.get_direction());
		Real C = dot(sphere2ray_pos, sphere2ray_pos) -  p.m_r * p.m_r;

		Real Disc = B*B-4*A*C;
		if (Disc >= 0) {
			Real t;

			if (Disc == 0) {
				return false;
			} else {
				// Calculate both intersections
				Real Disc_sqrt = sqrt(Disc);
				Real t0 = (- B + Disc_sqrt) / 2 * A;
				Real t1 = (- B - Disc_sqrt) / 2 * A;
				
				bool b0 = t0 < b.get_parameter();
				bool b1 = t1 < b.get_parameter();
				if (b0 || b1) {
					it.tr0 = min(t0,t1);
					it.tr1 = max(t0,t1);
					it.did_hit = true;
					return true;
				}
			}
		}
		return false;
	};

	/*******************************************************/
	// Point data - Short Beam query 3D
	template<int D>
	class PSB3D: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_vpoints_bvh, *m_caustic_vpoints_bvh;
	public:
		PSB3D():VolumetricRadianceEstimator(1){}
		PSB3D(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_vpoints_bvh = new BVH<Particle<D>, D>(particles); }
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> > p_isec;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> >::iterator pit;
			std::vector<const Particle<D>*>::const_iterator it;
			 
			m_global_vpoints_bvh->intersect_all(c, c.get_parameter(), found);
			// Trace nsamples in the ray, REUSE EACH ONE it in all the point particles within ray
			for (it=found.begin(); it != found.end(); it++)
			{	
				typename Particle<D>::PIsect isec;
				if (intersect(c, **it, isec))
					p_isec.push_back(std::make_pair(*it, isec));
			}		
			Real pr, pp, dist;
			for (int i=0; i < m_nb_samples_ray; i++)
			{
				const Medium<D>* m = c.get_medium();
				dist=m->sample_mean_free_path(c,pr);
				Ray<D> b(c); 
				b.set_parameter(dist); 
				for (pit=p_isec.begin(); pit != p_isec.end(); pit++)
					I += (*this)(b, *((*pit).first), (*pit).second, pp);
			}
			return I/m_nb_samples_ray;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			Spectrum f, Tr;
			const Medium<D>* m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			
			if (b.get_parameter() < isec.tr0)
				return Spectrum(0.);
			else if (b.get_parameter() > isec.tr1)
				Tr = Spectrum(isec.tr1 - isec.tr0);
			else Tr = Spectrum(b.get_parameter() - isec.tr0);

			f = m->f(b.get_position(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*Tr*K(b,p)*f;
		}

		Real K(const Ray<D> &b, const Particle<D> &p) const
		{
			return 3./(4.*M_PI*p.m_r*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<D> &b,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};
	
	template<>
	bool PSB3D<3>::intersect(const Ray<3> &b,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 sphere2ray_pos = b.get_origin() - p.m_v->m_position;

		Real A = 1;
		Real B = 2 * dot(sphere2ray_pos, b.get_direction());
		Real C = dot(sphere2ray_pos, sphere2ray_pos) -  p.m_r * p.m_r;

		Real Disc = B*B-4*A*C;
		if (Disc >= 0)
		{
			Real t;

			if (Disc == 0)
				return false;
			else
			{
				// Calculate both intersections
				Real Disc_sqrt = sqrt(Disc);
				Real t0 = (- B + Disc_sqrt) / 2 * A;
				Real t1 = (- B - Disc_sqrt) / 2 * A;
				
				bool b0 = t0 < b.get_parameter();
				bool b1 = t1 < b.get_parameter();
				if( b0||b1 )
				{
					it.tr0 = min(t0,t1);
					it.tr1 = max(t0,t1);
					it.did_hit = true;
					return true;
				}
			}
		}
		return false;
	};
	
	/*******************************************************/
	// Point data - Long Beam query 2D
	template<int D>
	class PLB2D: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_vpoints_bvh, *m_caustic_vpoints_bvh;
	public:
		PLB2D():VolumetricRadianceEstimator(1){}
		PLB2D(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_vpoints_bvh = new BVH<Particle<D>, D>(particles); }
		
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<const Particle<D>*>::const_iterator it;
			typename Particle<D>::PIsect isec; 
			Real pp;
			
			m_global_vpoints_bvh->intersect_all(c, c.get_parameter(), found);
			
			for (it=found.begin(); it != found.end(); it++)
				if (intersect(c, **it, isec))
					I += (*this)(c, **it, isec, pp);

			return I;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			Spectrum f, Tr;
			const Medium<D>* m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			typename Particle<D>::PIsect isec2;
			
			Tr = m->get_transmittance(isec.tr0);
			
			f = m->f(b.get_origin(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*Tr*K(b,p)*f;
		}

		Real K(const Ray<D> &b, const Particle<D> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<D> &b,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};

	template<>
	bool PLB2D<3>::intersect(const Ray<3> &b,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 b2p = p.m_v->m_position - b.get_origin();
		Real d = dot(b2p, b.get_direction());
		if (d < 0 || d > b.get_parameter()) 
			return false;
		
		Real s2 = cross(b2p,  b.get_direction()).length2();//dot(b2p, b2p) - d*d;
		if (s2 > p.m_r*p.m_r) return false;
		it.did_hit=true;
		it.tr0 = d;
		return true;
	};

	/*******************************************************/
	// Point data - Long Beam query 2D
	template<int D>
	class PSB2D: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_vpoints_bvh, *m_caustic_vpoints_bvh;
	public:
		PSB2D():VolumetricRadianceEstimator(1){}
		PSB2D(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_vpoints_bvh = new BVH<Particle<D>, D>(particles); }
		
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> > p_isec;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> >::iterator pit;
			std::vector<const Particle<D>*>::const_iterator it;
			 
			m_global_vpoints_bvh->intersect_all(c, c.get_parameter(), found);
			// Trace nsamples in the ray, REUSE EACH ONE it in all the point particles within ray
			for (it=found.begin(); it != found.end(); it++)
			{	
				typename Particle<D>::PIsect isec;
				if (intersect(c, **it, isec))
					p_isec.push_back(std::make_pair(*it, isec));
			}
			Real pr, pp, dist;
			for (int i=0; i < m_nb_samples_ray; i++)
			{
				const Medium<D>* m = c.get_medium();
				dist=m->sample_mean_free_path(c,pr);
				Ray<D> b(c); 
				b.set_parameter(dist); 
				for (pit=p_isec.begin(); pit != p_isec.end(); pit++)
					I += (*this)(b, *((*pit).first), (*pit).second, pp);
			}
			return I/m_nb_samples_ray;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			Spectrum f, Tr;
			const Medium<D>* m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			typename Particle<D>::PIsect isec2;
			
			if (b.get_parameter() > isec.tr0)
				Tr = Spectrum(1.);//m->get_transmittance(isec.tr0);
			else return Spectrum(0.);

			f = m->f(b.get_position(), p.m_v->m_direction, -b.get_direction());
			return p.m_power*Tr*K(b,p)*f;
		}

		Real K(const Ray<D> &b, const Particle<D> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}
		
		bool intersect(const Ray<D> &b,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};

	template<>
	bool PSB2D<3>::intersect(const Ray<3> &b,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 b2p = p.m_v->m_position - b.get_origin();
		Real d = dot(b2p, b.get_direction());
		if (d < 0 || d > b.get_parameter()) 
			return false;
		
		Real s2 = cross(b2p,  b.get_direction()).length2();//dot(b2p, b2p) - d*d;
		if (s2 > p.m_r*p.m_r) return false;
		it.did_hit=true;
		it.tr0 = d;
		return true;
	};
#endif
	/*******************************************************/
	// Long Beam data - Long Beam query 2D(1)
	template<unsigned D, class Radiance, class RadianceAttenuation>
	class LBLB2D1: public VolumetricRadianceEstimator<D, Radiance, RadianceAttenuation>
	{
protected:
		using VRE = VolumetricRadianceEstimator<D, Radiance, RadianceAttenuation>;
		using RadianceSampleR = typename VRE::RadianceSampleR;
		using ParticleR = typename VRE::ParticleR;
#ifdef _USE_EMBREE_
protected:
		using ParticleAtomR = ParticleAtom<D, Radiance, RadianceAttenuation>;
protected:
		std::vector<ParticleR> m_particles;

		RTCScene b_scene;
		RTCDevice b_device;
		EmbreeParticleData<Radiance, RadianceAttenuation> m_embree_data;
#else // _USE_EMBREE_
		BVH<ParticleR, D> m_global_longbeams_bvh;
#endif // _USE_EMBREE_
	public:
		LBLB2D1() :
			VRE(1)
#ifdef _USE_EMBREE_
			, b_scene(nullptr), b_device(nullptr)
#endif // _USE_EMBREE_
		{}
		LBLB2D1(std::vector<ParticleR> particles, unsigned int nb_samples) :
			VRE(nb_samples)
#ifdef _USE_EMBREE_
			, m_particles(particles)
#else // _USE_EMBREE_
			, m_global_longbeams_bvh(particles)
#endif // _USE_EMBREE_
		{
#ifdef _USE_EMBREE_
			for (const ParticleR& p : m_particles) {
				p.atomize(m_embree_data.atoms);
			}

			//EMBREE ACCELERATION 
			/* create new Embree device */
			b_device = rtcNewDevice(NULL);
			error_handler(nullptr, rtcDeviceGetError(b_device));

			/* set error handler */
			rtcDeviceSetErrorFunction2(b_device, error_handler, nullptr);

			/* create scene */
			b_scene = rtcDeviceNewScene(b_device, RTC_SCENE_STATIC | RTC_SCENE_INCOHERENT, RTC_INTERSECT1 | RTC_INTERPOLATE);
			
			unsigned int rtc_beams = rtcNewUserGeometry(b_scene, m_embree_data.atoms.size());
			rtcSetUserData(b_scene, rtc_beams, (void *)&(m_embree_data));

			rtcSetBoundsFunction(b_scene, rtc_beams, (RTCBoundsFunc)&rtc_particle_bounds);
			rtcSetIntersectFunction(b_scene, rtc_beams, (RTCIntersectFunc)&rtc_particle_intersect);
			rtcSetOccludedFunction(b_scene, rtc_beams, (RTCOccludedFunc)&rtc_particle_occluded);
			rtcCommit(b_scene);
#endif // _USE_EMBREE_
		}
		
		virtual ~LBLB2D1()
		{}

		// Steady
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			Real pp;
#ifdef _USE_EMBREE_
			RTCRay rtc_ray;
			/* Ray parameters */
			rtc_ray.org[0] = (float)c.get_origin()[0];
			rtc_ray.org[1] = (float)c.get_origin()[1];
			rtc_ray.org[2] = (float)c.get_origin()[2];
			rtc_ray.dir[0] = (float)c.get_direction()[0];
			rtc_ray.dir[1] = (float)c.get_direction()[1];
			rtc_ray.dir[2] = (float)c.get_direction()[2];
			rtc_ray.tnear = 0.0f;
			rtc_ray.tfar = c.get_parameter();
			rtc_ray.time = 0; /* Time for motion blur */
			rtc_ray.mask = -1;
			rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
			rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

			m_embree_data.clear();

			/* Intersect ray with beams */
			unsigned lastPrimID = rtc_ray.primID;
			rtcIntersect(b_scene, rtc_ray);

			/* Calculate beam contribution */
			using HitType = typename std::map<unsigned int, typename ParticleR::PIsect>::value_type;
			for (HitType& hit : m_embree_data.isects) {
				typename ParticleR::PIsect& isec = hit.second;

				I += (*this)(c, isec.atom.m_p, isec, pp);
			}
#else // _USE_EMBREE_
			std::vector<const ParticleR*> found;
			m_global_longbeams_bvh.intersect_all(c, c.get_parameter(), found);
			
			typename ParticleR::PIsect isec;
			for (const ParticleR* p : found) {
				if (intersect(c, *p, isec)) {
					I += (*this)(c, *p, isec, pp);
				}
			}
#endif // _USE_EMBREE_
			return I;
		}

		// Transient
		virtual void operator()(const Ray<D> &c, std::vector<RadianceSampleR> &samples) const override
		{
#ifdef _USE_EMBREE_
			RTCRay rtc_ray;
			// Ray origin
			rtc_ray.org[0] = (float)c.get_origin()[0];
			rtc_ray.org[1] = (float)c.get_origin()[1];
			rtc_ray.org[2] = (float)c.get_origin()[2];
			rtc_ray.dir[0] = (float)c.get_direction()[0];
			rtc_ray.dir[1] = (float)c.get_direction()[1];
			rtc_ray.dir[2] = (float)c.get_direction()[2];
			rtc_ray.tnear = 0.0f;
			rtc_ray.tfar = c.get_parameter();
			rtc_ray.time = 0; // Time for motion blur
			rtc_ray.mask = -1;
			rtc_ray.geomID = RTC_INVALID_GEOMETRY_ID;
			rtc_ray.primID = RTC_INVALID_GEOMETRY_ID;

			m_embree_data.clear();

			/* Intersect ray with beams */
			unsigned lastPrimID = rtc_ray.primID;
			rtcIntersect(b_scene, rtc_ray);

			/* Calculate beam contribution */
			using HitType = typename std::map<unsigned int, typename ParticleR::PIsect>::value_type;
			for (HitType& hit : m_embree_data.isects) {
				typename ParticleR::PIsect& isec = hit.second;

				(*this)(c, isec.atom.m_p, isec, samples);
			}
#else // _USE_EMBREE_
			std::vector<const ParticleR*> found;
			m_global_longbeams_bvh.intersect_all(c, c.get_parameter(), found);
			
			typename ParticleR::PIsect isec;
			for (const ParticleR* p : found) {
				if (intersect(c, *p, isec)) {
					(*this)(c, *p, isec, samples);
				}
			}
#endif // _USE_EMBREE_
		}

		// Steady
		virtual Radiance operator()(const Ray<D> &c, const ParticleR &p, const typename ParticleR::PIsect &isec, Real &pdf) const override
		{ 
//			const Medium<D> *m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			const Medium<D> *m = VRE::m_world.medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), num, den, f;
			Real abs_c_theta = abs(isec.cos_theta);

			f = m->f(c.get_origin() + isec.tr0*c.get_direction(), p.m_v->m_direction, -c.get_direction());
			pdf = 1.;

#if BEAM1D
			return u_s*f*p.m_power*exp(-u_t*isec.tr0)*exp(-u_t*isec.tp0) / (isec.cos_theta * 2 * p.m_r);
#else // BEAM1D
			num = exp(-u_t*(isec.tr0 - isec.tr1)*(abs_c_theta - 1.)) - 1.;
			den = exp(u_t*(isec.tr0 + isec.tp0))*u_t*(abs_c_theta - 1);

			return u_s*K(c, p)*p.m_power*f*num / den;
#endif // BEAM1D
		}

		// Transient
		virtual void operator()(const Ray<D> &c, const ParticleR &p, const typename ParticleR::PIsect &isec, std::vector<RadianceSampleR> &samples) const override
		{
//			const Medium<D> *m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			const Medium<D> *m = VRE::m_world.medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), num, den, f;
			Real abs_c_theta = abs(isec.cos_theta);
			Real time;
			Real tr0, tr1, tp0, tp1;
			tr0 = isec.tr0;
			tr1 = isec.tr1;
			tp0 = isec.tp0;
			tp1 = isec.tp1;

			Real tr, tp, trstep = (1. + isec.cos_theta)/VRE::m_time_exp;
			Real p_time = p.m_v->m_time;
			for (tr = tr0; tr < tr1; tr += trstep) {
				tp = tp0 + isec.cos_theta*(tr - tr0);
				num = exp(-u_t*(-trstep)*(abs_c_theta - 1.)) -1.;
				den = exp(u_t*(tr + tp))*u_t*(abs_c_theta - 1);
				f = m->f(c.get_origin() + tr*c.get_direction(), p.m_v->m_direction, -c.get_direction());

				Real total_tof = p_time + VRE::m_world->time_of_flight(tr + tp);
				Real step_tof = VRE::m_world->time_of_flight(tr);

				samples.push_back(RadianceSampleR(K(c,p)*p.m_power*f*num/den, total_tof, 0, step_tof));
			}
		}

		Real K(const Ray<D> &c, const ParticleR &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}

		bool intersect(const Ray<D> &c,  const ParticleR &p, typename ParticleR::PIsect &it) const;
	};

	template<class Radiance, class RadianceAttenuation>
	class LBLB2D1<3, Radiance, RadianceAttenuation>
	{
	protected:
		using VRE = VolumetricRadianceEstimator<3, Radiance, RadianceAttenuation>;
		using RadianceSampleR = typename VRE::RadianceSampleR;
		using ParticleR = typename VRE::ParticleR;
	public:
		bool intersect(const Ray<3> &c, const ParticleR &p, typename ParticleR::PIsect &it) const
		{
			Vector3 end = p.get_end();// p.m_v->m_position+p.m_dist*p.m_v->m_direction;
			Real t;

			// Intersect
			Vector3 rw, ro, bw, bo, u1, u2;
			Real k1, k2, u_A, u_B, u_C, discr;
			Real tr0, tr1, tb0, tb1;
			const Real eps = 1e-7;

			bo = p.get_init();
			bw = p.m_v->m_direction;
			ro = c.get_origin();
			rw = c.get_direction();

			k1 = dot(rw, bw);
			if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
				return false;
		
			k2 = dot(ro-bo,bw);
			u1 = ro - bo - k2*bw;
			u2 = rw - k1*bw;
			u_A = dot(u2, u2); //
			u_B = 2.*dot(u1, u2); //
			u_C = dot(u1, u1) - p.m_r*p.m_r;
			discr = u_B*u_B - 4.*u_A*u_C;
		
			if (discr < eps) {//there's no solution to the quadratic equation or there's only one solution => no intersection
				return false;
			} else {
				Real discr_sqrt = sqrt(discr);
				tr1 = (-u_B + discr_sqrt)/(2.*u_A);
				tr0 = (-u_B - discr_sqrt)/(2.*u_A);

				tb0 = tr0*k1 + k2;
				tb1 = tr1*k1 + k2;
				if (fabs(tr1 - tr0) < eps)
					return false;
			}
			//clamp caps
			Real c_theta = dot(rw, bw);
			if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
				|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
				|| (tb1 > p.m_dist && tb0 > p.m_dist))
				return false;


			if (tr0 < 0.0) {
				tb0 = tb0 - tr0*c_theta;
				tr0 = 0.0f;
			}

			if (tr1 > c.get_parameter()) {
				tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
				tr1 = c.get_parameter();
			}

			//CHECK BOUNDS OF BEAM DISTANCES
			if (tb0 < 0.0) {
				tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
				tb0 = 0.0f;
			} else if (tb0 > p.m_dist) {
				tr0 = tr0 - (tb0-p.m_dist)/c_theta;
				tb0 = p.m_dist;
			}

			if (tb1 < 0.0) {
				tr1 = tr1 - tb1/c_theta;
				tb1 = 0.0f;
			} else if (tb1 > p.m_dist) {
				tr1 = tr1 - (tb1-p.m_dist)/c_theta;
				tb1 = p.m_dist;
			}

			if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.)
				return false;

			it.tr0 = tr0;
			it.tr1 = tr1;
			it.tp0 = tb0;
			it.tp1 = tb1;
			it.cos_theta = c_theta;

			return true;
		}
	};

#if 0
	/*******************************************************/
	// Long Beam data - Long Beam query 1D(1)
	template <int D>
	class LBLB1D1: public VolumetricRadianceEstimator<2>
	{
		BVH<Particle<2>, 2> *m_global_longbeams_bvh, *m_caustic_longbeams_bvh;
	public:
		LBLB1D1():VolumetricRadianceEstimator(1){}
		LBLB1D1(std::vector<Particle<2> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_longbeams_bvh = new BVH<Particle<2>, 2>(particles); }
		
		// Radiance of all beams over a ray
		// Steady
		Spectrum operator()(const Ray<2> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<2>*> found;
			std::vector<const Particle<2>*>::const_iterator it;
			typename Particle<2>::PIsect isec; 
			Real pp;
			
			m_global_longbeams_bvh->intersect_all(c, c.get_parameter(), found);
			
			for (it=found.begin(); it != found.end(); it++)
				if (intersect(c, **it, isec))
					I += (*this)(c, **it, isec, pp);

			return I;
		}
		
		// Transient
		void operator()(const Ray<2> &c, std::list<RadianceSample> &samples) const
		{
			Spectrum I(0.);
			std::vector<const Particle<2>*> found;
			std::vector<const Particle<2>*>::const_iterator it;
			typename Particle<2>::PIsect isec; 
			
			m_global_longbeams_bvh->intersect_all(c, c.get_parameter(), found);
			
			for (it=found.begin(); it != found.end(); it++)
				if (intersect(c, **it, isec))
					(*this)(c, **it, isec, samples);
		}

		// Radiance of all beams over a point
		// Steady
		Spectrum operator()(const VectorN<2> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<2>*> found;
			std::vector<const Particle<2>*>::const_iterator it;
			typename Particle<2>::PIsect isec; 
			Real pp;
			
			m_global_longbeams_bvh->intersect_all(c, found);
			
			for (it=found.begin(); it != found.end(); it++)
				if (intersect(c, **it, isec))
					I += (*this)(c, **it, isec, pp);

			return I;
		}

		// Radiance of single beam over a ray
		// Steady
		Spectrum operator()(const Ray<2> &c, const Particle<2> &p, const typename Particle<2>::PIsect &isec, Real &pdf ) const
		{ 
			const Medium<2> *m = dynamic_cast<const Path<2>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), num, den, f;
			Real abs_c_theta = abs(isec.cos_theta);

			num = exp(-u_t*(isec.tr0 - isec.tr1)*(abs_c_theta - 1.)) -1.;
			den = exp(u_t*(isec.tr0+min(isec.tp0, isec.tp1)))*u_t*(abs_c_theta - 1);

			f = m->f(c.get_origin()+isec.tr0*c.get_direction(),p.m_v->m_direction,-c.get_direction());
			pdf = 1.;

			return K(c,p)*p.m_power*f*num/den;
		}

		// Transient
		void operator()(const Ray<2> &c, const Particle<2> &p, const typename Particle<2>::PIsect &isec, std::list<RadianceSample> &samples ) const
		{
			const Medium<2> *m = dynamic_cast<const Path<2>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), num, den, f;
			Real abs_c_theta = abs(isec.cos_theta);
			Real time;
			Real tr0, tr1, tp0, tp1;
			tr0 = isec.tr0;
			tr1 = isec.tr1;
			tp0 = isec.tp0;
			tp1 = isec.tp1;

#define ToF(d) (m_world->time_of_flight(d))
#define FoT(t) (1./(t))

			Real tr, tp, trstep = FoT(m_time_exp/(1+isec.cos_theta));
			Real p_time = p.m_v->m_time;
			for (tr=tr0; tr<tr1; tr+=trstep)
			{
				tp = tp0 + isec.cos_theta*(tr - tr0);
				num = exp(-u_t*(-trstep)*(abs_c_theta - 1.)) -1.;
				den = exp(u_t*(tr+tp))*u_t*(abs_c_theta - 1);	
				f = m->f(c.get_origin()+tr*c.get_direction(),p.m_v->m_direction,-c.get_direction());
				samples.push_back(RadianceSample(K(c,p)*p.m_power*f*num/den, p_time + ToF(tr+tp), 0, ToF(tr)));
			}
		}

		// Radiance of a single beam over a point
		// Steady
		Spectrum operator()(const VectorN<2> &c, const Particle<2> &p, const typename Particle<2>::PIsect &isec, Real &pdf ) const
		{ 
			const Medium<2> *m = dynamic_cast<const Path<2>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), num, den, f;
			return K(c,p)*p.m_power*exp(-u_t*isec.tp0)*u_s/u_t;
		}

		Real K(const Ray<2> &c, const Particle<2> &p) const
		{
			return 1./(2*p.m_r);
		}

		Real K(const VectorN<2> &c, const Particle<2> &p) const
		{
			return 1./(2*p.m_r);
		}

		bool intersect(const Ray<2> &c,  const Particle<2> &p, typename Particle<2>::PIsect &it) const
		{
			Vector2 end = p.m_v->m_position+p.m_dist*p.m_v->m_direction; 
			Real t;
		
			Real c_theta = dot(p.m_v->m_direction, c.get_direction()), abs_c_theta;
			Real rdx, rdy, bdx, bdy; //direction
			Real rox, roy, box, boy; //origin
			Real den, tr, tr0, tr1, tb, tb0, tb1, h;

			abs_c_theta = fabs(c_theta);
	
			rdx = c.get_direction()[0]; 
			rdy = c.get_direction()[1];
			rox = c.get_origin()[0]; 
			roy = c.get_origin()[1];
	
			bdx = p.m_v->m_direction[0];
			bdy = p.m_v->m_direction[1];
			box = p.m_v->m_position[0];
			boy = p.m_v->m_position[1];
	
			den = bdy*rdx - bdx*rdy;
	
			if (fabs(den) < 1e-5 || (1.-abs_c_theta*abs_c_theta) < 1e-5)
				return false; //no intersection
	
			//intersection between beam ray and ray
			tr = (roy*bdx - boy*bdx + box*bdy - rox*bdy)/den;
			//intersection with blur region
			h = p.m_r/sqrtf(1.-abs_c_theta*abs_c_theta);

			tr1 = tr+h;
			tr0 = tr-h;
	
			Vector2 rtr = c.get_origin() + tr*c.get_direction();
			//perpendicular projections of tr1 and tr0 over beam ray
			tb = sqrtf(dot(p.m_v->m_position-rtr, p.m_v->m_position-rtr));
			if (fabs(p.m_v->m_direction[0] + tb*p.m_v->m_direction[0] - rtr[0]) > 0.0001 
				|| fabs(p.m_v->m_direction[1] + tb*p.m_v->m_direction[1] - rtr[1]) > 0.0001)
				tb = -tb;
	
			tb0 = tb - h*c_theta;
			tb1 = tb + h*c_theta;

			//clamp caps
			if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
				|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
				|| (tb1 > p.m_dist && tb0 > p.m_dist))
				return false;


			if (tr0 < 0.0)
			{
				tb0 = tb0 - tr0*c_theta;
				tr0 = 0.0f;
			}
		
			if (tr1 > c.get_parameter())
			{
				tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
				tr1 = c.get_parameter();
			}
		
			//CHECK BOUNDS OF BEAM DISTANCES
			if (tb0 < 0.0)
			{	
				tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
				tb0 = 0.0f;
			}
			else if (tb0 > p.m_dist)
			{
				tr0 = tr0 - (tb0-p.m_dist)/c_theta;
				tb0 = p.m_dist;
			}
		
			if (tb1 < 0.0)
			{
				tr1 = tr1 - tb1/c_theta;
				tb1 = 0.0f;
			}
			else if (tb1 > p.m_dist)
			{
				tr1 = tr1 - (tb1-p.m_dist)/c_theta;
				tb1 = p.m_dist;
			}

			if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
				return false;

			it.tr0 = tr0; it.tr1 = tr1; it.tp0 = tb0; it.tp1 = tb1; it.cos_theta = c_theta;
			return true;
		}
		
		bool intersect(const VectorN<2> &c,  const Particle<2> &p, typename Particle<2>::PIsect &it) const
		{
			Vector2 p2c = c-p.m_v->m_position;
			it.tp0 = p.m_v->m_direction.dot(p2c);
			return (it.tp0 > 0.) && it.tp0 < p.m_dist && (p.m_v->m_direction.cross(p2c).length2() < p.m_r*p.m_r);
		}
	};

	/*******************************************************/
	//Short Beam data - Long Beam query
	template<int D>
	class SBLB2D1: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_longbeams_bvh, *m_caustic_longbeams_bvh;
	public:
		SBLB2D1():VolumetricRadianceEstimator(1){}
		SBLB2D1(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_longbeams_bvh = new BVH<Particle<D>, D>(particles); }
		
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<const Particle<D>*>::const_iterator it;
			typename Particle<D>::PIsect isec; 
			Real pp;
			
			m_global_longbeams_bvh->intersect_all(c, c.get_parameter(), found);
			
			for (it=found.begin(); it != found.end(); it++)
				if (intersect(c, **it, isec))
					I += (*this)(c, **it, isec, pp);

			return I;
		}

		Spectrum operator()(const Ray<D> &c, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			const Medium<D> *m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), u_s = m->get_scattering(p.m_v->m_position), Tr, f;
			Real abs_c_theta = abs(isec.cos_theta);

			Tr = (m->get_transmittance(isec.tr0) - m->get_transmittance(isec.tr1))/u_t;
			f = m->f(c.get_origin()+isec.tr0*c.get_direction(),p.m_v->m_direction,-c.get_direction());
			pdf = 1.;

			return K(c,p)*p.m_power*f*Tr;
		}

		Real K(const Ray<D> &c, const Particle<D> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}

		bool intersect(const Ray<D> &c,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};

	template<>
	bool SBLB2D1<3>::intersect(const Ray<3> &c,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 end = p.m_v->m_position+p.m_dist*p.m_v->m_direction; 
		Real t;
		
		//intersect
		Vector3 rw, ro, bw, bo, u1, u2;
		Real k1, k2, u_A, u_B, u_C, discr;
		Real tr0, tr1, tb0, tb1;
		const Real eps = 1e-7;
	
		bo = p.m_v->m_position;
		bw = p.m_v->m_direction;
		ro = c.get_origin();
		rw = c.get_direction();

		k1 = dot(rw, bw);
		if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
			return false; 
	
		k2 = dot(ro-bo,bw);
		u1 = ro - bo - k2*bw;
		u2 = rw - k1*bw;
		u_A = dot(u2, u2); //
		u_B = 2.*dot(u1, u2); //
		u_C = dot(u1, u1) - p.m_r*p.m_r;
		discr = u_B*u_B - 4.*u_A*u_C;
	
		if (discr < eps) //there's no solution to the quadratic equation or there's only one solution => no intersection
			return false;
		else 
		{
			Real discr_sqrt = sqrt(discr);
			tr1 = (-u_B + discr_sqrt)/(2.*u_A); 
			tr0 = (-u_B - discr_sqrt)/(2.*u_A);
		
			tb0 = tr0*k1 + k2;
			tb1 = tr1*k1 + k2;
			if (fabs(tr1 - tr0) < eps) 
				return false;
		}
		//clamp caps
		Real c_theta = dot(rw, bw);
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
			|| (tb1 > p.m_dist && tb0 > p.m_dist))
			return false;

		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > c.get_parameter())
		{
			tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
			tr1 = c.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
			tb0 = 0.0f;
		}
		else if (tb0 > p.m_dist)
		{
			tr0 = tr0 - (tb0-p.m_dist)/c_theta;
			tb0 = p.m_dist;
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > p.m_dist)
		{
			tr1 = tr1 - (tb1-p.m_dist)/c_theta;
			tb1 = p.m_dist;
		}

		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
			return false;

		it.tr0 = tr0; it.tr1 = tr1; it.tp0 = tb0; it.tp1 = tb1; it.cos_theta = c_theta;
		return true;
	}

	/*******************************************************/
	// Long Beam data - Short Beam query 2D(1)
	template<int D>
	class LBSB2D1: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_longbeams_bvh, *m_caustic_longbeams_bvh;
	public:
		LBSB2D1():VolumetricRadianceEstimator(1){}
		LBSB2D1(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_longbeams_bvh = new BVH<Particle<D>, D>(particles); }
		
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> > p_isec;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> >::iterator pit;
			std::vector<const Particle<D>*>::const_iterator it;
			
			m_global_longbeams_bvh->intersect_all(c, c.get_parameter(), found);
			// Trace nsamples in the ray, REUSE EACH ONE it in all the point particles within ray
			for (it=found.begin(); it != found.end(); it++)
			{	
				typename Particle<D>::PIsect isec;
				if (intersect(c, **it, isec))
					p_isec.push_back(std::make_pair(*it, isec));
			}
			Real pr, pp, dist;
			for (int i=0; i < m_nb_samples_ray; i++)
			{
				const Medium<D>* m = c.get_medium();
				dist=m->sample_mean_free_path(c,pr);
				Ray<D> b(c); 
				b.set_parameter(dist); 
				for (pit=p_isec.begin(); pit != p_isec.end(); pit++)
					I += (*this)(b, *((*pit).first), (*pit).second, pp);
			}

			return I;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			const Medium<D> *m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), f, Tr;
			Real abs_c_theta = abs(isec.cos_theta);

			if (b.get_parameter() < isec.tr0)
				return Spectrum(0.);
			else if (b.get_parameter() > isec.tr1)
			{
				Tr = Spectrum(isec.tr1 - isec.tr0) *(m->get_transmittance(min(isec.tp1,isec.tp0)) - m->get_transmittance(max(isec.tp1,isec.tp0)))/(u_t);
			}
			else
			{	
				Tr = Spectrum(b.get_parameter() - isec.tr0); 
				
				Real tp = min(isec.tp0, isec.tp1) - abs_c_theta*(b.get_parameter() - isec.tr0);
				//Tr *= (isec.tp1>isec.tp0)?(m->get_transmittance(isec.tp0) - m->get_transmittance(isec.tp1)):(m->get_transmittance(isec.tp1) - m->get_transmittance(isec.tp0))/(u_t);
				Tr *= (m->get_transmittance(tp) - m->get_transmittance(min(isec.tp1,isec.tp0)))/(u_t);
			}

			f = m->f(b.get_origin()+isec.tr0*b.get_direction(),p.m_v->m_direction,-b.get_direction());
			pdf = 1.;

			return K(b,p)*Tr*p.m_power*f;
		}

		Real K(const Ray<D> &c, const Particle<D> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}

		bool intersect(const Ray<D> &c,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};

	template<>
	bool LBSB2D1<3>::intersect(const Ray<3> &c,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 end = p.m_v->m_position+p.m_dist*p.m_v->m_direction; 
		Real t;
		
		//intersect
		Vector3 rw, ro, bw, bo, u1, u2;
		Real k1, k2, u_A, u_B, u_C, discr;
		Real tr0, tr1, tb0, tb1;
		const Real eps = 1e-7;
	
		bo = p.m_v->m_position;
		bw = p.m_v->m_direction;
		ro = c.get_origin();
		rw = c.get_direction();

		k1 = dot(rw, bw);
		if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
			return false; 
	
		k2 = dot(ro-bo,bw);
		u1 = ro - bo - k2*bw;
		u2 = rw - k1*bw;
		u_A = dot(u2, u2); //
		u_B = 2.*dot(u1, u2); //
		u_C = dot(u1, u1) - p.m_r*p.m_r;
		discr = u_B*u_B - 4.*u_A*u_C;
	
		if (discr < eps) //there's no solution to the quadratic equation or there's only one solution => no intersection
			return false;
		else 
		{
			Real discr_sqrt = sqrt(discr);
			tr1 = (-u_B + discr_sqrt)/(2.*u_A); 
			tr0 = (-u_B - discr_sqrt)/(2.*u_A);
		
			tb0 = tr0*k1 + k2;
			tb1 = tr1*k1 + k2;
			if (fabs(tr1 - tr0) < eps) 
				return false;
		}
		//clamp caps
		Real c_theta = dot(rw, bw);
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
			|| (tb1 > p.m_dist && tb0 > p.m_dist))
			return false;

		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > c.get_parameter())
		{
			tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
			tr1 = c.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
			tb0 = 0.0f;
		}
		else if (tb0 > p.m_dist)
		{
			tr0 = tr0 - (tb0-p.m_dist)/c_theta;
			tb0 = p.m_dist;
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > p.m_dist)
		{
			tr1 = tr1 - (tb1-p.m_dist)/c_theta;
			tb1 = p.m_dist;
		}

		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
			return false;

		it.tr0 = tr0; it.tr1 = tr1; it.tp0 = tb0; it.tp1 = tb1; it.cos_theta = c_theta;
		return true;
	}


	/*******************************************************/
	// Short Beam data - Short Beam query 2D(1)
	template<int D>
	class SBSB2D1: public VolumetricRadianceEstimator<D>
	{
		BVH<Particle<D>, D> *m_global_longbeams_bvh, *m_caustic_longbeams_bvh;
	public:
		SBSB2D1():VolumetricRadianceEstimator(1){}
		SBSB2D1(std::vector<Particle<D> > particles, unsigned int nb_samples):VolumetricRadianceEstimator(nb_samples)
		{ m_global_longbeams_bvh = new BVH<Particle<D>, D>(particles); }
		
		Spectrum operator()(const Ray<D> &c) const
		{
			Spectrum I(0.);
			std::vector<const Particle<D>*> found;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> > p_isec;
			std::vector<std::pair<const Particle<D>*, typename Particle<D>::PIsect> >::iterator pit;
			std::vector<const Particle<D>*>::const_iterator it;
			
			m_global_longbeams_bvh->intersect_all(c, c.get_parameter(), found);
			// Trace nsamples in the ray, REUSE EACH ONE it in all the point particles within ray
			for (it=found.begin(); it != found.end(); it++)
			{	
				typename Particle<D>::PIsect isec;
				if (intersect(c, **it, isec))
					p_isec.push_back(std::make_pair(*it, isec));
			}
			Real pr, pp, dist;
			for (int i=0; i < m_nb_samples_ray; i++)
			{
				const Medium<D>* m = c.get_medium();
				dist=m->sample_mean_free_path(c,pr);
				Ray<D> b(c); 
				b.set_parameter(dist); 
				for (pit=p_isec.begin(); pit != p_isec.end(); pit++)
					I += (*this)(b, *((*pit).first), (*pit).second, pp);
			}

			return I;
		}

		Spectrum operator()(const Ray<D> &b, const Particle<D> &p, const typename Particle<D>::PIsect &isec, Real &pdf) const
		{ 
			const Medium<D> *m = dynamic_cast<const Path<D>::VolumeVertex*>(p.m_v)->m_medium;
			Spectrum u_t = m->get_extinction(p.m_v->m_position), f, Tr;
			Real abs_c_theta = abs(isec.cos_theta);

			if (b.get_parameter() < isec.tr0)
				return Spectrum(0.);
			else if (b.get_parameter() > isec.tr1)
				Tr = Spectrum(isec.tr1 - isec.tr0);
			else
				Tr = Spectrum(b.get_parameter() - isec.tr0); 

			f = m->f(b.get_origin()+isec.tr0*b.get_direction(),p.m_v->m_direction,-b.get_direction());
			pdf = 1.;

			return K(b,p)*Tr*p.m_power*f;
		}

		Real K(const Ray<D> &c, const Particle<D> &p) const
		{
			return 1./(M_PI*p.m_r*p.m_r);
		}

		bool intersect(const Ray<D> &c,  const Particle<D> &p, typename Particle<D>::PIsect &it) const;
	};

	template<>
	bool SBSB2D1<3>::intersect(const Ray<3> &c,  const Particle<3> &p, typename Particle<3>::PIsect &it) const
	{
		Vector3 end = p.m_v->m_position+p.m_dist*p.m_v->m_direction; 
		Real t;
		
		//intersect
		Vector3 rw, ro, bw, bo, u1, u2;
		Real k1, k2, u_A, u_B, u_C, discr;
		Real tr0, tr1, tb0, tb1;
		const Real eps = 1e-7;
	
		bo = p.m_v->m_position;
		bw = p.m_v->m_direction;
		ro = c.get_origin();
		rw = c.get_direction();

		k1 = dot(rw, bw);
		if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
			return false; 
	
		k2 = dot(ro-bo,bw);
		u1 = ro - bo - k2*bw;
		u2 = rw - k1*bw;
		u_A = dot(u2, u2); //
		u_B = 2.*dot(u1, u2); //
		u_C = dot(u1, u1) - p.m_r*p.m_r;
		discr = u_B*u_B - 4.*u_A*u_C;
	
		if (discr < eps) //there's no solution to the quadratic equation or there's only one solution => no intersection
			return false;
		else 
		{
			Real discr_sqrt = sqrt(discr);
			tr1 = (-u_B + discr_sqrt)/(2.*u_A); 
			tr0 = (-u_B - discr_sqrt)/(2.*u_A);
		
			tb0 = tr0*k1 + k2;
			tb1 = tr1*k1 + k2;
			if (fabs(tr1 - tr0) < eps) 
				return false;
		}
		//clamp caps
		Real c_theta = dot(rw, bw);
		if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
			|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
			|| (tb1 > p.m_dist && tb0 > p.m_dist))
			return false;

		if (tr0 < 0.0)
		{
			tb0 = tb0 - tr0*c_theta;
			tr0 = 0.0f;
		}
		
		if (tr1 > c.get_parameter())
		{
			tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
			tr1 = c.get_parameter();
		}
		
		//CHECK BOUNDS OF BEAM DISTANCES
		if (tb0 < 0.0)
		{	
			tr0 = tr0 - tb0/c_theta; //tb0 is negative, so it's actually increasing tr0
			tb0 = 0.0f;
		}
		else if (tb0 > p.m_dist)
		{
			tr0 = tr0 - (tb0-p.m_dist)/c_theta;
			tb0 = p.m_dist;
		}
		
		if (tb1 < 0.0)
		{
			tr1 = tr1 - tb1/c_theta;
			tb1 = 0.0f;
		}
		else if (tb1 > p.m_dist)
		{
			tr1 = tr1 - (tb1-p.m_dist)/c_theta;
			tb1 = p.m_dist;
		}

		if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) 
			return false;

		it.tr0 = tr0; it.tr1 = tr1; it.tp0 = tb0; it.tp1 = tb1; it.cos_theta = c_theta;
		return true;
	}
#endif
};
#endif
