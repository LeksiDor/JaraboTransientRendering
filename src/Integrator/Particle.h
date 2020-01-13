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

#ifndef _PARTICLE_H_
#define _PARTICLE_H_

#include "bunnykiller.h"

#include <map>

#include "LinearAlgebra/VectorN.h"
#include "RayTracing/Intersection.h"
#include "Integrator/Path.h"

#ifdef _USE_EMBREE_
#define EMBREE_STATIC_LIB
#include "External/Embree/include/embree2/rtcore_geometry_user.h"
#endif

/** Struct used to store the particles in the ParticleTracing class. 
	It stores all the information of the hit point, allowing  to know 
	the origin of the particle and its hit point, plus the energy and 
	the path travelled. 
	This allows using the particle for all particle-tracing-based algorithms,
	such as Photon Mapping [Jensen01], Instant Radiosity [Keller97], or 
	Photon Beams [Jarosz11].
	In order to include heterogeneous participating media for the Photon 
	Beams algorithm, some modifications must be performed to include the
	information of the travelled medium. */
template<unsigned D, class Radiance, class RadianceAttenuation>
class ParticleAtom;

template<unsigned D, class Radiance, class RadianceAttenuation>
class Particle
{
protected:
	using PathR = Path<D, Radiance, RadianceAttenuation>;
	using ParticleAtomR = ParticleAtom<D, Radiance, RadianceAttenuation>;
public:
	class PIsect{
	public:
		PIsect(ParticleAtomR& atom) :
			atom(atom), tr0(0.), tr1(0.), tp0(0.), tp1(0.), dist(0.), cos_theta(0.)
		{}
	public:
		ParticleAtomR& atom;
		Real tr0, tr1, tp0, tp1, dist;
		Real cos_theta;
	};
	enum ParticleType {
		NONE = 0,
		SurfPoint = 1,
		VolPoint = 2,
		ShortBeam = 4,
		LongBeam = 8
	};
	ParticleType m_type;
	unsigned char m_level;
	typename PathR::Vertex m_v;
	bool visited = false;
	Real m_r; //radius of particle
	Radiance m_power;
	Real m_dist; // beam range inside media
	Real m_init_dist;
	unsigned int m_id;
	
	Particle() :
		m_power(0.), m_level(0), m_v(NULL), m_type(SurfPoint)
	{}

	Particle(const typename PathR::SurfaceVertex *v, int level, ParticleType type, unsigned int id = 0) :
		m_v(new typename PathR::SurfaceVertex(*v)),
		m_level(level),
		m_type(type),
		m_init_dist(0.),
		m_power((_finite(v->get_vertex_value().avg()/v->get_subpath_pdf()))?(v->get_vertex_value())/v->get_subpath_pdf():(Radiance(0.))),
		m_id(id)
	{}

	Particle(const typename PathR::VolumeVertex *v, int level, ParticleType type, Real radius = 0.01, Real dist = 0., Real init_dist = 0., unsigned int id = 0)
		:m_v(new typename PathR::VolumeVertex(*v)),
		m_level(level),
		m_type(type),
		m_power((_finite(v->get_vertex_value().avg()/v->get_subpath_pdf()))?(v->get_vertex_value()/v->get_subpath_pdf()):(Radiance(0.))),
		m_r(radius),
		m_dist (dist),
		m_init_dist(init_dist),
		m_id(id)
		
	{}

	Real dist(const VectorN<D> &p) const
	{
		switch (m_type){
			case SurfPoint:
			case VolPoint:
				return (p - m_v->m_position).length();
			case LongBeam:
			case ShortBeam:
				return m_v->m_direction.cross(p - m_v->m_position).length();
			default:
				return 0.0f;
		}
	}

	unsigned int atomize(std::vector<ParticleAtomR>& atoms) const {
		switch (m_type) {
			case SurfPoint:
			case VolPoint:
			case ShortBeam:
				return 1;
			case LongBeam: {
				// TODO: Do this scene-dependent
				const Real alength = std::max(20.0f*m_r, 2.0f);// max(dot_abs(VectorN<D>(0, 1), m_v->m_direction)*m_dist, 3 * m_r);

				unsigned int k = 0;
				for (Real d = 0.0; d < m_dist; d += alength, k++) {
					Real init_dist = d;
					Real dist = std::min(m_dist - d, alength);

					atoms.push_back(ParticleAtomR(*this, init_dist, dist));
				}
				return k;
			}
			default:
				return 0.0f;
		}
	}

	AABB<D> get_bounding_box() const;

	bool intersect(Ray<D> &r, Intersection<D> &it, Real max_distance) const {
		return true;
		//return get_bounding_box().intersect(r, it, max_distance); 
	}

	bool intersect(const Ray<D> &r, Real max_distance) const { 
		return true;
		//return get_bounding_box().intersect(r, max_distance); 
	}

	bool intersect(const VectorN<D> &p) const { 
		return true;
		//return get_bounding_box().intersect(p);
	}

	Real get_intersection_cost() const {
		switch (m_type){
			case SurfPoint:
			case VolPoint:
				return 10.;
			case LongBeam:
			case ShortBeam:
				return 100.;
			default:
				return 0.0f;
		}
	}

	VectorN<D> get_init() const
	{
		return m_v->m_position + m_init_dist*m_v->m_direction;
	}

	VectorN<D> get_end() const
	{
		return m_v->m_position + (m_init_dist + m_dist)*m_v->m_direction;
	}

	VectorN<D> get_center() const
	{
		return m_v->m_position;
	}
}; //Particle

template<class Radiance, class RadianceAttenuation>
class Particle<3, Radiance, RadianceAttenuation>
{
	AABB<3> get_bounding_box() const
	{
		switch (Particle::m_type){
			case Particle::ParticleType::SurfPoint:
			case Particle::ParticleType::VolPoint:
				return AABB3D(Vector3(Particle::m_v->m_position[0] - Particle::m_r,
									  Particle::m_v->m_position[1] - Particle::m_r,
									  Particle::m_v->m_position[2] - Particle::m_r),
							  Vector3(Particle::m_v->m_position[0] + Particle::m_r,
									  Particle::m_v->m_position[1] + Particle::m_r,
									  Particle::m_v->m_position[2] + Particle::m_r));
			case Particle::ParticleType::LongBeam:
			case Particle::ParticleType::ShortBeam:
			default:
				Vector3 init = Particle::get_init();
				Vector3 end = Particle::get_end();
				return AABB3D(Vector3(std::min(init[0], end[0]) - Particle::m_r,
						              std::min(init[1], end[1]) - Particle::m_r,
									  std::min(init[2], end[2]) - Particle::m_r),
					          Vector3(std::max(init[0], end[0]) + Particle::m_r,
					        		  std::max(init[1], end[1]) + Particle::m_r,
									  std::max(init[2], end[2]) + Particle::m_r));
		}
	}
};

template<class Radiance, class RadianceAttenuation>
class Particle<2, Radiance, RadianceAttenuation>
{
	AABB<2> get_bounding_box() const
	{
		switch (Particle::m_type) {
			case Particle::ParticleType::SurfPoint:
			case Particle::ParticleType::VolPoint:
				return AABB2D(Vector2(Particle::m_v->m_position[0] - Particle::m_r,
									  Particle::m_v->m_position[1] - Particle::m_r),
							  Vector2(Particle::m_v->m_position[0] + Particle::m_r,
									  Particle::m_v->m_position[1] + Particle::m_r));
			case Particle::ParticleType::LongBeam:
			case Particle::ParticleType::ShortBeam:
			default:
				Vector2 end = Particle::m_v->m_position + Particle::m_dist*Particle::m_v->m_direction;
				return AABB2D(Vector2(std::min(Particle::m_v->m_position[0], end[0]) - Particle::m_r,
						              std::min(Particle::m_v->m_position[1], end[1]) - Particle::m_r),
							  Vector2(std::max(Particle::m_v->m_position[0], end[0]) + Particle::m_r,
									  std::max(Particle::m_v->m_position[1], end[1]) + Particle::m_r));
		}
	}
};

template<unsigned D, class Radiance, class RadianceAttenuation>
class ParticleAtom
{
protected:
	using ParticleR = Particle<D, Radiance, RadianceAttenuation>;
public:
	const ParticleR& m_p;
	Real m_init_dist;
	Real m_dist; // beam range inside media

public:
	ParticleAtom(const ParticleR& p, Real init_dist, Real dist) :
		m_p(p), m_init_dist(init_dist), m_dist(dist)
	{}

	AABB<D> get_bounding_box() const;

	VectorN<D> get_init() const
	{
		return m_p->m_v->m_position + m_init_dist * m_p->m_v->m_direction;
	}

	VectorN<D> get_end() const
	{
		return m_p->m_v->m_position + (m_init_dist + m_dist)*m_p->m_v->m_direction;
	}

	VectorN<D> get_center() const
	{
		return m_p->m_v->m_position;
	}
};

template<class Radiance, class RadianceAttenuation>
class ParticleAtom<3, Radiance, RadianceAttenuation>
{
	AABB<3> get_bounding_box() const
	{
		switch (ParticleAtom::m_p->m_type) {
			case ParticleAtom::ParticleR::SurfPoint:
			case ParticleAtom::ParticleR::VolPoint:
				return AABB3D(Vector3(ParticleAtom::m_p->m_v->m_position[0] - ParticleAtom::m_p->m_r,
									  ParticleAtom::m_p->m_v->m_position[1] - ParticleAtom::m_p->m_r,
									  ParticleAtom::m_p->m_v->m_position[2] - ParticleAtom::m_p->m_r),
							  Vector3(ParticleAtom::m_p->m_v->m_position[0] + ParticleAtom::m_p->m_r,
									  ParticleAtom::m_p->m_v->m_position[1] + ParticleAtom::m_p->m_r,
									  ParticleAtom::m_p->m_v->m_position[2] + ParticleAtom::m_p->m_r));
			case ParticleAtom::ParticleR::LongBeam:
			case ParticleAtom::ParticleR::ShortBeam:
			default:
				Vector3 init = ParticleAtom::get_init();
				Vector3 end = ParticleAtom::get_end();
				return AABB3D(Vector3(std::min(init[0], end[0]) - ParticleAtom::m_p->m_r,
						              std::min(init[1], end[1]) - ParticleAtom::m_p->m_r,
									  std::min(init[2], end[2]) - ParticleAtom::m_p->m_r),
					          Vector3(std::max(init[0], end[0]) + ParticleAtom::m_p->m_r,
					        		  std::max(init[1], end[1]) + ParticleAtom::m_p->m_r,
									  std::max(init[2], end[2]) + ParticleAtom::m_p->m_r));
		}
	}
};

template<class Radiance, class RadianceAttenuation>
class ParticleAtom<2, Radiance, RadianceAttenuation>
{
	AABB<2> get_bounding_box() const
	{
		switch (ParticleAtom::m_p->m_type) {
			case ParticleAtom::ParticleR::SurfPoint:
			case ParticleAtom::ParticleR::VolPoint:
				return AABB2D(Vector2(ParticleAtom::m_p->m_v->m_position[0] - ParticleAtom::m_p->m_r,
						              ParticleAtom::m_p->m_v->m_position[1] - ParticleAtom::m_p->m_r),
							  Vector2(ParticleAtom::m_p->m_v->m_position[0] + ParticleAtom::m_p->m_r,
									  ParticleAtom::m_p->m_v->m_position[1] + ParticleAtom::m_p->m_r));
			case ParticleAtom::ParticleR::LongBeam:
			case ParticleAtom::ParticleR::ShortBeam:
			default:
				Vector2 init = ParticleAtom::m_p->m_v->m_position;
				Vector2 end = init + ParticleAtom::m_dist*ParticleAtom::m_p->m_v->m_direction;
				return AABB2D(Vector2(std::min(init[0], end[0]) - ParticleAtom::m_p->m_r,
						              std::min(init[1], end[1]) - ParticleAtom::m_p->m_r),
					          Vector2(std::max(init[0], end[0]) + ParticleAtom::m_p->m_r,
					        		  std::max(init[1], end[1]) + ParticleAtom::m_p->m_r));
		}
	}
};

#ifdef _USE_EMBREE_

template<class Radiance, class RadianceAttenuation>
class EmbreeParticleData {
public:
	std::vector<ParticleAtom<3, Radiance, RadianceAttenuation>> atoms;
	std::map<unsigned int, typename Particle<3, Radiance, RadianceAttenuation>::PIsect> isects;
public:
	void clear()
	{
		atoms.clear();
		isects.clear();
	}
};

template<class Radiance, class RadianceAttenuation>
void rtc_particle_bounds(EmbreeParticleData<Radiance, RadianceAttenuation>* data, /*!< pointer to user data */
						size_t item,                                              /*!< item to calculate bounds for */
						RTCBounds& bounds_o                                       /*!< returns calculated bounds */)
{
	AABB<3> bb = data->atoms[item].get_bounding_box();
	bounds_o.lower_x = bb._min[0];
	bounds_o.lower_y = bb._min[1];
	bounds_o.lower_z = bb._min[2];

	bounds_o.upper_x = bb._max[0];
	bounds_o.upper_y = bb._max[1];
	bounds_o.upper_z = bb._max[2];
}

template<class Radiance, class RadianceAttenuation>
void rtc_particle_occluded(EmbreeParticleData<Radiance, RadianceAttenuation>* data, /*!< pointer to user data */
	                       RTCRay& ray,                                             /*!< ray to intersect */
	                       size_t item                                              /*!< item to intersect */)
{
	Ray<3> c(Vector3(ray.org[0], ray.org[1], ray.org[2]),
		     Vector3(ray.dir[0], ray.dir[1], ray.dir[2]), ray.tfar);

	const Particle<3, Radiance, RadianceAttenuation>* p = data->atoms[item].m_p;
	Vector3 end = p->get_end();
	Real t;

	//intersect
	Vector3 rw, ro, bw, bo, u1, u2;
	Real k1, k2, u_A, u_B, u_C, discr;
	Real tr0, tr1, tb0, tb1;
	const Real eps = 1e-7;

	bo = p->get_init();
	bw = p->m_v->m_direction;
	ro = c.get_origin();
	rw = c.get_direction();

	k1 = dot(rw, bw);
	if ((1. - std::abs(k1)) < eps) { //Beam and ray are parallel
		ray.geomID = 0;
		return;
	}

	k2 = dot(ro - bo, bw);
	u1 = ro - bo - k2*bw;
	u2 = rw - k1*bw;
	u_A = dot(u2, u2);
	u_B = 2.*dot(u1, u2);
	u_C = dot(u1, u1) - p->m_r*p->m_r;
	discr = u_B*u_B - 4.*u_A*u_C;

	if (discr < eps) { //there's no solution to the quadratic equation or there's only one solution => no intersection
		ray.geomID = 0;
		return;
	} else {
		Real discr_sqrt = std::sqrt(discr);
		tr1 = (-u_B + discr_sqrt) / (2.*u_A);
		tr0 = (-u_B - discr_sqrt) / (2.*u_A);

		tb0 = tr0*k1 + k2;
		tb1 = tr1*k1 + k2;
		if (fabs(tr1 - tr0) < eps) {
			ray.geomID = 0;
			return;
		}
	}
	// clamp caps
	Real c_theta = dot(rw, bw);
	if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
		|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
		|| (tb1 > p->m_dist && tb0 > p->m_dist)) {
		ray.geomID = 0;
		return;
	}

	if (tr0 < 0.0) {
		tb0 = tb0 - tr0*c_theta;
		tr0 = 0.0f;
	}

	if (tr1 > c.get_parameter()) {
		tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
		tr1 = c.get_parameter();
	}

	// CHECK BOUNDS OF BEAM DISTANCES
	if (tb0 < 0.0) {
		tr0 = tr0 - tb0 / c_theta; //tb0 is negative, so it's actually increasing tr0
		tb0 = 0.0f;
	} else if (tb0 > p->m_dist) {
		tr0 = tr0 - (tb0 - p->m_dist) / c_theta;
		tb0 = p->m_dist;
	}

	if (tb1 < 0.0) {
		tr1 = tr1 - tb1 / c_theta;
		tb1 = 0.0f;
	} else if (tb1 > p->m_dist) {
		tr1 = tr1 - (tb1 - p->m_dist) / c_theta;
		tb1 = p->m_dist;
	}

	if (tr0 < 0. || tb0 < 0. || tr1 < 0. || tb1 < 0.) {
		ray.geomID = 0;
		return;
	}
}

template<class Radiance, class RadianceAttenuation>
void rtc_particle_intersect(EmbreeParticleData<Radiance, RadianceAttenuation>* data, /*!< pointer to user data */
							RTCRay& ray,                                             /*!< ray to intersect */
							size_t item                                              /*!< item to intersect */)
{
	Ray<3> c(Vector3(ray.org[0], ray.org[1], ray.org[2]), 
			 Vector3(ray.dir[0], ray.dir[1], ray.dir[2]), ray.tfar);

	const ParticleAtom<3, Radiance, RadianceAttenuation> atom = data->atoms[item];
	const Particle<3, Radiance, RadianceAttenuation>* p = atom.m_p;
	Vector3 end = p->get_end();
	Real t;
	int i = 0;
	
	auto iter = data->isects.lower_bound(p->m_id);
	if (iter != data->isects.end() && !(data->isects.key_comp()(p->m_id, data->isects->first))) {
	    // Already intersected
		return;
	}

	// Intersect
	Vector3 rw, ro, bw, bo, u1, u2;
	Real k1, k2, u_A, u_B, u_C, discr;
	Real tr0, tr1, tb0, tb1;
	const Real eps = 1e-7;

	bo = p->get_init();
	bw = p->m_v->m_direction;
	ro = c.get_origin();
	rw = c.get_direction();

#define BEAM1D 1

#if BEAM1D
	Real minT1 = ray.tnear, minT2 = 0.0;
	Real maxT1 = ray.tfar, maxT2 = p->m_dist;
	Real oT1, oT2;
	const Vector3  d1d2c = cross(rw, bw);
	// Square of the sine between the two lines (||cross(rw, bw)|| = sinTheta).
	const Real sinThetaSqr = dot(d1d2c, d1d2c);

	// Slower code to test if the lines are too far apart.
	// oDistance = absDot((bo - ro), d1d2c) / d1d2c.size();
	// if(oDistance*oDistance >= maxDistSqr) return false;

	const Real ad = dot((bo - ro), d1d2c);
	// Lines too far apart.
	if (ad*ad >= (p->m_r*p->m_r)*sinThetaSqr)
		return ;
		
	// Cosine between the two lines.
	const Real d1d2 = dot(rw, bw);
	const Real d1d2Sqr = d1d2*d1d2;
	const Real d1d2SqrMinus1 = d1d2Sqr - 1.0f;

	// Parallel lines?
	if ((d1d2SqrMinus1 < 1e-5f && d1d2SqrMinus1 > -1e-5f))
		return ;

	const Real d1O1 = dot(rw, ro);
	const Real d1O2 = dot(rw, bo);

	oT1 = (d1O1 - d1O2 - d1d2 * (dot(bw, ro) - dot(bw, bo))) / d1d2SqrMinus1;

	// Out of range on ray 1.
	if ((oT1 <= minT1 || oT1 >= maxT1))
		return;

	oT2 = (oT1 + d1O1 - d1O2) / d1d2;
	// Out of range on ray 2.
	if (oT2 <= minT2 || oT2 >= maxT2 || !std::isfinite(oT2))
		return;

	const Real sinTheta = std::sqrt(sinThetaSqr);
	
	typename Particle<3, Radiance, RadianceAttenuation>::PIsect isec(atom);
	isec.tr0 = oT1;
	isec.tr1 = oT1;
	isec.tp0 = oT2;
	isec.tp1 = oT2;
	isec.cos_theta = std::sqrt(1 - sinThetaSqr);
	isec.dist = std::abs(ad) / sinTheta;
	isec.id = item;

	// Save intersection
	data->isects.insert(iter, std::make_pair(p->m_id, isec));

#else
	k1 = dot(rw, bw);

	if ((1. - fabs(k1)) < eps) //Beam and ray are parallel
		return;

	k2 = dot(ro - bo, bw);
	u1 = ro - bo - k2*bw;
	u2 = rw - k1*bw;
	u_A = dot(u2, u2); //
	u_B = 2.*dot(u1, u2); //
	u_C = dot(u1, u1) - p->m_r*p->m_r;
	discr = u_B*u_B - 4.*u_A*u_C;

	if (discr < eps) //there's no solution to the quadratic equation or there's only one solution => no intersection
		return;
	else
	{
		Real discr_sqrt = sqrt(discr);
		tr1 = (-u_B + discr_sqrt) / (2.*u_A);
		tr0 = (-u_B - discr_sqrt) / (2.*u_A);

		tb0 = tr0*k1 + k2;
		tb1 = tr1*k1 + k2;
		if (fabs(tr1 - tr0) < eps)
			return;
	}
	//clamp caps
	Real c_theta = dot(rw, bw);

	if ((tr0 < 0.0 && tr1 < 0.0) || (tb1 < 0.0 && tb0 < 0.0)
		|| (tr1 > c.get_parameter() && tr0 > c.get_parameter())
		|| (tb1 > p->m_dist && tb0 > p->m_dist))
		return;

	if (tr0 < 0.0)
	{
		tb0 = tb0 - tr0*c_theta;
		tr0 = 0;
	}

	if (tr1 > c.get_parameter())
	{
		tb1 = tb1 - (tr1 - c.get_parameter())*c_theta;
		tr1 = c.get_parameter();
	}

	//CHECK BOUNDS OF BEAM DISTANCES
	if (tb0 < 0.0)
	{
		tr0 = tr0 - tb0 / c_theta; //tb0 is negative, so it's actually increasing tr0
		tb0 = 0.0f;
	}
	else if (tb0 > p->m_dist)
	{
		tr0 = tr0 - (tb0 - p->m_dist) / c_theta;
		tb0 = p->m_dist;
	}

	if (tb1 < 0.0)
	{
		tr1 = tr1 - tb1 / c_theta;
		tb1 = 0.0f;
	}
	else if (tb1 > p->m_dist)
	{
		tr1 = tr1 - (tb1 - p->m_dist) / c_theta;
		tb1 = p->m_dist;
	}


	if (tr0 < 0.0 || tb0 < 0. || tr1 < 0. || tb1 < 0.)
		return;

	ray.tnear = tr0;
	ray.tfar = tr1;
	ray.instID = item;
	ray.primID = p->m_id;
	ray.geomID = item;
	ray.u = tb0;
	ray.v = tb1;
	ray.Ng[0] = c_theta;
	//it.tr0 = tr0; it.tr1 = tr1; it.tp0 = tb0; it.tp1 = tb1; it.cos_theta = c_theta;
#endif
}

template<class Radiance, class RadianceAttenuation>
void rtc_particle_filter(EmbreeParticleData<Radiance, RadianceAttenuation>* data, /*!< pointer to user data */
						 RTCRay& ray                                              /*!< intersection to filter */)
{
	// Reject all the beams intersections, all are transparent
	ray.geomID = RTC_INVALID_GEOMETRY_ID;
}
#endif

#endif //_PARTICLE_H_
