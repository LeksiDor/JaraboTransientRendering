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

#ifndef _PATH_H_
#define _PATH_H_

#include "bunnykiller.h"

#include <algorithm>

#include "LinearAlgebra/VectorN.h"
#include "Material/Material.h"
#include "Media/Medium.h"

template<unsigned D, class Radiance, class RadianceAttenuation>
class Path
{
public:
	class Vertex
	{
	public:
		enum Type {
			UNKNOWN = 0,
			SURFACE = 1,
			MEDIUM  = 2
		} m_type;
	protected:
		VectorN<D> m_position;
		VectorN<D> m_direction;

		RadianceAttenuation m_f;
		Real m_p;
		Real m_pi;

		Real m_time;
		Real m_distance;

        union VertexData {
			/* Data */
            struct SurfaceData {
                VectorN<D> m_normal;
                VectorN<D-1> m_uv;
                const Material<D>* m_material;
                const Medium<D>* m_medium;

                SurfaceData(VectorN<D> normal, VectorN<D-1> uv,
						const Material<D>* material, const Medium<D>* medium) :
					m_normal(normal),
					m_uv(uv),
					m_material(material),
					m_medium(medium)
				{}
            } surface;

            struct VolumeData {
                Spectrum m_sigma_s;
                Spectrum m_sigma_t;
                const Medium<D>* m_medium;

                VolumeData(VectorN<D> position, const Medium<D>* medium) :
                	m_sigma_s(medium->get_scattering(position)),
					m_sigma_t(medium->get_extinction(position)),
					m_medium(medium)
                {}
            } volume;

            /* Constructors */
            VertexData() {}

            VertexData(VectorN<D> normal, VectorN<D-1> uv, const Material<D>* material,
					const Medium<D>* medium) :
				surface(normal, uv, material, medium)
			{}

			VertexData(VectorN<D> position, const Medium<D>* medium) :
				volume(position, medium)
			{}

			~VertexData() {}
		} m_vertexData;
	public:
		Vertex() :
			m_type(Type::UNKNOWN), m_position(0.), m_direction(0.),
			m_f(0.), m_p(0.), m_pi(0.), m_time(0.), m_distance(0.),
			m_vertexData()
		{}

		Vertex(const Vertex &v) :
			m_type(v.m_type), m_position(v.m_position), m_direction(v.m_direction),
			m_f(v.m_f), m_p(v.m_p), m_pi(v.m_pi), m_time(v.m_time), m_distance(v.m_distance)
		{
			switch (m_type) {
				case SURFACE: {
					m_vertexData.surface = v.m_vertexData.surface;
				}
				break;
				case MEDIUM: {
					m_vertexData.volume = v.m_vertexData.volume;
				}
				break;
				default:
				break;
			}
		}

		Vertex& operator=(Vertex other)
		{
		    std::swap(*this, other);
		    return *this;
		}

		Vertex(const VectorN<D> &position, const VectorN<D> &direction, const VectorN<D> &normal,
				const VectorN<D-1> &uv, const Material<D> *material, const Medium<D> *medium,
				const RadianceAttenuation &f, Real p, Real pi, Real time = 0., Real distance = 0.) :
			m_type(Type::SURFACE),
			m_position(position),
			m_direction(direction),
			m_f(f),
			m_p(p),
			m_pi(pi),
			m_time(time),
			m_distance(distance),
			m_vertexData(normal, uv, material, medium)
		{}

		Vertex(const VectorN<D> &position, const VectorN<D> &direction, const Medium<D> *medium,
				const RadianceAttenuation &f, Real p, Real pi, Real time = 0., Real distance = 0.) :
			m_type(Type::MEDIUM),
			m_position(position),
			m_direction(direction),
			m_f(f),
			m_p(p),
			m_pi(pi),
			m_time(time),
			m_distance(distance),
			m_vertexData(position, medium)
		{}

		void compute_scattering(const VectorN<D> &dir1, RadianceAttenuation &f1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->f(m_direction, dir1, surf.m_normal, surf.m_uv, f1);

					if (surf.m_material->is_type(Reflectance::TWO_SIDED)) {
						f1 *= dot_abs(dir1, surf.m_normal);
					} else {
						f1 *= dot_clamped(dir1, surf.m_normal);
					}
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->f(m_position, m_direction, dir1, f1);

					f1 *= vol.m_sigma_s;
				}
				break;
				default:
				break;
			}
		}

		void compute_scattering(const VectorN<D> &dir1, RadianceAttenuation &f1, Real &t1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->f(m_direction, dir1, surf.m_normal, surf.m_uv, f1, t1, p1);

					if (surf.m_material->is_type(Reflectance::TWO_SIDED)) {
						f1 *= dot_abs(dir1, surf.m_normal);
					} else {
						f1 *= dot_clamped(dir1, surf.m_normal);
					}
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->f(m_position, m_direction, dir1, f1, t1, p1);

					f1 *= vol.m_sigma_s;
				}
				break;
				default:
				break;
			}
		}

		Spectrum compute_attenuation(const Ray<D> &r) const {
			switch (m_type) {
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					return vol.m_medium->get_transmittance(r);
				}
				break;
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;

					if (surf.m_medium) {
						return surf.m_medium->get_transmittance(r);
					}
				}
				default:
					return Spectrum(1.);
				break;
			}
		}

		Spectrum compute_attenuation(const VectorN<D> &p0, const VectorN<D> &p1) const {
			switch (m_type) {
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;

					VectorN<D> dir = p1 - p0;
					Real length = dir.length();
					dir /= length;

					const Ray<D> r(length, p0, dir);
					return vol.m_medium->get_transmittance(r);
				}
				break;
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					if (surf.m_medium) {
						VectorN<D> dir = p1 - p0;
						Real length = dir.length();
						dir /= length;

						const Ray<D> r(length, p0, dir);
						return surf.m_medium->get_transmittance(r);
					}
				}
				default:
					return Spectrum(1.);
				break;
			}
		}

		void compute_photon_scattering(const VectorN<D> &dir1, RadianceAttenuation &f1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->f(m_direction, dir1, surf.m_normal, surf.m_uv, f1);

					if (!surf.m_material->is_type(Reflectance::TWO_SIDED)
							&& dot(dir1, surf.m_normal) < 0) {
						f1 *= 0.;
					}
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->f(m_position, m_direction, dir1, f1);

					f1 *= vol.m_sigma_s;
				}
				break;
				default:
				break;
			}
		}

		void compute_photon_scattering(const VectorN<D> &dir1, RadianceAttenuation &f1, Real &t1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->f(m_direction, dir1, surf.m_normal, surf.m_uv, f1, t1, p1);

					if (!surf.m_material->is_type(Reflectance::TWO_SIDED)
							&& dot(dir1, surf.m_normal) < 0) {
						f1 *= 0.;
					}
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->f(m_position, m_direction, dir1, f1, t1, p1);
					f1 *= vol.m_sigma_s;
				}
				break;
				default:
				break;
			}
		}

		void compute_probability(const VectorN<D> &dir1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					p1 = surf.m_material->p(m_direction, dir1, surf.m_normal, surf.m_uv);
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					p1 = vol.m_medium->p(m_position, m_direction, dir1);
				}
				break;
				default:
				break;
			}
		}

		void compute_probability(const VectorN<D> &dir1, Real t1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					p1 = surf.m_material->p(m_direction, dir1, surf.m_normal, surf.m_uv, t1);
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					p1 = vol.m_medium->p(m_position, m_direction, dir1, t1);
				}
				break;
				default:
				break;
			}
		}

		void sample_direction(VectorN<D> &dir1, RadianceAttenuation &f1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->sample_direction(m_direction, dir1, surf.m_normal, surf.m_uv, f1, p1);
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->sample_direction(m_position, m_direction, dir1, f1, p1);
				}
				break;
				default:
				break;
			}
		}

		void sample_direction(VectorN<D> &dir1, RadianceAttenuation &f1, Real &t1, Real &p1) const {
			switch (m_type) {
				case SURFACE: {
					const typename VertexData::SurfaceData& surf = m_vertexData.surface;
					surf.m_material->sample_direction(m_direction, dir1, surf.m_normal, surf.m_uv, f1, t1, p1);
				}
				break;
				case MEDIUM: {
					const typename VertexData::VolumeData& vol = m_vertexData.volume;
					vol.m_medium->sample_direction(m_position, m_direction, dir1, f1, t1, p1);
				}
				break;
				default:
				break;
			}
		}

		inline const VectorN<D> &get_vertex_position() const
		{
			return m_position;
		}

		inline const VectorN<D> &get_vertex_direction() const
		{
			return m_direction;
		}

		inline const RadianceAttenuation &get_vertex_value() const
		{
			return m_f;
		}

		inline Real get_vertex_pdf() const
		{
			return m_pi;
		}

		inline Real get_subpath_pdf() const
		{
			return m_p;
		}

		inline Real get_subpath_distance() const
		{
			return m_distance;
		}

		inline Real get_subpath_delay() const
		{
			return m_time;
		}

		inline Real operator-(const Vertex &v1)
		{
			return (v1.m_position - Vertex::m_position).length();
		}

		inline Type get_type()const
		{
			return m_type;
		}
	}; // Vertex
private:
	std::vector<Vertex> m_vertices;
	Radiance m_value;

public:
	Path(Radiance value = Radiance(1.), size_t N = 32) :
		m_value(value)
	{
		// Allocate a minimum space
		m_vertices.reserve(N);
	}

	Path(size_t N) :
		m_value(Radiance(1.))
	{
		// Allocate a minimum space
		m_vertices.reserve(N);
	}

	~Path()
	{
		clear();
	}

	inline unsigned int size() const
	{
		return m_vertices.size();
	}
	
	inline const Vertex& operator[](unsigned int i) const
	{
		return m_vertices[i];
	}

	inline void add_surface_vertex(const VectorN<D> &position, const VectorN<D> &direction,
		const VectorN<D> &normal, const VectorN<D-1> &uv, const Material<D>* material, const Medium<D>* medium,
		const RadianceAttenuation &f, Real p, Real pi, Real time = 0., Real distance = 0.)
	{
		m_vertices.push_back(Vertex(position, direction, normal, uv, material, medium, f, p, pi, time, distance));
	}

	inline void add_media_vertex(const VectorN<D> &position, const VectorN<D> &direction, const Medium<D> *medium,
		const RadianceAttenuation &f, Real p, Real pi, Real time = 0., Real distance = 0.)
	{
		m_vertices.push_back(Vertex(position, direction, medium, f, p, pi, time, distance));
	}

	void clear()
	{
		m_vertices.clear();
	}

	void set_value(Radiance value = 1.)
	{
		m_value = value;
	}

	inline Radiance get_value() const
	{
		return m_value;
	}
}; //Path

#endif //_PATH_H_
