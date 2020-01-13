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

#ifndef _PHOTON_MAPPING_INTEGRATOR_H_
#define _PHOTON_MAPPING_INTEGRATOR_H_

//#define PHOTONS_REMOVE_SPURIOUS

#include "bunnykiller.h"

#include "DensityEstimation/DensityEstimationKernel.h"
#include "Integrator/BidirectionalPathTracing.h"
#include "DataStructures/KDTree.h"

template<unsigned D, class Radiance, class RadianceAttenuation>
class PhotonMapping : public BidirectionalPathTracing<D, Radiance, RadianceAttenuation>
{
protected:
	typedef BidirectionalPathTracing<D, Radiance, RadianceAttenuation> BDPTR;
	typedef typename BDPTR::RadianceSampleR RadianceSampleR;
	typedef typename BDPTR::RadSampleList RadSampleList;
	typedef typename BDPTR::FilmR FilmR;
protected:
	typedef Path<D, Radiance, RadianceAttenuation> PathR;
	typedef typename PathR::Vertex VertexR;
	typedef typename PathR::Vertex::Type VertexType;
protected:
	typedef DensityEsimationKernels::Box<2> Kernel2D;
	typedef DensityEsimationKernels::Box<3> Kernel3D;

	struct Photon
	{
		VectorN<D> m_position;
		VectorN<D> m_direction;

		Radiance m_R;
		Real m_p;
		Real m_pi;

		Real m_time, m_distance;
		unsigned int m_bounce; // NEED TO FIX THIS... Not crucial...

		Photon() :
			m_position(0.), m_direction(0.), m_R(0.), m_p(1.), m_pi(1.),
			m_time(0.), m_distance(0.), m_bounce(0)
		{}

		Photon(const VertexR &v, const Radiance &R) :
			m_position(v.get_vertex_position()), m_direction(v.get_vertex_direction()),
			m_R(v.get_vertex_value()*R), m_p(v.get_subpath_pdf()), m_pi(v.get_vertex_pdf()),
			m_time(v.get_subpath_delay()), m_distance(v.get_subpath_distance()), m_bounce(0)
		{}
	};

	struct PhotonInfo
	{
		Radiance f;
		float p;
		float t;
		Real maxChannel;
		PhotonInfo(Radiance fi, float pi, float ti, Real maxChanneli) :
			f(fi), p(pi), t(ti), maxChannel(maxChanneli)
		{}
		PhotonInfo() :
			f(0.), p(0.), t(0.), maxChannel(0.)
		{}
	};

	typedef KDTree<Photon, D> PhotonTree;
	typedef typename PhotonTree::Node PhotonTreeNode;
	typedef typename std::vector<const PhotonTreeNode*> VPhotons;
	typedef typename std::vector<const PhotonTreeNode*>::iterator ItPhoton;

	PhotonTree m_surface_tree, m_medium_tree;

	unsigned int m_nb_shots, m_nb_photons_de;
	Real m_max_bandwidth, m_bandwidth_reduction;

	Kernel2D m_kernel_surface;
	Kernel3D m_kernel_medium;

	void merge_vertices(const VertexR & veye, const Photon &photon, Real distance, VertexType type,
		Radiance &f, Real &p, Real &t)const;

	Radiance compute_radiance(const PathR &eye_path) const;
	void compute_radiance(const PathR &eye_path, RadSampleList &samples) const;

	mutable PhotonInfo* photonInfoSingleton;

	void print_info(FILE *f)const;
public:
	PhotonMapping(World<D, Radiance>& w, unsigned int nb_shots, unsigned int nb_photons_de, Real bandwidth_shrink = 1.,
			unsigned int incoming_samples = 1, FILE *_f_log = nullptr, FilmR *_film = nullptr) :
		BDPTR(w, incoming_samples, _f_log, _film), m_nb_shots(nb_shots), m_nb_photons_de(nb_photons_de),
		m_max_bandwidth(std::numeric_limits<Real>::infinity()), m_bandwidth_reduction(bandwidth_shrink),
		photonInfoSingleton(nullptr)
	{
#ifdef PHOTONS_REMOVE_SPURIOUS
		photonInfoSingleton = new PhotonInfo[nb_photons_de];
#endif
	}
	~PhotonMapping(){
#ifdef PHOTONS_REMOVE_SPURIOUS
		delete[] photonInfoSingleton;
#endif
	}

	/** Preprocess function for the integrator. */
	void preprocess();

	/**	Operator that integrates the radiance incoming from Ray r.
	This function includes the light incoming from the environment,
	and also the light scattered by the medium through the ray.	*/
	Radiance operator()(const Ray<D> &r)const;

	/**	Operator that integrates the radiance incoming to
	r.origin(), from the direction r.direction(). This is similar
	to operator(r), but taking into account different time samples.
	This more addecuate for transient rendering.	*/
	void operator()(const Ray<D> &r, RadSampleList &samples, const Real delta_time = 0., const unsigned int nb_time_samples = 1) const;
}; //PhotonMapping


template<unsigned D, class Radiance, class RadianceAttenuation>
void PhotonMapping<D, Radiance, RadianceAttenuation>::print_info(FILE *f_log) const
{
	if (f_log) {
		fprintf(f_log, "Photon Mapping\n");
		fprintf(f_log, "- Photon Shots: %d\n", m_nb_shots);
		fprintf(f_log, "- Photons @ Surfaces: %d\n", m_surface_tree.nb_elements());
		fprintf(f_log, "- Photons @ Media: %d\n", m_medium_tree.nb_elements());
		fprintf(f_log, "- Samples: %d\n", BDPTR::m_incoming_samples);
		fprintf(f_log, "- Photons Density Estimation: %d\n", m_nb_photons_de);
		fprintf(f_log, "- Kernel Bandwidth Reduction: %f\n", m_bandwidth_reduction);
		fprintf(f_log, "- Next Event: ");
		switch (this->m_next_event_sampling_technique) {
			case BDPTR::NextEventSampling::Time:
				fprintf(f_log, "Time Sampling - #Vertices %d\n", BDPTR::m_avg_nb_vertices_path);
				break;
			case BDPTR::NextEventSampling::MeanFreePath:
			default:
				fprintf(f_log, "Mean Free Path Sampling\n");
		}
		fprintf(f_log, "- Angular Sampling: ");
		switch (BDPTR::m_angular_sampling_technique) {
			case BDPTR::AngularSampling::PhaseFunction:
				fprintf(f_log, "Phase Function\n");
				break;
			case BDPTR::AngularSampling::ShortestTime:
				fprintf(f_log, "Shortest Time\n");
				break;
		}
		fprintf(f_log, "- Path Termination (Max Path Size %d): ", BDPTR::m_max_path_size);
		switch (BDPTR::m_termination_sampling_technique) {
			case BDPTR::TerminationSampling::TotalTime:
				fprintf(f_log, "Total Time\n");
				break;
			case BDPTR::TerminationSampling::RussianRoulette:
				fprintf(f_log, "Russian Roulette\n");
				break;
		}
		fprintf(f_log, "\n");
	}
}

template<unsigned D, class Radiance, class RadianceAttenuation>
void PhotonMapping<D,Radiance,RadianceAttenuation>::preprocess()
{
	// Preallocate a conservative amount of vertices
    unsigned int n_path = std::min(BDPTR::m_max_path_size, (unsigned int)(MAX_BPT_VERTICES));
	PathR light_path(n_path);

	Real m_inv_nb_shots = 1. / (Real)m_nb_shots;

	typename BDPTR::AngularSampling curr_angular_sampling = BDPTR::m_angular_sampling_technique;
	BDPTR::m_angular_sampling_technique = BDPTR::AngularSampling::PhaseFunction;

	m_surface_tree.clear();
	m_medium_tree.clear();

	for (size_t s = 0; s < m_nb_shots; ++s) {
		LightSample<D, Radiance> light_sample;
		Real pl1, pl2;
		BDPTR::m_world.sample_light(pl1)->sample(light_sample, pl2);

		BDPTR::generate_path(Ray<D>(light_sample.pos, light_sample.dir, false,
			BDPTR::m_world.get_ior(), BDPTR::m_world.get_medium()),
			light_sample.irradiance, pl1*pl2,
			light_sample.instant,
			VectorN<D>(0), TraceDirection::FROM_LIGHT, light_path);

		bool in_media = false;
		for (size_t il = 0; il < light_path.size(); ++il) {
			if (light_path[il].get_type() == VertexR::Type::MEDIUM) {
				in_media = true;
			}
		}
		
		in_media = true;
		if (in_media) {
			for (size_t il = 0; il < light_path.size(); ++il) {
				std::vector<Real> p(light_path[il].get_vertex_position().m_data,
					light_path[il].get_vertex_position().m_data + D);

				if (light_path[il].get_type() == VertexR::Type::SURFACE) {
					m_surface_tree.store(p, Photon(light_path[il], light_path.get_value() * (m_inv_nb_shots)));
				} else if (light_path[il].get_type() == VertexR::Type::MEDIUM) {
					m_medium_tree.store(p, Photon(light_path[il], light_path.get_value() * (m_inv_nb_shots)));
				}
			}
		}
	
		light_path.clear();
	}

	m_surface_tree.balance();
	m_medium_tree.balance();

	printf("# Photons: %d surface, %d media\n", m_surface_tree.nb_elements(), m_medium_tree.nb_elements());

	BDPTR::m_angular_sampling_technique = curr_angular_sampling;

	print_info(BDPTR::f_log);
}

// -----------------------------------------------------------------------------------------------
template<unsigned D, class Radiance, class RadianceAttenuation>
void PhotonMapping<D, Radiance, RadianceAttenuation>::merge_vertices(const VertexR &veye, const Photon &photon, Real bandwidth, VertexType type,
	Radiance &f, Real &p, Real &t) const
{
	RadianceAttenuation r;

	Real d = (veye.get_vertex_position() - photon.m_position).length();
	if (d > bandwidth) {
		f = Radiance(0.);
		t = veye.get_subpath_delay() + photon.m_time;
		p = 1.;

		return;
	}

	Real p1;
	veye.compute_photon_scattering(-photon.m_direction, r, t, p1);
	
	RadianceAttenuation fe(veye.get_vertex_value());
	set_attenuation_tracing_direction(fe, TraceDirection::FROM_LIGHT);
	set_attenuation_tracing_direction(r, TraceDirection::FROM_EYE);

	try {
		f = fe*(r * photon.m_R);
	} catch (...) {
		printf("Problem on Product...\n");
		photon.m_R.print();
		r.print();
		fe.print();
	}

	f *= (type == VertexR::Type::SURFACE) ?
		m_kernel_surface(d, bandwidth) :
		m_kernel_medium(d, bandwidth);
	p = p1 * veye.get_subpath_pdf() * photon.m_p;

	//printf("[%f, %f, %f] - %f - [%f - %f - %f]\n", f[0], f[1], f[2], p, p1, veye.get_subpath_pdf(), photon.m_p);
	t += veye.get_subpath_delay() + photon.m_time;
}

// -----------------------------------------------------------------------------------------------
template<unsigned D, class Radiance, class RadianceAttenuation>
Radiance PhotonMapping<D, Radiance, RadianceAttenuation>::compute_radiance(const PathR &eye_path) const
{
	Radiance radiance(0.);
	Radiance f;

	Real p, t;
	//Real w; // Not MIS for now...

	int nb_eye_vertices = std::min((unsigned int)(2), eye_path.size());
	// Iterate over all eye vertices
	for (size_t ie = 0; ie < 1; ++ie) {
		std::vector<Real> x(eye_path[ie].get_vertex_position().m_data,
						eye_path[ie].get_vertex_position().m_data + D);

		VPhotons photons;
		Real max_distance;

		// For now, ignore ray-traced direct illumination... TO FIX
		if (eye_path[ie].get_type() == VertexR::Type::SURFACE)
			m_surface_tree.find(x, m_nb_photons_de, photons, max_distance);
		else if (eye_path[ie].get_type() == VertexR::Type::MEDIUM)
			m_medium_tree.find(x, m_nb_photons_de, photons, max_distance);

		ItPhoton ph;
#ifdef PHOTONS_REMOVE_SPURIOUS
		float mean = 0;
		int cont = 0;
#endif
		for (ph = photons.begin(); ph != photons.end(); ++ph) {
			merge_vertices(eye_path[ie], (*ph)->data(), max_distance*m_bandwidth_reduction, eye_path[ie].get_type(), f, p, t);

#ifdef PHOTONS_REMOVE_SPURIOUS
			Real maxChannel = f.maxChannel();
			if (maxChannel > 0.0) {
				//cout << f.avg() << " ";
				photonInfoSingleton[cont] = PhotonInfo(f, p, t, maxChannel);
				mean += maxChannel;
				cont++;
			}
#else
			// NaN or infinity check. At worst, only a few energy will be lost, not the whole sample.
			if (f.is_valid()) {
				radiance += f / p;
			}
#endif
		}
		//cout << std::endl;

#ifdef PHOTONS_REMOVE_SPURIOUS
		if (cont > 0) {
			mean = mean / cont;
			//cout << mean <<std::endl;
			Real stdev = 0;
			for (int i = 0; i < cont; i++) {
				stdev += (photonInfoSingleton[i].maxChannel - mean)*(photonInfoSingleton[i].maxChannel - mean);
			}
			stdev = sqrt(stdev / cont);
			//cout << stdev <<" " << cont<<std::endl;
			//Inclusion
			for (int i = 0; i < cont; i++) {
				if (abs(photonInfoSingleton[i].maxChannel - mean) < 3 * stdev) {
					radiance += photonInfoSingleton[i].f / photonInfoSingleton[i].p;
				}
			}
		}
#endif
	}

	return radiance;
}

template<unsigned D, class Radiance, class RadianceAttenuation>
void PhotonMapping<D, Radiance, RadianceAttenuation>::compute_radiance(const PathR &eye_path, RadSampleList &samples) const
{
	Radiance radiance(0.);
	Radiance f;

	Real p, t;
	Real w; // Not MIS for now...

	int nb_eye_vertices = std::min((unsigned int)(2), eye_path.size());
	// Iterate over all eye vertices
	for (size_t ie = 0; ie < 1; ++ie) {
		std::vector<Real> x(eye_path[ie].get_vertex_position().m_data,
			eye_path[ie].get_vertex_position().m_data + D);

		VPhotons photons;
		Real max_distance;

		// For now, ignore ray-traced direct illumination... TO FIX
		if (eye_path[ie].get_type() == VertexR::Type::SURFACE)
			m_surface_tree.find(x, m_nb_photons_de, photons, max_distance);
		else if (eye_path[ie].get_type() == VertexR::Type::MEDIUM)
			m_medium_tree.find(x, m_nb_photons_de, photons, max_distance);

		ItPhoton ph;
#ifdef PHOTONS_REMOVE_SPURIOUS
		float mean = 0;
		int cont = 0;
#endif
		for (ph = photons.begin(); ph != photons.end(); ++ph) {
			merge_vertices(eye_path[ie], (*ph)->data(), max_distance*m_bandwidth_reduction, eye_path[ie].get_type(), f, p, t);
#ifdef PHOTONS_REMOVE_SPURIOUS
			Real maxChannel = f.maxChannel();
			if (maxChannel > 0) {
				photonInfoSingleton[cont] = PhotonInfo(f, p, t, maxChannel);
				mean += maxChannel;
				cont++;
			}
#else
			// NaN or infinity check
			if (f.is_valid()) {
				samples.push_back(RadianceSampleR(f * (BDPTR::m_inv_incoming_samples / p), t, ie + 1)); //TO FIX the bounce ID!
			}
#endif
		}

#ifdef PHOTONS_REMOVE_SPURIOUS
		if (cont > 0) {
			mean = mean / cont;
			Real stdev = 0;
			for (int i = 0; i < cont; i++) {
				stdev += (photonInfoSingleton[i].maxChannel - mean)*(photonInfoSingleton[i].maxChannel - mean);
			}
			stdev = sqrt(stdev / cont);

			// Inclusion
			for (int i = 0; i < cont; i++) {
				if (abs(photonInfoSingleton[i].maxChannel - mean) < 3 * stdev) {
					samples.push_back(RadianceSampleR(
						Radiance(photonInfoSingleton[i].f) * (m_inv_incoming_samples / photonInfoSingleton[i].p),
						photonInfoSingleton[i].t, ie + 1)); //TO FIX the bounce ID!
				}
			}
		}
#endif
	}
}

template<unsigned D, class Radiance, class RadianceAttenuation>
Radiance PhotonMapping<D, Radiance, RadianceAttenuation>::operator()(const Ray<D> &r) const
{
	Radiance reflected_radiance(0.);

	for (size_t sample = 0; sample < BDPTR::m_incoming_samples; ++sample) {
		PathR eye_path;

		// Sample light source. We do this first since it can be used
		// to guide the eye path when:
		// m_angular_sampling_technique = AngularSampling::ShortestTime.
		LightSample<D, Radiance> light_sample;
		Real pl1, pl2;
		BDPTR::m_world.sample_light(pl1)->sample(light_sample, pl2);

		// Generate eye path
		BDPTR::generate_path(r, Radiance(1.0), 1.0, 0.0, light_sample.pos, TraceDirection::FROM_EYE, eye_path);
		if (eye_path.size() == 0)
			continue;
				
		Radiance reye(compute_radiance(eye_path));
		
		try {
			reflected_radiance += reye;
		} catch (...) {
			printf("Ray: [%f %f %f]\n", r.get_direction()[0], r.get_direction()[1], r.get_direction()[2]);

			printf("Vertex Direction: [%f %f %f]\n", eye_path[0].get_vertex_direction()[0],
					eye_path[0].get_vertex_direction()[1], eye_path[0].get_vertex_direction()[2]);
		}
	}

	return reflected_radiance * BDPTR::m_inv_incoming_samples;
}

template<unsigned D, class Radiance, class RadianceAttenuation>
void PhotonMapping<D, Radiance, RadianceAttenuation>::operator()(const Ray<D> &r, RadSampleList &samples, const Real delta_time, const unsigned int nb_time_samples) const
{
	for (size_t sample = 0; sample < BDPTR::m_incoming_samples; ++sample) {
		PathR eye_path;

		// Sample light source. We do this first since it can be used
		// to guide the eye path when:
		// m_angular_sampling_technique = AngularSampling::ShortestTime.
		LightSample<D, Radiance> light_sample;
		Real pl1, pl2;
		BDPTR::m_world.sample_light(pl1)->sample(light_sample, pl2);

		// Generate eye path
		BDPTR::generate_path(r, Radiance(1), 1, 0, light_sample.pos, TraceDirection::FROM_EYE, eye_path);
		if (eye_path.size() == 0)
			continue;

		compute_radiance(eye_path, samples);
	}
}

#endif // _PHOTON_MAPPING_INTEGRATOR_H_
