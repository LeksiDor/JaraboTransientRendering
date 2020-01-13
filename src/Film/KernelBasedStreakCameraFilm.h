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

#ifndef _KERNEL_BASED_STREAK_CAMERA_FILM_H_
#define _KERNEL_BASED_STREAK_CAMERA_FILM_H_

#include "bunnykiller.h"

#include <vector>

#include "Film/StreakCameraFilm.h"
#include "DataStructures/KDTree.h"
#include "DensityEstimation/DensityEstimationKernel.h"

template<unsigned D, class Radiance>
class KernelBasedStreakCameraFilm : public StreakCameraFilm<D, Radiance>
{
protected:
	using SCFTR = StreakCameraFilm<D, Radiance>;
	using RadianceSampleR = typename SCFTR::RadianceSampleR;
    using RadianceSampleRecordR = typename SCFTR::RadianceSampleRecordR;
	using RadianceSampleRecordVectorR = typename SCFTR::RadianceSampleRecordVectorR;
	using PhotonKDTree = KDTree<Radiance, 1>;
	using PhotonKDTreeNode = typename KDTree<Radiance, 1>::Node;
protected:
	// Fixed attributes
	unsigned int m_nearest_photons;
	bool m_filter_direct;
	unsigned int m_it;
	Real m_max_radius;
	Real m_alpha;
protected:
	// Attributes that vary each iteration
	struct IterationInfo {
		unsigned int m_it;
		Real m_radius_scale;

		IterationInfo(unsigned int it, Real radius_scale) :
			m_it(it), m_radius_scale(radius_scale)
		{}
	};
	std::vector<IterationInfo> m_iter_info;
protected:
	Radiance get_radiance(const Real ts, const PhotonKDTree &kd_tree,
	                      const int nearest_photons, const Real radius_scale,
						  const Real max_radius)
	{
		Radiance radiance(0.);
		DensityEsimationKernels::Perlin<1>  Kt;

		if (kd_tree.is_empty())
			return Radiance();

		if (nearest_photons > 0) {
			// Compute initial radius by NN search
			std::vector<const PhotonKDTreeNode*> photons;
			Real rad;

			std::vector<Real> v {ts};
			kd_tree.find(v, nearest_photons, photons, rad);

			// Update radius
			rad *= radius_scale;

			if (max_radius < rad)
				rad = max_radius*radius_scale;

			// And compute contribution
			for (const PhotonKDTreeNode* photon : photons) {
				Real tp = photon->point()[0];

				Real dist = std::abs(tp - ts);
				if (dist < rad) {
					radiance += photon->data() * Kt(dist, rad);
				}
			}
		} else {
			std::vector<const PhotonKDTreeNode*> photons;

			Real rad = max_radius*radius_scale;
			std::vector<Real> v {ts};
			kd_tree.find(v, rad, photons);

			// And compute contribution
			for (const PhotonKDTreeNode* photon : photons) {
				Real tp = photon->point()[0];

				Real dist = std::abs(tp - ts);
				radiance += photon->data() * Kt(dist, rad);
			}
		}

		return radiance;
	}
	
	Real compute_radius_scale(unsigned int it, Real alpha)
	{
		Real radius_scale = 1.;
		for (size_t i = 1; i <= it; i++) {
			radius_scale *= (Real(i - 1) + alpha) / Real(i);
		}
		return radius_scale;
	}

	Real next_radius_scale(Real radius_scale, unsigned int it, Real alpha)
	{
		return radius_scale * (Real(it) + alpha) / (Real(it) + 1.0);
	}
	
	void init()
	{
		m_iter_info = std::vector<IterationInfo>(SCFTR::width,
			IterationInfo(m_it, compute_radius_scale(m_it, m_alpha)));
	}
public:
	KernelBasedStreakCameraFilm(int w, int h, int t, Real exposure, Filter *_filter,
			Real kernel_radius, bool _camera_unwarp = false,
			Real alpha = Real(4.0/7.0), unsigned int it = 0,
			FilmComponents comp = FilmComponents::RADIANCE) :
		SCFTR(w, h, t, exposure, _camera_unwarp, _filter, comp),
			m_nearest_photons(0), m_filter_direct(false), m_it(it),
			m_max_radius(kernel_radius), m_alpha(alpha)
	{
		init();
	}

	KernelBasedStreakCameraFilm(int w, int h, int t, Real exposure, Filter *_filter,
			int nn_photons, bool _camera_unwarp = false,
			Real alpha = Real(4.0/7.0), Real max_radius = 1024.0,
			unsigned int it = 0, FilmComponents comp = FilmComponents::RADIANCE) :
		SCFTR(w, h, t, exposure, _camera_unwarp, _filter, comp),
			m_nearest_photons(nn_photons), m_filter_direct(false), m_it(it),
			m_max_radius(max_radius), m_alpha(alpha)
	{
		init();
	}

	virtual ~KernelBasedStreakCameraFilm()
	{
	}

	virtual void add_samples(const Sample &sample, const RadianceSampleRecordVectorR& rec) override
	{
//		// First the film's samples are found...
//		std::vector<Vector2> film_samples;
//		std::vector<Real> temporal_samples;
//		SCFTR::FTR::get_film_samples(sample.position, film_samples);
//
//		// We iterate over all film samples
//		for (const Vector2& s : film_samples) {
//			unsigned int x = (unsigned int)(s[0]);
//			unsigned int y = (unsigned int)(s[1]);
//
//			//Check if new scanline is needed...
//			while (y > SCFTR::m_last_available_slice) {
//				// ... if so, get it!
//				// First store the unused scanline...
//				SCFTR::fix_and_store_slice(SCFTR::m_first_available_slice);
//
//				// ... and remove it from memory.
////				delete[] SCFTR::m_available_slices.front();
//				SCFTR::m_available_slices.erase(SCFTR::m_available_slices.begin());
//
//				// And then, add the new scanline
//				SCFTR::add_new_scanline();
//
//				// Also, clean progressive info
//				m_iter_info = std::vector<IterationInfo>(SCFTR::width,
//					IterationInfo(m_it, compute_radius_scale(m_it, m_alpha)));
//			}
//
//			// Get the current slice, and the temporal samples in the streak film
//			int idx = y - SCFTR::m_first_available_slice;
//
//			// We first get the weight of the pixel...
//			// Filter in the y dimension...
//			Real w = SCFTR::m_filter->evaluate(s - sample.position)*sample.weight;
//			// ... and store it on weights matrix!
//			SCFTR::m_weight.add(x, y, &w);
//
//			// Get density estimation info
//			IterationInfo& info = m_iter_info[x];
//
//			// Acceleration structure for density estimation
//			PhotonKDTree photons;
//			Real min_t = Real(SCFTR::m_time_resolution);
//			Real max_t = Real(0.);
//
//			// Now, we can iterate over all samples
//			for (const RadianceSampleR& r : rec.samples) {
//				// Get time_coordinate
//				Real pos_sensor = ((SCFTR::camera_unwarp ? (r.time - rec.distance) : r.time) -
//				                  SCFTR::m_offset) / SCFTR::m_exposure_time;
//
//				// If out of the range of the streak sensor, then discard.
//				if (!(pos_sensor < SCFTR::m_time_resolution && pos_sensor >= 0.0))
//					continue;
//
//				min_t = std::min(pos_sensor, min_t);
//				max_t = std::max(pos_sensor, max_t);
//
//				if (r.bounce > 0 || m_filter_direct) {
//					// Store radiance
//					std::vector<Real> v {pos_sensor};
//					photons.store(v, r.radiance);
//				} else {
//					// If we do not apply the filtering to the direct component, store it as
//					// histogram
//					SCFTR::get_film_samples(pos_sensor, 2, temporal_samples);
//
//					// Finally, just draw the radiances on the streak film, multiplied by the
//					// spatio-temporal filtering weight!. Note that the temporal weight doesn't act
//					// as a filtering weight, but as a PSF. Thus, is additive.
//					for (const Real& t : temporal_samples) {
//						SCFTR::add_temporal_sample(r, w, x, idx, t, pos_sensor);
//					}
//				}
//				// Store radiance into the accumulated image
//				SCFTR::draw_pixel(x, y, r.radiance*w);
//			}
//
//			// Store depth if requested
//			SCFTR::FTR::add_components(x, y, rec, w);
//
//			// Prepare kd-tree for searching
//			photons.balance();
//
//			min_t = std::max(Real(0.), min_t -  m_max_radius);
//			max_t = std::min(Real(SCFTR::m_time_resolution - 1), max_t + m_max_radius);
//
//			// Now we apply density estimation at each pixel in the scanline
////			for (unsigned int t = 0; t < SCFTR::m_time_resolution; t++) {
//			for (unsigned int t = (unsigned int)(min_t); t <= (unsigned int)(max_t); t++) {
//				Real pos_sensor = Real(t) + Real(0.5);
//
//				Radiance r = get_radiance(pos_sensor, photons, m_nearest_photons,
//						info.m_radius_scale, m_max_radius);
//
//				if (!r.is_zero()) {
//					SCFTR::add_temporal_sample(r, w, x, idx, t, pos_sensor, false);
//				}
//			}
//
//			// Update density estimation info
//			info.m_it++;
//			info.m_radius_scale = next_radius_scale(info.m_radius_scale, info.m_it, m_alpha);
//		}
	}
}; // KernelBasedStreakCameraFilm

#endif // _STREAK_CAMERA_FILM_H_
