/*
 * Copyright (C) 2018, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _SAMPLES_STREAK_CAMERA_FILM_H_
#define _SAMPLES_STREAK_CAMERA_FILM_H_

#include "bunnykiller.h"

#include <cstdio>
#include <vector>

#include "Film/StreakCameraFilm.h"

template<unsigned D, class Radiance>
class SamplesStreakCameraFilm : public StreakCameraFilm<D, Radiance>
{
protected:
	using SCFTR = StreakCameraFilm<D, Radiance>;
	using Components = typename SCFTR::Components;
	using RadianceSampleR = typename SCFTR::RadianceSampleR;
protected:
	struct TimeSample
	{
		Real m_t;
		Radiance m_value;

		TimeSample(Real t, Radiance value) :
			m_t(t), m_value(value)
		{}
	};

	struct FilmSample
	{
		std::vector<TimeSample> samples;

		FilmSample()
		{}
	};

	std::vector< std::vector<FilmSample> > m_pixel_samples;
	unsigned int m_curr_pixel;
protected:
	void save_samples(unsigned int y, unsigned int x)
	{
		if (!m_pixel_samples[x].empty()) {
#if 0
			char name[1024];
			sprintf(name, "%s_%04d_%04d.csv", SCFTR::m_name.c_str(), y, x);

			FILE* sample_file = fopen(name, "w");

			fprintf(sample_file, "Iteration,Time");
			for (size_t i = 0; i < Radiance::spectral_samples; i++) {
				fprintf(sample_file, ",Channel %zu", i);
			}
			fprintf(sample_file, "\n");

			size_t n = 1;
			for (FilmSample& film_sample: m_pixel_samples[x]) {
				for (TimeSample& time_sample : film_sample.samples) {
					fprintf(sample_file, "%zu,%f", n, time_sample.m_t);
					for (size_t j = 0; j < Radiance::spectral_samples; j++) {
						fprintf(sample_file, ",%f", time_sample.m_value[j]);
					}
					fprintf(sample_file, "\n");
				}
				n++;
			}
#else
//			char name[1024];
//			sprintf(name, "%s_%04d_%04d.bin", SCFTR::m_name.c_str(), y, x);
//
//			FILE* sample_file = fopen(name, "w");
//
//			size_t n = 1;
//			for (FilmSample& film_sample: m_pixel_samples[x]) {
//				for (TimeSample& time_sample : film_sample.samples) {
//					fwrite(&n, sizeof(n), 1, sample_file);
//					fwrite(&time_sample.m_t, sizeof(time_sample.m_t), 1, sample_file);
//
//					for (size_t j = 0; j < Radiance::spectral_samples; j++) {
//						fwrite(&(time_sample.m_value[j]), sizeof(time_sample.m_value[j]), 1, sample_file);
//					}
//				}
//				n++;
//			}
#endif
			fflush(sample_file);
			fclose(sample_file);

			m_pixel_samples[x].clear();
		}
	}

	void init()
	{
		m_curr_pixel = 0;
		m_pixel_samples = std::vector< std::vector<FilmSample> >(SCFTR::width);
	}
public:
	SamplesStreakCameraFilm(int w, int h, int t, Real exposure, Filter *_filter,
			bool _camera_unwarp = false,	Components comp = SCFTR::RADIANCE) :
		SCFTR(w, h, t, exposure, _filter, _camera_unwarp, comp)
	{
		init();
	}

	virtual ~SamplesStreakCameraFilm()
	{
		clear();
	}

	virtual void add_samples(const Sample &sample, const std::vector<RadianceSampleR> &radiances)
		override
	{
		// First the film's samples are found...
		std::vector<Vector2> film_samples;
		std::vector<Real> temporal_samples;
		SCFTR::FTR::get_film_samples(sample.position, film_samples);

		// We iterate over all film samples
		for (const Vector2& s : film_samples) {
			unsigned int x = (unsigned int)(std::floor(s[0]));
			unsigned int y = (unsigned int)(std::floor(s[1]));

			// Check if new scanline is needed...
			while (y > SCFTR::m_last_available_slice) {
				// Save samples before discarding scanline
				save_samples(y, m_curr_pixel);
				m_curr_pixel = 0;

				for (std::vector<FilmSample>& samples: m_pixel_samples) {
					samples.clear();
				}

				// ... if so, get it!
				// First store the unused scanline...
				SCFTR::fix_and_store_slice(SCFTR::m_first_available_slice);

				// ... and remove it from memory.
				delete[] SCFTR::m_available_slices.front();
				SCFTR::m_available_slices.erase(SCFTR::m_available_slices.begin());

				// And then, add the new scanline
				SCFTR::add_new_scanline();
			}

			// Don't bother commenting this
			while (x > m_curr_pixel) {
				save_samples(y, m_curr_pixel);

				m_curr_pixel++;
			}

			// Else, get the current slice, and the temporal samples in the streak film
			int idx = y - SCFTR::m_first_available_slice;

			// We first get the weight of the pixel...
			// Filter in the y dimension...
			Real w = SCFTR::filter->evaluate(s - sample.position)*sample.weight;
			// ... and store it on weights matrix!
			SCFTR::weight[(y*SCFTR::width + x)] += w;

			// Add new FilmSample
			FilmSample& film_sample = *(m_pixel_samples[x].emplace(m_pixel_samples[x].end(), FilmSample()));

			// Now, we can iterate over all samples
			// Histogram!
			for (const RadianceSampleR& r : radiances) {
				// Get time_coordinate
				Real pos_sensor = ((SCFTR::camera_unwarp ? (r.time - r.distance) : r.time) -
								                  SCFTR::m_offset) / SCFTR::m_exposure_time;

				// If out of the range of the streak sensor, then discard.
				if (!(pos_sensor < SCFTR::m_time_resolution && pos_sensor >= 0.0))
					continue;

				SCFTR::get_film_samples(pos_sensor, 2, temporal_samples);

				// Finally, just draw the radiances on the streak film, multiplied by the
				// spatio-temporal filtering weight!. Note that the temporal weight doesn't act as
				// a filtering weight, but as a PSF. Thus, is additive.
				for (const Real& t : temporal_samples) {
					SCFTR::add_temporal_sample(r, w, x, idx, t, pos_sensor);

					film_sample.samples.push_back(TimeSample(t, r.radiance));
				}

				// And store it into the accumulated image
				SCFTR::draw_pixel(x, y, r.radiance*w, -1.);
			}
		}
	}

	void clear()
	{
		SCFTR::clear();
	}
}; // SamplesBasedStreakCameraFilm

#endif // _SAMPLES_STREAK_CAMERA_FILM_H_
