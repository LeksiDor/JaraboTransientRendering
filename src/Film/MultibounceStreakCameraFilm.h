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

#ifndef _MULTIBOUNCE_STREAK_CAMERA_FILM_H_
#define _MULTIBOUNCE_STREAK_CAMERA_FILM_H_

#include "bunnykiller.h"

#include <memory>
#include <vector>

#include "Film/Film.h"
#include "Film/StreakCameraFilm.h"
#include "Film/KernelBasedStreakCameraFilm.h"
#include "Filter/UberFilter.h"

template<unsigned D, class Radiance>
class MultibounceStreakFilm: public Film<D, Radiance>
{
protected:
    using FTR = Film<D, Radiance>;
    using RadianceSampleR = typename FTR::RadianceSampleR;
    using RadianceSampleRecordR = typename FTR::RadianceSampleRecordR;
	using RadianceSampleRecordVectorR = typename FTR::RadianceSampleRecordVectorR;
protected:
	typedef StreakCameraFilm<D, Radiance> StreakFilm;
	typedef KernelBasedStreakCameraFilm<D, Radiance> KernelBasedStreakFilm;

	Real m_exposure_time;

	unsigned int m_nb_bounces;
	std::vector<std::unique_ptr<StreakFilm>> m_bounce_films;

	void init_kernel_streaks(unsigned int _w, unsigned int _h, unsigned int _t, Real exposure,
			int nn_photons, Real kernel_radius, Real alpha, Filter *_filter,
			bool _camera_unwarp)
	{
		for (unsigned int b = 0; b < m_nb_bounces; ++b) {
			m_bounce_films.push_back(std::make_unique<KernelBasedStreakFilm>(_w, _h, _t, exposure,
					_filter, nn_photons, _camera_unwarp, alpha, kernel_radius));

			char name[1024];
			sprintf(name, "frame_b%02d", b+1);
			m_bounce_films.back()->set_name(name, "hdr");
		}
	}
public:
	MultibounceStreakFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
			unsigned int nb_bounces, Filter* filter = nullptr, bool camera_unwarp = false,
			FilmComponents comp = FilmComponents::RADIANCE) :
		FTR(w, h, filter, comp),
		m_exposure_time(exposure),
		m_nb_bounces(nb_bounces)
	{	
		for (unsigned int b = 0; b < m_nb_bounces; ++b) {
			m_bounce_films.emplace_back(std::make_unique<StreakFilm>(w, h, t, exposure, camera_unwarp, filter, comp));

			char name[1024];
			sprintf(name, "frame_b%02d", b+1);

			m_bounce_films.back()->set_name(name, "hdr");
		}
	}

	MultibounceStreakFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
			unsigned int nb_bounces, Filter *_filter, Real kernel_radius,
			Real alpha = Real(4.0/7.0), bool _camera_unwarp = false,
			FilmComponents comp = FilmComponents::RADIANCE) :
		FTR(w, h, _filter, comp), m_exposure_time(exposure), m_nb_bounces(nb_bounces)
	{
		init_kernel_streaks(w, h, t, exposure, 0, kernel_radius, alpha, _filter, _camera_unwarp);
	}

	MultibounceStreakFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
			unsigned int nb_bounces, Filter *_filter, int nn_photons, Real alpha = Real(4.0/7.0),
			Real kernel_radius = Real(1024.0), bool _camera_unwarp = false,
			FilmComponents comp = FilmComponents::RADIANCE) :
		FTR(w, h, _filter, comp), m_exposure_time(exposure), m_nb_bounces(nb_bounces)
	{
		init_kernel_streaks(w, h, t, exposure, nn_photons, kernel_radius, alpha, _filter, _camera_unwarp);
	}

	virtual ~MultibounceStreakFilm()
	{
	}
	
    virtual Real get_time_length() const override
    {
    	if (m_nb_bounces)
    		return m_bounce_films.front()->get_time_length();

    	return 0;
    }

    virtual Real get_time_offset() const override
    {
    	if (m_nb_bounces)
    		return m_bounce_films.front()->get_time_offset();

    	return 0;
    }
    
	// Returns time resolution of StreakFilm
    virtual unsigned int get_time_resolution() const override
    {
		if (m_nb_bounces)
			return m_bounce_films.front()->get_time_resolution();

		return 0;
	}

	// Allows setting the name of the streak images
    virtual void set_name(const char *name, const char *ext) override
	{
    	char nname[1024];

		for (unsigned int b = 0; b < m_nb_bounces; b++) {
			sprintf(nname, "%s_b%02d", name, b+1);

			m_bounce_films[b]->set_name(nname, ext);
		}

		FTR::set_name(name, ext);
	}

	virtual void add_sample(const Sample &sample, const RadianceSampleRecordR& rec) override
	{
		add_samples(sample,
				RadianceSampleRecordVectorR( { rec.sample }, rec.distance, rec.pos, rec.normal));
	}

	virtual void add_samples(const Sample &sample, const RadianceSampleRecordVectorR& rec) override
	{
		std::vector<RadianceSampleRecordVectorR> bounce_recs;
		bounce_recs.resize(m_nb_bounces, RadianceSampleRecordVectorR(rec.distance, rec.pos, rec.normal));

		for (const RadianceSampleR& s : rec.samples) {
			if (s.bounce < m_nb_bounces) {
				bounce_recs[s.bounce].add_sample(s);
			} else {
				bounce_recs.back().add_sample(s);
			}
		}

		for (size_t b = 0; b < m_bounce_films.size(); ++b) {
			if (bounce_recs[b].has_samples()) {
				m_bounce_films[b]->add_samples(sample, bounce_recs[b]);
			}
		}

		FTR::add_samples(sample, rec);
	}

	virtual Real get_exposure_time ()
	{
		return m_exposure_time;
	}

	void set_offset(Real offset)
	{
		for (size_t b = 0; b < m_nb_bounces; ++b) {
			m_bounce_films[b]->set_offset(offset);
		}
	}

	void set_streak(unsigned int y)
	{
		for (size_t b = 0; b < m_nb_bounces; ++b) {
			m_bounce_films[b]->set_streak(y);
		}
	}

public:
	Imaging::Image<Real>& bounce_streak_image(unsigned int b = 0, unsigned int y = 0, unsigned int i = 0)
	{
		return m_bounce_films[b]->streak_image(y);
	}
}; // MultibounceStreakFilm

#endif // _MULTIBOUNCE_STREAK_CAMERA_FILM_H_
