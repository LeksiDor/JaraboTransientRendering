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

#ifndef _STREAK_CAMERA_FILM_H_
#define _STREAK_CAMERA_FILM_H_

#include "bunnykiller.h"

#include <algorithm>
#include <array>
#include <memory>
#include <vector>

#include "Color/Spectrum.h"
#include "Film/Film.h"
#include "Filter/UberFilter.h"
#include "Image/Image.h"
#include "Image/ImageIO.h"

template<unsigned D, class Radiance>
class StreakCameraFilm : public Film<D, Radiance>
{
protected:
	using FTR = Film<D, Radiance>;
	using RadianceSampleR = typename FTR::RadianceSampleR;
	using RadianceSampleRecordR = typename FTR::RadianceSampleRecordR;
	using RadianceSampleRecordVectorR = typename FTR::RadianceSampleRecordVectorR;
protected:
	size_t m_time_resolution;
	Real m_exposure_time;
	Real m_offset;

	bool camera_unwarp;

	Filter* time_filter;
protected:
	/* Buffers containing the streak images */
	std::vector<std::array<Imaging::Image<Real>, Radiance::components>> m_available_slices;
	size_t m_first_available_slice;
	size_t m_last_available_slice;

	/* Buffers for optional variance estimation */
	std::vector<Imaging::Image<Real>> m_available_slices_variance;
	std::vector<Imaging::Image<Real>> m_available_slices_nsamples;
protected:
	void advance_scanline();

	void average_and_save_slice(const unsigned int y, const char* name = nullptr);

	void get_film_samples(const Real x, const unsigned int d, std::vector<Real>& samples) const
	{
		if (d < 2) {
			FTR::get_film_samples(x, d, samples);
		} else {
			samples.clear();

			int size_kernel = FTR::m_filter->get_size();
			int half_size = size_kernel / 2;

			int t = (int) (std::floor(x));

			unsigned int init_t = std::max<int>(0, t - half_size);
			unsigned int end_t = std::min<int>(t + half_size, int(m_time_resolution - 1));

			for (size_t i = init_t; i <= end_t; i++) {
				samples.push_back(Real(i) + .5);
			}
		}
	}

public:
	StreakCameraFilm(unsigned int w, unsigned int h, unsigned int t, Real exposure,
		bool camera_unwarp = false, Filter*_filter = nullptr,
		FilmComponents comp = FilmComponents::RADIANCE) :
			FTR(w, h, _filter, comp),
			m_time_resolution(t),
			m_exposure_time(exposure),
			m_offset(0.),
			camera_unwarp(camera_unwarp),
			time_filter(FTR::m_filter->get_subfilter(2)),
			m_first_available_slice(0),
			m_last_available_slice(0)
	{
		size_t available_slices = ((FTR::m_filter) ? FTR::m_filter->get_size() : 1);

		for (size_t i = 0; i < available_slices; i++) {

			m_available_slices.emplace_back(std::array<Imaging::Image<Real>,
					Radiance::components>());
			m_available_slices.back().fill(Imaging::Image<Real>(m_time_resolution,
					FTR::width, Radiance::spectral_samples));

			if (FTR::m_components && FilmComponents::VARIANCE) {
				m_available_slices_variance.emplace_back(Imaging::Image<Real>(
						m_time_resolution, FTR::width, Radiance::spectral_samples));
				m_available_slices_nsamples.emplace_back(Imaging::Image<Real>(
						m_time_resolution, FTR::width, Radiance::spectral_samples));
			}
		}

		m_last_available_slice = m_first_available_slice + available_slices - 1;
	}

	virtual ~StreakCameraFilm()
	{
	}

	// Returns time resolution of StreakFilm
	virtual unsigned int get_time_resolution() const override
	{
		return m_time_resolution;
	}

	virtual void add_sample(const Sample& sample, const RadianceSampleRecordR& rec) override
	{
		add_samples(sample,
				RadianceSampleRecordVectorR( { rec.sample }, rec.distance, rec.pos, rec.normal));
	}

	virtual void add_samples(const Sample& sample, const RadianceSampleRecordVectorR& rec) override
	{
		FTR::averaged = false;

		/* First the pixel coverage is calculated */
		std::vector<Vector2> film_samples;
		FTR::get_film_samples(sample.position, film_samples);

		/* We iterate over all pixels */
		std::vector<Real> temporal_samples;

		for (const Vector2& s : film_samples) {
			unsigned int x = (unsigned int) (std::floor(s[0]));
			unsigned int y = (unsigned int) (std::floor(s[1]));

			/* Check if new scanline is needed */
			while (y > m_last_available_slice) {
				average_and_save_slice(m_first_available_slice, FTR::m_name.c_str());
				advance_scanline();
			}

			/* Apply spatial filter */
			Real w = FTR::m_filter->evaluate(s - sample.position) * sample.weight;
			FTR::m_weight.add(x, y, &w);

			/* Iterate over samples */
			for (const RadianceSampleR& r : rec.samples) {
				/* Calculate sample's time coordinate */
				Real cam_time = camera_unwarp ? (r.time - rec.distance) : r.time;
				Real pos_sensor = (cam_time - m_offset)/ m_exposure_time;

				/* If out of the range of the sensor, discard */
				if (!(pos_sensor < m_time_resolution && pos_sensor >= 0.0))
					continue;

				/* Get positions on the temporal line */
				get_film_samples(pos_sensor, 2, temporal_samples);

				/* Finally, just draw the radiances on the streak film, multiplied by the
				 * spatio-temporal filtering weight!. Note that the temporal weight doesn't act as
				 * a filtering weight, but as a PSF. Thus, is additive.
				 */
				for (const Real& t : temporal_samples) {
					add_temporal_sample(r, w, x, y, t, pos_sensor);
				}

				/* And store it into the accumulated image */
				FTR::draw_pixel(x, y, r.radiance * w, -1.);
			}

			/* Store components if requested */
			FTR::add_components(x, y, rec, w);
		}
	}

	void add_temporal_sample(const RadianceSampleR& r, Real w, unsigned int x, unsigned int y,
		Real tr, Real ts, bool filtered = true);

	virtual Real get_time_length() const override
	{
		return Real(m_time_resolution) * m_exposure_time + m_offset;
	}

	virtual Real get_exposure_time() const override
	{
		return m_exposure_time;
	}

	void set_offset(Real offset)
	{
		m_offset = offset;
	}

	virtual Real get_time_offset() const override
	{
		return m_offset;
	}

	void set_streak(unsigned int y)
	{
		ptrdiff_t delta = m_last_available_slice - m_first_available_slice;
		m_first_available_slice = y;
		m_last_available_slice = y + delta;
	}

	virtual void average(FilmComponents comp = FilmComponents::ALL) override
	{
		for (size_t i = m_first_available_slice; i <= m_last_available_slice; i++) {
			average_and_save_slice(i);
		}

		FTR::average(comp);
	}

	virtual void write(const char* name) override
	{
		name = (name) ? name : FTR::m_name.c_str();

		for (size_t i = m_first_available_slice; i <= m_last_available_slice; i++) {
			average_and_save_slice(i, name);
		}

		FTR::write(name);
	}

public:
	Imaging::Image<Real>& streak_image(unsigned int y = 0, unsigned int i = 0)
	{
		size_t index = (size_t)std::floor(y) - m_first_available_slice;

		return m_available_slices[index][i];
	}
}; // StreakCameraFilm

template<unsigned D, class Radiance>
void StreakCameraFilm<D, Radiance>::advance_scanline()
{
	/* Clean the streak and move it to the back */
	for (size_t i = 0; i < Radiance::components; i++) {
		(*m_available_slices.begin())[i].clean();
	}

	std::rotate(m_available_slices.begin(),
		m_available_slices.begin()++, m_available_slices.end());

	/* If present, update the variance estimation too */
	if (FTR::m_components && FilmComponents::VARIANCE) {
		(*m_available_slices_variance.begin()).clean();
		std::rotate(m_available_slices_variance.begin(),
				m_available_slices_variance.begin()++, m_available_slices_variance.end());

		(*m_available_slices_nsamples.begin()).clean();
		std::rotate(m_available_slices_nsamples.begin(),
				m_available_slices_nsamples.begin()++, m_available_slices_nsamples.end());
	}

	m_first_available_slice = std::min<size_t>(m_first_available_slice + 1, FTR::height - 1);
	m_last_available_slice = std::min<size_t>(m_last_available_slice + 1, FTR::height - 1);
}

template<unsigned D, class Radiance>
void StreakCameraFilm<D, Radiance>::add_temporal_sample(const RadianceSampleR& r, Real w,
	unsigned int x, unsigned int y, Real tr, Real ts, bool filtered)
{
	size_t index = (size_t)std::floor(y) - m_first_available_slice;
	size_t t = (size_t)std::floor(ts);

	Radiance rad = r.radiance;

	if (w > 0.0) {
		rad *=  w;
	}

	if (filtered) {
		rad *= time_filter->evaluate(ts - tr);
	};

	m_available_slices[index][0].add(t, x, &rad[0], Radiance::spectral_samples);

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance */
		for (size_t c = 0; c < Radiance::spectral_samples; c++) {
			Real curr_val = m_available_slices[index][0](t, x, 1);
			Real& nsamples = m_available_slices_nsamples[index](t, x, 1);
			Real& var = m_available_slices_variance[index](t, x, 1);

			/* mean = current_value/nsamples */
			Real mean = (nsamples) ? curr_val/nsamples : Real(0.0);
			nsamples++;

			Real delta_rad = rad[c] - mean;
			mean = mean + delta_rad / nsamples;

			var += delta_rad * (rad[c] - mean);
		}
	}
}

template<>
void StreakCameraFilm<3, PolarizedLight<3>>::add_temporal_sample(const RadianceSampleR& r, Real w,
		unsigned int x, unsigned int y, Real tr, Real ts, bool filtered)
{
	size_t index = (size_t)std::floor(y) - m_first_available_slice;
	size_t t = (size_t)std::floor(ts);

	Real wt = 1.0;
	if (filtered) {
		wt = time_filter->evaluate(ts - tr);
	};

	for (size_t c = 0; c < PolarizedLight<3>::components; c++) {
		Spectrum rad = (w > 0.0) ? r.radiance[c] * w * wt : r.radiance[c] * wt;

		m_available_slices[index][c].add(t, x, &rad[0], PolarizedLight<3>::spectral_samples);
	}

	if (FTR::m_components && FilmComponents::VARIANCE) {
		/* Estimate online variance, over visible light */
		Spectrum rad = (w > 0.0) ? r.radiance[0] * w * wt : r.radiance[0] * wt;

		for (size_t c = 0; c < PolarizedLight<3>::spectral_samples; c++) {
			Real curr_val = m_available_slices[index][0](t, x, 1);
			Real& nsamples = m_available_slices_nsamples[index](t, x, 1);
			Real& var = m_available_slices_variance[index](t, x, 1);

			/* mean = current_value/nsamples */
			Real mean = (nsamples) ? curr_val/nsamples : Real(0.0);
			nsamples++;

			Real delta_rad = rad[c] - mean;
			mean = mean + delta_rad / nsamples;

			var += delta_rad * (rad[c] - mean);
		}
	}
}

template<unsigned D, class Radiance>
void StreakCameraFilm<D, Radiance>::average_and_save_slice(const unsigned int y, const char* name)
{
	size_t index = (size_t)std::floor(y) - m_first_available_slice;

	std::array<Imaging::Image<Real>, Radiance::components>& out_image = m_available_slices[index];

	Imaging::Image<Real>* out_image_var =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_available_slices_variance[index] : nullptr;
	Imaging::Image<Real>* out_image_ns =
			(FTR::m_components && FilmComponents::VARIANCE) ?
					&m_available_slices_nsamples[index] : nullptr;

	if (!FTR::averaged) {
		/* Apply spatial filter */
		for (size_t i = 0; i < FTR::width; i++) {
			Real w = FTR::m_weight(i, y, 1);

			for (size_t c = 0; c < Radiance::components; c++) {
				Imaging::ImageRow<> row = out_image[c].row(i);
				row.weight(w);

				if (c != 0) {
					/* Normalize value to 0..1, nonnegative */
					Imaging::ImageRow<> row_rad = out_image[0].row(i);

					row.weight(row_rad);
					row.weight(2.0);
					row.add(0.5);
				}
			}

			if (FTR::m_components && FilmComponents::VARIANCE) {
				Imaging::ImageRow<> row_var = out_image_var->row(i);
				row_var.weight(w);

				Imaging::ImageRow<> row_ns = out_image_ns->row(i);
				row_ns.weight(w);
			}
		}
	}

	/* Buffer to store the final name */
	char nimg[2048];

	/* Store streak */
	if (name) {
		for (size_t c = 0; c < Radiance::components; c++) {
			if (Radiance::components > 1) {
				sprintf(nimg, "%s_s%d_%04d.%s", name, unsigned(c), y, FTR::m_extension.c_str());
			} else {
				sprintf(nimg, "%s_%04d.%s", name, y, FTR::m_extension.c_str());
			}

			/* Workaround for hdr images less than 8 pixels wide (matlab fails reading them),
			 * store them transposed
			 */
			if (out_image[c].width() < 8) {
				out_image[c].transpose();
			}

			Imaging::save(out_image[c], nimg);
		}

		if (out_image_var) {
			sprintf(nimg, "%s_%04d_variance.%s", name, y, FTR::m_extension.c_str());
			if (out_image_var->width() < 8) {
				out_image_var->transpose();
			}
			Imaging::save(*out_image_var, nimg);
		}

		if (out_image_ns) {
			sprintf(nimg, "%s_%04d_nsamples.%s", name, y, FTR::m_extension.c_str());
			if (out_image_ns->width() < 8) {
				out_image_ns->transpose();
			}
			Imaging::save(*out_image_ns, nimg);
		}
	}
}

#endif // _STREAK_CAMERA_FILM_H_
