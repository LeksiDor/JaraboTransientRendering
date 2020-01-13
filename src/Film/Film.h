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

#ifndef _FILM_H_
#define _FILM_H_

#include "bunnykiller.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "LinearAlgebra/Vector2.h"
#include "Filter/BoxFilter.h"
#include "Sampling/Sample.h"
#include "Integrator/RadianceSample.h"
#include "Color/Spectrum.h"
#include "Color/PolarizedSpectrum.h"
#include "Image/Image.h"
#include "Image/ImageIO.h"
#include "Film/FilmComponents.h"

/** A render Film for RGB images in Real format. */
template<unsigned D, class Radiance>
class Film
{
protected:
	using RadianceSampleR = RadianceSample<D, Radiance>;
	using RadianceSampleRecordR = RadianceSampleRecord<D, Radiance>;
	using RadianceSampleRecordVectorR = RadianceSampleRecordVector<D, Radiance>;
protected:
	/* The components stored by the film */
	FilmComponents m_components;

	/* Reconstruction filter */
	Filter* m_filter;

	/* Buffer containing the image */
	std::array<Imaging::Image<Real>, Radiance::components> m_image;
	Imaging::Image<Real> m_weight;

	/* Buffers with additional components */
	std::unique_ptr<Imaging::Image<Real>> m_depth_comp;
	std::unique_ptr<Imaging::Image<Real>> m_positions_comp;
	std::unique_ptr<Imaging::Image<Real>> m_normals_comp;

	/* Width and height of film */
	size_t width, height;
	size_t shift_x, shift_y;
	Real ratio;

	/* Half width and height of film */
	Real xw2, yw2;

	bool averaged;

	/* Film name */
	std::string m_name;
	std::string m_extension;

	void get_film_samples(const Vector2 &position, std::vector<Vector2>& samples) const
	{
		samples.clear();

		int size_kernel = m_filter->get_size();
		int half_size = size_kernel / 2;
		
		int x = (int) (position[0]);
		int y = (int) (position[1]);

		size_t init_x = std::max(0, x - half_size);
		size_t init_y = std::max(0, y - half_size);

		size_t end_x = std::min(x + half_size, int(width - 1));
		size_t end_y = std::min(y + half_size, int(height - 1));

		for (size_t j = init_y; j <= end_y; j++) {
			for (size_t i = init_x; i <= end_x; i++) {
				samples.push_back(Vector2(Real(i) + .5, Real(j) + .5));
			}
		}
	}

	void get_film_samples(const Real x, const unsigned int d, std::vector<Real>& samples) const
	{
		samples.clear();

		int size_kernel = m_filter->get_size();
		int half_size = size_kernel / 2;
		
		int x_int = (int) (std::floor(x));

		size_t init_x = std::max(0, x_int - half_size);
		size_t end_x = std::min(x_int + half_size, (d ? int(width - 1) : int(height - 1)));

		for (size_t i = init_x; i <= end_x; i++) {
			samples.push_back(Real(i) + .5);
		}
	}

	inline void draw_pixel(unsigned int x, unsigned int y, const Radiance& cv, Real w = -1.);

	template<class Record>
	inline void add_components(unsigned int x, unsigned int y, const Record& rec, Real w)
	{
		if (m_components && FilmComponents::DEPTH) {
			Real dis = rec.distance * w;
			m_depth_comp->add(x, y, &dis, 1);
		}

		if (m_components && FilmComponents::POSITIONS) {
			VectorN<D> pos = rec.pos * w;
			m_positions_comp->add(x, y, &pos[0], D);
		}

		if (m_components && FilmComponents::NORMALS) {
			VectorN<D> normal = rec.normal * w;
			m_normals_comp->add(x, y, &normal[0], D);
		}
	}

	void average_and_write(FilmComponents comp, const char* name = nullptr);
public:
	/* Constructs a Film with specified dimensions */
	Film(unsigned w, unsigned h, Filter* _filter = nullptr,
		FilmComponents comp = FilmComponents::RADIANCE) :
			m_components(comp),
			m_filter(_filter ? _filter : (new BoxFilter(1))),
			m_image(),
			m_weight(w, h, 1),
			m_depth_comp(),
			m_positions_comp(),
			m_normals_comp(),
			width(w),
			height(h),
			shift_x(0),
			shift_y(0),
			ratio(Real(w) / Real(h)),
			xw2(Real(w) / 2.0),
			yw2(Real(h) / 2.0),
			averaged(false),
			m_name("image"),
			m_extension("hdr")
	{
		m_image.fill(Imaging::Image<Real>(width, height, Radiance::spectral_samples));

		if (m_components && FilmComponents::DEPTH) {
			m_depth_comp = std::make_unique<Imaging::Image<Real>>(width, height, 1);
		}

		if (m_components && FilmComponents::POSITIONS) {
			m_positions_comp = std::make_unique<Imaging::Image<Real>>(width, height, D);
		}

		if (m_components && FilmComponents::NORMALS) {
			m_normals_comp = std::make_unique<Imaging::Image<Real>>(width, height, D);
		}
	}

	virtual ~Film()
	{
	}

	virtual Real get_time_length() const
	{
		return 0.;
	}

	virtual Real get_time_offset() const
	{
		return 0.;
	}

	virtual unsigned int get_time_resolution() const
	{
		return 0;
	}

	virtual Real get_exposure_time() const
	{
		return 0.;
	}

	/* Returns width of Film */
	int get_width() const
	{
		return width;
	}

	/* Returns height of Film */
	int get_height() const
	{
		return height;
	}

	/** 2D transform from window coordinates to image coordinates.
	 Image coordinates are in the [-0.5, 0.5] x [-0.5/ratio, 0.5/ratio] range */
	const Vector2 window_coords2image_coords(const Vector2& p) const
	{
		return Vector2((p[0] - xw2 + shift_x) / width, (p[1] - yw2 + shift_y) / height / ratio);
	}
	
	void set_shift(const unsigned int sx, const unsigned int sy)
	{
		shift_x = sx;
		shift_y = sy;
	}

	void set_aspect_ratio(Real aspectRatioY)
	{
		ratio *= aspectRatioY;
	}
	
	// Allows setting the name of the streak images
	virtual void set_name(const char *name, const char *ext)
	{
		m_name = std::string(name);
		m_extension = std::string(ext);
	}

	virtual void add_sample(const Sample& sample, const RadianceSampleRecordR& rec)
	{
		averaged = false;

		std::vector<Vector2> film_samples;
		get_film_samples(sample.position, film_samples);
		
		for (const Vector2& s : film_samples) {
			Real w = m_filter->evaluate(s - sample.position) * sample.weight;
			unsigned int sx = (unsigned int) (s[0]);
			unsigned int sy = (unsigned int) (s[1]);

			draw_pixel(sx, sy, rec.sample.radiance, w);

			add_components(sx, sy, rec, w);
		}
	}

	virtual void add_samples(const Sample& sample, const RadianceSampleRecordVectorR& rec)
	{
		averaged = false;

		std::vector<Vector2> film_samples;
		get_film_samples(sample.position, film_samples);

		for (const Vector2& s : film_samples) {
			Real w = m_filter->evaluate(s - sample.position) * sample.weight;
			unsigned int sx = (unsigned int) (s[0]);
			unsigned int sy = (unsigned int) (s[1]);
			
			m_weight.add(sy, sx, &w);

			// Now, we can iterate over all samples
			for (const RadianceSampleR& s : rec.samples) {
				draw_pixel(sx, sy, s.radiance * w, -1.);
			}

			add_components(sx, sy, rec, w);
		}
	}

	virtual void average(FilmComponents comp = FilmComponents::ALL)
	{
		Film::average_and_write(comp & m_components);

		averaged = true;
	}

	virtual void write(const char* name = nullptr)
	{
		name = (name) ? name : m_name.c_str();

		Film::average_and_write(m_components, name);

		averaged = true;
	}

	void clean()
	{
		for (size_t i = 0; i < Radiance::components; i++) {
			m_image[i].clean();
		}

		if (m_components && FilmComponents::DEPTH) {
			m_depth_comp->clean();
		}

		if (m_components && FilmComponents::POSITIONS) {
			m_positions_comp->clean();
		}

		if (m_components && FilmComponents::NORMALS) {
			m_normals_comp->clean();
		}
	}

public:
	Imaging::Image<Real>& image(unsigned int i = 0)
	{
		return m_image[i];
	}
}; /* class Film */

template<unsigned D, class Radiance>
void Film<D, Radiance>::draw_pixel(unsigned int x, unsigned int y, const Radiance& cv, Real w)
{
	if (w > 0.) {
		Radiance data = cv * w;
		m_image[0].add(x, y, &data[0], Radiance::spectral_samples);
		m_weight.add(x, y, &w, 1);
	} else {
		Radiance data = cv;
		m_image[0].add(x, y, &data[0], Radiance::spectral_samples);
	}
}

template<>
void Film<3, PolarizedLight<3>>::draw_pixel(unsigned int x, unsigned int y,
	const PolarizedLight<3>& cv, Real w)
{
	Real wp = 1.;
	if (w > 0.) {
		m_weight.add(x, y, &w, 1);
		wp = w;
	}

	for (size_t c = 0; c < PolarizedLight<3>::components; c++) {
		Spectrum data = cv[c] * wp;
		m_image[c].add(x, y, &data[0], PolarizedLight<3>::spectral_samples);
	}
}

template<unsigned D, class Radiance>
void Film<D, Radiance>:: average_and_write(FilmComponents comp, const char* name)
{
	/* Image information */
	Imaging::Image<Real>* out_image;
	unsigned channels = 0, rad_components = 0;

	/* Buffer to store the name */
	char nimg[2048];

	const char* ext;

	/* Buffer for the name for each sub-component */
	char nimg_comp[2048];

	while (true) {
		int nimg_size = 0;

		if (name)
			nimg_size += sprintf(nimg, "%s", name);

		/* Component-specific information */
		if (comp && FilmComponents::RADIANCE) {
			comp = comp ^ FilmComponents::RADIANCE;

			channels = Radiance::spectral_samples;
			rad_components = Radiance::components;

			ext = "hdr";

			out_image = &m_image[0];
		} else if (comp && FilmComponents::DEPTH) {
			comp = comp ^ FilmComponents::DEPTH;

			channels = 1;
			rad_components = 1;

			if (name)
				nimg_size += sprintf(nimg + nimg_size, "_depth");
			ext = "hdr";

			out_image = m_depth_comp.get();
		} else if (comp && FilmComponents::POSITIONS) {
			comp = comp ^ FilmComponents::POSITIONS;

			channels = D;
			rad_components = 1;

			if (name)
				nimg_size += sprintf(nimg + nimg_size, "_positions");
			ext = "raw";

			out_image = m_positions_comp.get();
		} else if (comp && FilmComponents::NORMALS) {
			comp = comp ^ FilmComponents::NORMALS;

			channels = D;
			rad_components = 1;

			if (name)
				nimg_size += sprintf(nimg + nimg_size, "_normals");
			ext = "raw";

			out_image = m_normals_comp.get();
		} else {
			/* If not an individual component, ignore */
			return;
		}

		/* Store all the sub-components */
		for (size_t c = 0; c < rad_components; c++) {
			if (!averaged) {
				/* Compensate by filter weight */
				out_image[c].weight(m_weight);

				if (c != 0) {
					/* Normalize value to 0..1, nonnegative */
					out_image[c].weight(out_image[0]);
					out_image[c].weight(2.0);
					out_image[c].add(0.5);
				}
			}

			if (name) {
				int nimg_comp_size = 0;

				if (rad_components > 1) {
					nimg_comp_size += sprintf(nimg_comp, "%s_s%d",
							nimg, unsigned(c));
				} else {
					nimg_comp_size += sprintf(nimg_comp, "%s", nimg);
				}

				/* Store all */
				if (channels > 3) {
					/*  Using one image per channel */
					for (size_t channel = 0; channel < channels; channel++) {
						Imaging::Image<Real> channel_image =
								out_image[c].channels(channel, channel);

						sprintf(nimg_comp + nimg_comp_size, "channel%d.%s",
							unsigned(channel), ext);

						Imaging::save(channel_image, nimg_comp);
					}
				} else {
					/* All RGB channels in a single image */
					sprintf(nimg_comp + nimg_comp_size, ".%s", ext);

					Imaging::save(out_image[c], nimg_comp);
				}
			}
		}
	}
}
#endif // _FILM_H_
