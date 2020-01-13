/*
 * Copyright (C) 2017, Ibon Guillen (http://giga.cps.unizar.es/~ibon/)
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

#ifndef _PROGRAMPARAMS_H_
#define _PROGRAMPARAMS_H_

#include "bunnykiller.h"

#include <cstdint>
#include <vector>

#include "LinearAlgebra/Vector3.h"
#include "LinearAlgebra/Vector2.h"
#include "Color/Spectrum.h"
#include "Integrator/Integrator.h"
#include "Utils/RandomNumbers.h"

struct ProgramParams {
	/* CAMERA and SENSOR PARAMETERS */
	VectorN<DIM> camera_position;
	VectorN<DIM> camera_looking_at;
	VectorN<DIM> camera_up;
	Real camera_view_plane_dist;
	Real camera_w;
	Real camera_h;

	/* VIRTUAL EMITTER PARAMETERS */
	bool emitter;
	VectorN<DIM> emitter_position;
	VectorN<DIM> emitter_looking_at;
	VectorN<DIM> emitter_up;
	Real emitter_view_plane_dist;
	bool emitter_othographic;
	Real emitter_w;
	Real emitter_h;

	/* MATERIAL */
	Spectrum grey;

	/* FILM */
	int width;
	int height;
	int time;
	int scanline;
	Real time_resolution;
	Real film_offset;
	Real aspectRatioY;
	int sqrt_spp;
	bool single_pixel;
	int single_pixel_x;
	int single_pixel_y;
	unsigned int nb_multibounces;

	/* FILTER */
	bool streak_kernel_based;
	int streak_nn;
	float streak_alpha;
	float streak_kernel_radius;

	/* FILM COMPONENTS */
	bool film_depth;
	bool film_positions;
	bool film_normals;

	std::vector<VectorN<DIM>> image_coordinates_vpl;

	/* FILES AND FOLDERS */
	char* image_file;
	char* log_file;

	bool transient;
	bool is_multibounce_streak;
	bool camera_unwarp;
	bool time_sampling;
	bool film_usereveal;

	/* MEDIUM PARAMETERS */
	Spectrum medium_sigma_a;
	Spectrum medium_sigma_s;
	Real hg_g;
	Scattering::Level scat_level;

	/* RNG */
	uint64_t new_seed;

	/* BDPT */
	unsigned int vol_path_tracing_samples;
	bool path_tracing;

	/* PPM */
	int max_nb_bounces;
	unsigned int pm_photon_shots;
	unsigned int pm_photons_nn;
	Real pm_kernel_shrinkage;

	ProgramParams(int dim) :
		camera_position(0., 0., 0.),
		camera_looking_at(1., 0., 0.),
		camera_up(0., 1., 0.),
		camera_view_plane_dist(1.),
		camera_w(0.),
		camera_h(0.),
		emitter(false),
		emitter_position(0., 0., 0.),
		emitter_looking_at(1., 0., 0.),
		emitter_up(0., 1., 0.),
		emitter_view_plane_dist(1.),
		emitter_othographic(false),
		emitter_w(1.),
		emitter_h(1.),
		grey(0.3f),
		width(256),
		time(512),
		scanline(-1),
		time_resolution(100.),
		film_offset(0.0),
		aspectRatioY(1.),
		sqrt_spp(2),
		single_pixel(false),
		single_pixel_x(-1),
		single_pixel_y(-1),
		nb_multibounces(3),
		streak_kernel_based(false),
		streak_nn(0),
		streak_alpha(4.0/7.0),
		streak_kernel_radius(1024.0),
		film_depth(false),
		film_positions(false),
		film_normals(false),
		image_coordinates_vpl(),
		image_file((char*)"name_file"),
		log_file(nullptr),
		transient(true),
		is_multibounce_streak(false),
		camera_unwarp(false),
		time_sampling(false),
		film_usereveal(false),
		medium_sigma_a(0.),
		medium_sigma_s(0.),
		hg_g(0.),
		scat_level(Scattering::ALL),
		new_seed(Random::RNG::default_state),
		vol_path_tracing_samples(1024),
		path_tracing(true),
		max_nb_bounces(20),
		pm_photon_shots(1000000),
		pm_photons_nn(25),
		pm_kernel_shrinkage(1.)
	{
		height = (dim == 2) ? 1 : width;
	}
}; /* ProgramParams */

#endif
