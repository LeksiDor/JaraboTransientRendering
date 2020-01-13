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

#ifndef _RENDERENGINE_H_
#define _RENDERENGINE_H_

#include "bunnykiller.h"

#include <cmath>

#include "Utils/Timer.h"
#include "LinearAlgebra/Vector3.h"
#include "Sampling/Sampler.h"
#include "Film/Film.h"
#include "RayTracing/World.h"
#include "Camera/Camera.h"
#include "Integrator/Integrator.h"
#include "Integrator/RadianceSample.h"

/** Base class implementing the simplest render engine. 
	It contains the rendering main loop, and supports 
	non-progressive, steady-state rendering.		*/
template<unsigned D, class Radiance, class RadianceAttenuation>
class RenderEngine
{
public:
	using FilmR = Film<D, Radiance>;
	using RadianceSampleR = RadianceSample<D, Radiance>;
	using RadianceSampleRecordR = RadianceSampleRecord<D, Radiance>;
	using RadianceSampleRecordVectorR = RadianceSampleRecordVector<D, Radiance>;
	using IntegratorR = Integrator<D, Radiance, RadianceAttenuation>;
	using WorldR = World<D, Radiance>;
protected:
	FILE* f_log;
public:
	/** Constructor. */
 	RenderEngine(FILE *_f_log) :
 		f_log(_f_log)
 	{}

 	virtual ~RenderEngine() {}
	
	/** Function that implements the rendering main loop.
		It gets as input the name of the image to be rendered,
		and all the different elements needed to render the 
		scene, including the scene itself, the integrator, and 
		the camera, film and sampler. */
	virtual void render(const char *name, WorldR& world, IntegratorR *integrator,
						const Camera<D>& camera, FilmR *film, Sampler *sampler) const;
}; // RenderEngine

template<unsigned D, class Radiance, class RadianceAttenuation>
void RenderEngine<D, Radiance, RadianceAttenuation>::render(const char *name, WorldR& world,
		IntegratorR *integrator, const Camera<D>& camera, FilmR *film, Sampler *sampler) const
{
	if (f_log) {
		fprintf(f_log, "RENDER\n");
		fprintf(f_log, "======\n");
	}

	// ----------------------------------------------------------------------
	// If no integrator has been defined, set error...	
	if (!integrator) {
		throw std::runtime_error("Error: Need to define an integrator for rendering\n");
	}

	// ----------------------------------------------------------------------
	// If needed, preprocess integrator ...
	printf("Preprocessing integrator...\r");
	Timer timer;
	timer.start();

	integrator->preprocess();
	{
		Real secs = timer.get_secs();
		unsigned hours = (unsigned)(secs)/3600;
		secs -= hours*3600;
		unsigned minutes = (unsigned)(secs)/60;
		secs -= minutes*60;

		printf("Preprocessed integrator:  [%d:%d:%f]\n", hours, minutes, secs);

		if (f_log) {
			fprintf(f_log, "Preprocessed integrator: %d:%d:%f\n", hours, minutes, secs);
		}
	}
	
	// ----------------------------------------------------------------------
	// Go rendering!
	printf("Rendering...\r");
	timer.start();

	// Main Loop: Loop over all samples in the image
	Sample film_sample;
	size_t nb_samples = sampler->get_nb_samples(),
		   curr_sample = 0;
	double tstamp = -1.;
	
	const PolarizationFrame<D> camera_frame = camera.get_frame();

	while (sampler->get_next_sample(film_sample)) {
		// Update timer
		double secs = timer.get_secs();
		if ((secs - tstamp) > .5) {
			tstamp = secs;
			unsigned hours = (unsigned)(secs)/3600;
			secs -= hours*3600;
			unsigned minutes = (unsigned)(secs)/60;
			secs -= minutes*60;

			double perc = double(curr_sample)/double(nb_samples)*100;
			printf("Rendering ... %f%%:   [%d:%d:%f]             \r", perc, hours, minutes, secs);
			fflush(stdout);
		}
		curr_sample++;
				
		// Get the camera sample ray...		
		VectorN<D-1> ic;
		Vector2 meh = film->window_coords2image_coords(film_sample.position);
		
		if (D == 2) {
			ic[0] = meh[0];
		} else {
			ic[0] = meh[0];
			ic[1] = meh[1];
		}
		
		Real pdf_camera = 0.;
		Ray<D> r = camera.get_ray(ic, pdf_camera);

		if (pdf_camera == 0.)
			continue;

		film_sample.weight /= pdf_camera;
		
		// ... trace the ray and compute incoming radiance...
		RadianceSampleRecordR radiance_rec = world.Li(r, integrator);

		// ...discard NaNs and infinities
		if (!radiance_rec.sample.radiance.is_finite())
			continue;

		align_to_frame(camera_frame, radiance_rec.sample.radiance);

		// ... and store it in the film
 		film->add_sample(film_sample, radiance_rec);
	}

	//----------------------------------------------------------------------
	// End rendering, stop timer
	{
		double secs = timer.get_secs();
		unsigned hours = (unsigned)(secs)/3600;
		secs -= hours*3600;
		unsigned minutes = (unsigned)(secs)/60;
		secs -= minutes*60;

		printf("Rendering ... [DONE]:   [%d:%d:%f]\n", hours, minutes, secs);
		fflush(stdout);

		if (f_log) {
			fprintf(f_log, "Total time: %d:%d:%f\n", hours, minutes, secs);
			fflush(f_log);
		}
	}
	
	// Save the image
	film->write(name);
}

#endif //_RENDERENGINE_H_
