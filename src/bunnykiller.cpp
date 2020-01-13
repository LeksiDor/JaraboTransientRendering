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

#include "bunnykiller.h"

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "External/TinyXML/tinyxml2.h"

#include "Utils/ProgramParams.h"
#include "Utils/Timer.h"

#include "Camera/Camera.h"
#include "Camera/PinholeCamera.h"
#include "Camera/OrthographicCamera.h"

#include "RayTracing/World.h"

#include "Renderer/RenderEngine.h"
#include "Renderer/TransientRenderer.h"

#include "Film/MultibounceStreakCameraFilm.h"
#include "Film/KernelBasedStreakCameraFilm.h"

#include "Geometry/Primitives/3D/Plane.h"
#include "Geometry/Primitives/3D/Sphere.h"
#include "Geometry/Aggregate/Mesh.h"
#include "Geometry/Aggregate/BVH.h"
#include "Geometry/ObjectPointer.h"

#include "Material/UberMaterial.h"
#include "Material/Reflectance/Phong.h"
#include "Material/Reflectance/Lambertian.h"
#include "Material/Reflectance/Ward.h"
#include "Material/FresnelDielectric.h"
#include "Material/FresnelConductor.h"
#include "Material/FluorescentTrue.h"
#include "Material/Transparent.h"
#include "Material/LambertianTextured.h"
#include "Material/TwoSided.h"

#include "LightSource/PointLightSource.h"
#include "LightSource/DirectionalLightSource.h"
#include "LightSource/CosineLightSource.h"
#include "LightSource/HemisphericalLightSource.h"
#include "LightSource/SpotLightSource.h"
#include "LightSource/RectangularAreaLightSource.h"

#include "Filter/GaussianFilter.h"
#include "Filter/UberFilter.h"

#include "Sampling/StratifiedSampler.h"

#include "Integrator/BidirectionalPathTracing.h"
#include "Integrator/PhotonMapping.h"

#include "Media/HomogeneousMedium.h"
#include "Media/PhaseFunction/IsotropicPhaseFunction.h"
#include "Media/PhaseFunction/HenyeyGreensteinPhaseFunction.h"
#include "Media/PhaseFunction/RayleighPhaseFunction.h"
#include "Media/PhaseFunction/TabulatedPhaseFunction.h"

#include "Color/PolarizationFrame.h"
#include "Color/PolarizedAttenuation.h"
#include "Color/PolarizedSpectrum.h"
#include "Color/FluorescentAttenuation.h"
#include "Color/FluorescentPolarizedAttenuation.h"

#include "Utils/RandomNumbers.h"
#include "Utils/Filesystem.h"

#ifdef _POLARIZATION_
#ifdef _FLUORESCENCE_
typedef PolarizedLight<3> Radiance;
typedef FluorescentPolarizedAttenuation<3> RadianceAttenuation;
typedef PolarizedLight<3> PolarizedRadiance;
#else
typedef PolarizedLight<3> Radiance;
typedef PolarizedAttenuation<3> RadianceAttenuation;
typedef PolarizedLight<3> PolarizedRadiance;
#endif
#else
#ifdef _FLUORESCENCE_
typedef Spectrum Radiance;
typedef FluorescentAttenuation<3> RadianceAttenuation;
#else
typedef Spectrum Radiance;
typedef Spectrum RadianceAttenuation;
#endif
#endif

typedef Film<3, Radiance> FilmS;
typedef StreakCameraFilm<3, Radiance> StreakCameraFilmS;
typedef KernelBasedStreakCameraFilm<3, Radiance> KernelBasedStreakCameraFilmS;
typedef MultibounceStreakFilm<3, Radiance> MultibounceStreakFilmS;

//typedef ParticleTracing<3, Radiance, RadianceAttenuation> ParticleTracingS;

bool parse_args_from_file(const char* filename, std::vector<char*>& args, char*& buffer)
{
	std::ifstream inf(filename);

	if (inf.good()) {
		// Calculate file length
		inf.seekg(0, std::ios::end);
		unsigned int length = inf.tellg();
		inf.seekg(0, std::ios::beg);

		buffer = new char[length + 1];
		inf.read(buffer, length);
		length = inf.gcount(); // Sometimes tellg() returns more bytes than the real size of the file.
							   // Just in case, we correct it using the real number of read bytes
		buffer[length] = '\0';

		// Parse data and create args
		bool inArg = false, inQuot = false;
		for (unsigned int i = 0; i < length; i++) {
			bool sep = buffer[i] == ' ' || buffer[i] == '\n' || buffer[i] == '\r';
			if (inQuot) {
				if (buffer[i] == '"') {
					inQuot = false;
				}
			} else {
				if (buffer[i] == '"') {
					inQuot = true;
				}
				if (inArg) {
					if (sep) {
						buffer[i] = '\0';
						inArg = false;
					}
				} else {
					if (!sep) {
						args.push_back(&buffer[i]);
						inArg = true;
					}
				}
			}
		}
	} else {
		printf("Error: The args file specified is not valid\n");

		return false;
	}

	inf.close();

	return true;
}

Medium<DIM>* load_medium(int &i, int argc, char* argv[])
{
	if (i >= argc) {
		return nullptr;
	}

	if (!strcmp("-homogeneous-medium", argv[i])) {
		i++;

		Spectrum medium_sigma_a(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;
		Spectrum medium_sigma_s(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		Real hg_g = 0, pf_time = 0;
		int type = 0;

		if (!strcmp("-pf-isotropic", argv[i])) {
			i++;
			type = 0;
			hg_g = 0;
		} else if (!strcmp("-pf-hg", argv[i])) {
			i++;
			type = 1;
			hg_g = atof(argv[i]);
			i++;
			if (hg_g == 0.)
				type = 0;
		} else if (!strcmp("-pf-rayleigh", argv[i])) {
			i++;
			type = 2;
		}

		if (i < (argc - 2) && !strcmp("-pf-time", argv[i])) {
			i++;
			pf_time = atof(argv[i]);
			i++;
		}

		switch (type) {
			case 0:
				return new HomogeneousMedium<DIM, PhaseFunction::Isotropic<DIM>>(medium_sigma_a,
						medium_sigma_s, PhaseFunction::Isotropic<DIM>(pf_time));
			case 1:
				return new HomogeneousMedium<DIM, PhaseFunction::HenyeyGreenstein<DIM>>(
						medium_sigma_a, medium_sigma_s,
						PhaseFunction::HenyeyGreenstein<DIM>(hg_g, pf_time));
			case 2:
				return new HomogeneousMedium<DIM, PhaseFunction::Rayleigh<DIM>>(medium_sigma_a,
						medium_sigma_s, PhaseFunction::Rayleigh<DIM>());
		}
	} else if (!strcmp("-mie-medium", argv[i])) {
		i++;
		Spectrum medium_sigma_t(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3; // This way of handling cmd arguments make child Jesus cry, I'm certain of it

		Real pf_time = 0.;
		Spectrum hg_g = Spectrum(0.);
		int type = 0;

		// Read Mie parameters from file
		Spectrum cross_t, cross_s, cross_a;

		char* filename = argv[i];
		try {
			std::ifstream source;
			source.open(argv[i], std::ios_base::in);
			if (!source) {
				throw std::runtime_error("Can't open " + std::string(filename));
			}
			i++;

			// Read bulk parameters
			std::string line;
			std::getline(source, line);
			{
				std::istringstream in(line);
				float r, g, b;

				// Read extinction
				in >> r >> g >> b;
				cross_t = Spectrum(r, g, b);

				// Read scattering
				in >> r >> g >> b;
				cross_s = Spectrum(r, g, b);

				// Read absortion
				in >> r >> g >> b;
				cross_a = Spectrum(r, g, b);

				// Read asymmetry
				in >> r >> g >> b;
				hg_g = Spectrum(r, g, b);
			}

			if (!strcmp("-pf-isotropic", argv[i])) {
				i++;
				type = 0;
			} else if (!strcmp("-pf-hg", argv[i])) {
				i++;
				type = 1;
				if (hg_g == 0.)
					type = 0;
			} else if (!strcmp("-pf-rayleigh", argv[i])) {
				i++;
				type = 2;
			} else if (!strcmp("-pf-tabulated", argv[i])) {
				i++;
				type = 3;
			}

			if (i < (argc - 2) && !strcmp("-pf-time", argv[i])) {
				i++;
				pf_time = atof(argv[i]);
				i++;
			}

			// Scattering and absorption coeff from extinction
			Spectrum medium_sigma_a = medium_sigma_t * (cross_a / cross_t.get_max());
			Spectrum medium_sigma_s = medium_sigma_t * (cross_s / cross_t.get_max());

			switch (type) {
				case 0:
					return new HomogeneousMedium<DIM, PhaseFunction::Isotropic<DIM>>(medium_sigma_a,
							medium_sigma_s, PhaseFunction::Isotropic<DIM>(pf_time));
				case 1:
					return new HomogeneousMedium<DIM, PhaseFunction::HenyeyGreenstein<DIM>>(
							medium_sigma_a, medium_sigma_s,
							PhaseFunction::HenyeyGreenstein<DIM>(hg_g, pf_time));
				case 2:
					return new HomogeneousMedium<DIM, PhaseFunction::Rayleigh<DIM>>(medium_sigma_a,
							medium_sigma_s, PhaseFunction::Rayleigh<DIM>());
				case 3:
					return new HomogeneousMedium<DIM, PhaseFunction::Tabulated<DIM>>(medium_sigma_a,
							medium_sigma_s, PhaseFunction::Tabulated<DIM>(source, pf_time));
			}
		} catch (...) {
			throw std::runtime_error("Error reading Mie parameters from " + std::string(filename));
		}
	}

	return nullptr;
}

// Material loader:
Material3D* load_material(int& i, int argc, char* argv[], Material3D* default_material)
{
	Material3D* mesh_mat = default_material;

	if (!strcmp("-lambertian", argv[i])) {
		i++;

		Spectrum absorp(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		UberMaterial3D *um = new UberMaterial3D();
		um->add_bsdf(new Lambertian3D(absorp));
		mesh_mat = um;
	} else if (!strcmp("-phong", argv[i])) {
		i++;

		Spectrum absorp(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		Real n = atof(argv[i++]);

		UberMaterial3D *um = new UberMaterial3D();
		um->add_bsdf(new Phong(absorp, n));
		mesh_mat = um;
	} else if (!strcmp("-ward", argv[i])) {
		i++;

		Spectrum absorp(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		Real ax = atof(argv[i++]);
		Real ay = atof(argv[i++]);

		UberMaterial3D *um = new UberMaterial3D();
		um->add_bsdf(new Ward(absorp, ax, ay));
		mesh_mat = um;
	} else if (!strcmp("-transparent", argv[i])) {
		i++;

		Spectrum absorp(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		Medium<DIM>* medium = load_medium(i, argc, argv);

		Transparent3D* um;
		if (medium != nullptr) {
			um = new Transparent3D(absorp, medium, nullptr);
		} else {
			um = new Transparent3D(absorp, nullptr, nullptr);
		}
		mesh_mat = um;
	} else if (!strcmp("-dielectric", argv[i])) {
		i++;

		Spectrum absorp(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
		i += 3;

		Real n1 = atof(argv[i++]);
		Real n0 = atof(argv[i++]);

		FresnelDielectric3D *um;
		Medium<DIM>* medium = load_medium(i, argc, argv);

		if (medium != nullptr) {
			um = new FresnelDielectric3D(n1, n0, medium);
		} else {
			um = new FresnelDielectric3D(absorp, n1, n0);
		}
		mesh_mat = um;
	} else if (!strcmp("-conductor", argv[i])) {
		i++;

		Real n = atof(argv[i++]);
		Real k = atof(argv[i++]);

		FresnelConductor3D* um = new FresnelConductor3D(n, k);
		mesh_mat = um;
	} else if (!strcmp("-lambertian-textured", argv[i])) {
		i++;

		LambertianTextured3D *lt = new LambertianTextured3D(Spectrum(atof(argv[i])), argv[i + 1]);
		i += 2;

		mesh_mat = lt;
	} else if (!strcmp("-two-sided", argv[i])) {
		i++;

		mesh_mat = new TwoSidedMaterial3D(load_material(i, argc, argv, default_material));
	}

	return mesh_mat;
}

Spectrum parse_mitsuba_RGB(tinyxml2::XMLElement* e)
{
	float r, g, b;
	const char* intensityValue = e->Attribute("value");

	std::stringstream ss(intensityValue);
	std::string xs, ys, zs;
	std::getline(ss, xs, ',');
	std::getline(ss, ys, ',');
	std::getline(ss, zs);
	r = atof(xs.c_str());
	g = atof(ys.c_str());
	b = atof(zs.c_str());

	return Spectrum(r, g, b);
}

Spectrum parse_mitsuba_spectrum(tinyxml2::XMLElement* e)
{
	float r, g, b = 1;
	//intensity it is
	const char* intensityValue = e->Attribute("value");
	if (intensityValue[0] == '#') {
		unsigned int fullIntensity;
		std::stringstream ss;
		ss << std::hex << &(intensityValue[1]);
		ss >> fullIntensity;
		r = ((fullIntensity >> 16) & 0xFF) / 255.0;
		g = ((fullIntensity >> 8) & 0xFF) / 255.0;
		b = ((fullIntensity) & 0xFF) / 255.0;
	} else {
		r = g = b = atof(intensityValue);
	}

	return Spectrum(r, g, b);
}

Medium<DIM>* load_medium_from_mitsuba(tinyxml2::XMLElement* mediumTag)
{
	if (!strcmp("homogeneous", mediumTag->Attribute("type"))) {
		Spectrum medium_sigma_a(0.);
		Spectrum medium_sigma_s(0.);
		for (tinyxml2::XMLElement* e = mediumTag->FirstChildElement("spectrum"); e != NULL;
				e = e->NextSiblingElement("spectrum")) {
			if (!strcmp("sigmaA", e->Attribute("name"))) {
				medium_sigma_a = parse_mitsuba_spectrum(e);
			} else if (!strcmp("sigmaS", e->Attribute("name"))) {
				medium_sigma_s = parse_mitsuba_spectrum(e);
			}
		}
		for (tinyxml2::XMLElement* e = mediumTag->FirstChildElement("rgb"); e != NULL;
				e = e->NextSiblingElement("rgb")) {
			if (!strcmp("sigmaA", e->Attribute("name"))) {
				medium_sigma_a = parse_mitsuba_RGB(e);
			} else if (!strcmp("sigmaS", e->Attribute("name"))) {
				medium_sigma_s = parse_mitsuba_RGB(e);
			}
		}

		Real hg_g = 0, pf_time = 0;
		int type = 0;

		tinyxml2::XMLElement* phaseTag = mediumTag->FirstChildElement("phase");
		if (phaseTag != NULL) {
			if (!strcmp(phaseTag->Attribute("type"), "isotropic")) {
				type = 0;
				hg_g = 0;
			} else if (!strcmp(phaseTag->Attribute("type"), "hg")) {
				type = 1;
				hg_g = 0;
				tinyxml2::XMLElement* hgfloat = phaseTag->FirstChildElement("float");
				if (hgfloat != NULL && !strcmp(hgfloat->Attribute("name"), "g"))
					hg_g = atof(hgfloat->Attribute("value"));
				if (hg_g == 0)
					type = 0;
			} else if (!strcmp(phaseTag->Attribute("type"), "rayleigh")) {
				type = 2;
			}
		}

		switch (type) {
			case 0:
				fprintf(stdout, "Isotropic medium with sigma a ");
				medium_sigma_a.print(stdout);
				fprintf(stdout, " , sigma t ");
				medium_sigma_s.print(stdout);
				return new HomogeneousMedium<DIM, PhaseFunction::Isotropic<DIM>>(medium_sigma_a,
						medium_sigma_s, PhaseFunction::Isotropic<DIM>(pf_time));
			case 1:
				fprintf(stdout, "Henyey Greenstein medium with hg_g %f sigma a ", hg_g);
				medium_sigma_a.print(stdout);
				fprintf(stdout, " , sigma t ");
				medium_sigma_s.print(stdout);
				return new HomogeneousMedium<DIM, PhaseFunction::HenyeyGreenstein<DIM>>(
						medium_sigma_a, medium_sigma_s,
						PhaseFunction::HenyeyGreenstein<DIM>(hg_g, pf_time));
			case 2:
				fprintf(stdout, "Rayleigh medium with sigma a ");
				medium_sigma_a.print(stdout);
				fprintf(stdout, " , sigma t ");
				medium_sigma_s.print(stdout);
				return new HomogeneousMedium<DIM, PhaseFunction::Rayleigh<DIM>>(medium_sigma_a,
						medium_sigma_s, PhaseFunction::Rayleigh<DIM>());
		}
	}

	return nullptr;
}

Material3D* load_material_from_mitsuba(tinyxml2::XMLElement *bsdf, UberMaterial3D *m_grey)
{
	if (!strcmp("diffuse", bsdf->Attribute("type"))) {
		//Without textures
		Spectrum diffuseColor(0.6);
		if (bsdf->FirstChildElement("srgb") != nullptr) {
			std::string val = bsdf->FirstChildElement("srgb")->Attribute("value");
			int r = strtol(val.substr(1, 2).c_str(), nullptr, 16);
			int g = strtol(val.substr(3, 2).c_str(), nullptr, 16);
			int b = strtol(val.substr(5, 2).c_str(), nullptr, 16);
			diffuseColor = Spectrum(r / 255.0, g / 255.0, b / 255.0);
		}

		if (bsdf->FirstChildElement("rgb") != nullptr) {
			diffuseColor = parse_mitsuba_RGB(bsdf->FirstChildElement("rgb"));
		}

		Vector3 diffRGB = diffuseColor.to_rgb();

		printf("Diffuse material loaded with RGB: %f %f %f\n", diffRGB[0], diffRGB[2], diffRGB[3]);

		UberMaterial3D *um = new UberMaterial3D();
		um->add_bsdf(new Lambertian3D(diffuseColor));

		return um;
	} else if (!strcmp("dielectric", bsdf->Attribute("type"))) {
		//|TODO esto en mitsuba no es asi
		float n0 = 1., n1 = 1.;
		Spectrum absorption;
		for (tinyxml2::XMLElement* e = bsdf->FirstChildElement("float"); e != nullptr;
				e = e->NextSiblingElement("float")) {
			if (!strcmp(e->Attribute("name"), "n0")) {
				n0 = atof(e->Attribute("value"));
			} else if (!strcmp(e->Attribute("name"), "n1")) {
				n1 = atof(e->Attribute("value"));
			} else if (!strcmp(e->Attribute("name"), "absorption")) {
				absorption = Spectrum(atof(e->Attribute("value")));
			}
		}

		tinyxml2::XMLElement* mediumTag = bsdf->FirstChildElement("medium");

		Medium<DIM> *medium = nullptr;
		if (mediumTag != nullptr) {
			medium = load_medium_from_mitsuba(mediumTag);
		}

		if (medium != nullptr) {
			printf("Dielectric material loaded with n0: %f, n1: %f, and the last medium.\n", n0,
					n1);

			return new FresnelDielectric3D(n1, n0, medium);
		} else {
			printf("Dielectric material loaded with n0: %f, n1: %f, absortion:", n0, n1);
			absorption.print(stdout);
			return new FresnelDielectric3D(absorption, n1, n0);
		}
#ifdef _FLUORESCENCE_
	} else if (!strcmp("fluorescent", bsdf->Attribute("type"))) {
		FluorescentMatrix direct(1.);
		FluorescentMatrix reemitted(0.);
		Real reemisionTime = 0;
		if (bsdf->FirstChildElement("float") != nullptr && !strcmp(bsdf->FirstChildElement("float")->Attribute("name"), "reemitTime")) {
			reemisionTime = atof(bsdf->FirstChildElement("float")->Attribute("value"));
		}
		for (tinyxml2::XMLElement* e = bsdf->FirstChildElement("matrix"); e != nullptr; e = e->NextSiblingElement("matrix")) {
			Real contents[9]; // HARD-CODED FOR SPECTRUM 3!!!! |TODO
			std::stringstream ss(e->Attribute("value"));
			std::string content;
			for (int i = 0; i < 9; i++) {
				ss >> content;
				contents[i] = atof(content.c_str());
			}

			if (!strcmp(e->Attribute("name"), "direct")) {
				direct = FluorescentMatrix(contents);
			} else if (!strcmp(e->Attribute("name"), "reemitted")) {
				reemitted = FluorescentMatrix(contents);
			}
		}

		printf("New fluorescent material loaded with time %f  and matrices \n", reemisionTime);
		printf("Direct: \n");
		direct.print(stdout);
		printf("Reemitted: \n");
		reemitted.print(stdout);

		Fluorescent3D *fluorMaterial= new Fluorescent3D(direct,reemitted, reemisionTime);

		return fluorMaterial;
#endif
	} else if (!strcmp("conductor", bsdf->Attribute("type"))) {
		Real n = 1., k = 1.;

		for (tinyxml2::XMLElement* e = bsdf->FirstChildElement("spectrum"); e != NULL;
				e = e->NextSiblingElement("spectrum")) {
			Radiance s = parse_mitsuba_spectrum(e);

			// UNSUPPORTED SPECTRUM ETA AND K IN THIS COMPILATION, AVERAGING
			if (!strcmp("eta", e->Attribute("name"))) {
				n = s.avg();
			} else if (!strcmp("k", e->Attribute("name"))) {
				k = s.avg();
			}
		}

		printf("New conductor material loaded with n: %f and k: %f", n, k);

		return new FresnelConductor3D(n, k);
	}

	return m_grey;
}

void parse_args_from_mitsuba(tinyxml2::XMLDocument& doc, ProgramParams& p, World<DIM, Radiance>& w,
	UberMaterial3D& m_grey, std::string& filename)
{
	// ASSUMING EMBREE!!!!
	tinyxml2::XMLElement* base = doc.FirstChildElement("scene");
	tinyxml2::XMLElement* element;

	/********************** BSDFS *********************/
	std::map<std::string, Material3D*> idMaterials;
	for (tinyxml2::XMLElement* e = base->FirstChildElement("bsdf"); e != nullptr;
			e = e->NextSiblingElement("bsdf")) {
		Material3D *mat = load_material_from_mitsuba(e, &m_grey);
		std::string id = e->Attribute("id");
		printf("Material assigned to ID %s\n", id.c_str());

		idMaterials.insert(std::pair<std::string, Material3D*>(id, mat));
	}

	/********************** SHAPES *********************/
	for (tinyxml2::XMLElement* e = base->FirstChildElement("shape"); e != nullptr;
			e = e->NextSiblingElement("shape")) {
		// Load OBJ mesh
		if (!strcmp("obj", e->Attribute("type"))) {
			Material3D *mat = &m_grey;
			// A bsdf can be also declared into the shape
			tinyxml2::XMLElement* bsdf = e->FirstChildElement("bsdf");
			bool bsdfDeclared = false;
			if (bsdf != nullptr) {
				bsdfDeclared = true; //if a ref is declared too, a warning will be written on screen
				mat = load_material_from_mitsuba(bsdf, &m_grey);
			}

			tinyxml2::XMLElement* ref = e->FirstChildElement("ref");
			if (ref != nullptr) {
				if (ref->Attribute("id") != nullptr) {
					if (bsdfDeclared)
						printf(
								"WARNING: Bsdf AND ref both declared on a shape. Using ref by default.\n");

					std::string id = ref->Attribute("id");
					std::map<std::string, Material3D*>::iterator it = idMaterials.find(id);
					if (it != idMaterials.end()) {
						printf("Using material ID %s\n", id.c_str());

						mat = it->second;
					} else {
						printf("WARNING: Shape using unexisting reference\n");
					}
				}
			}

			for (tinyxml2::XMLElement* eobj = e->FirstChildElement("string"); eobj != nullptr;
					eobj = eobj->NextSiblingElement("string")) {
				if (!strcmp("filename", eobj->Attribute("name"))) {
					std::string fname = eobj->Attribute("value");
					printf("Loading mesh %s\n", fname.c_str());

					//Transforms not enabled on this raytracer (extend?)
					//meshes as light emitters not implemented on this raytracer
#ifdef MITSUBA_ABSOLUTE_ROUTE
					w.add_triangle_mesh(fname, mat);
#else
					std::stringstream ss;
					ss << filename << fname;
					std::string realRoute = ss.str();
#ifdef _USE_EMBREE_
					w.add_triangle_mesh(realRoute, mat);
#else
					Mesh *m = new Mesh(realRoute, mat);
					w.add_object(Object3DPointer(m));
#endif
#endif
				}
			}
		} else if (!strcmp("sphere", e->Attribute("type"))) {
			// Sphere
			throw std::runtime_error(
					"Unsupported in the Embree-based ray tracing. Transform it in a triangle mesh.");
		}
	}

	/******************** INTEGRATORS ******************/

	// The integrator is right now defined as a constant, so this gets commented
	element = base->FirstChildElement("integrator");

	if (element != nullptr) {
		std::string integratorType = element->Attribute("type");
#ifdef _USE_PM_
		if (!strcmp(integratorType.c_str(), "photonmapper") ||
				!strcmp(integratorType.c_str(), "ppm")) {
			p.vol_path_tracing_samples = 1;

			for (tinyxml2::XMLElement* e = element->FirstChildElement("integer");
					e != nullptr; e = e->NextSiblingElement("integer")) {
				if (!strcmp("lookupSize", e->Attribute("name"))) {
					p.pm_photons_nn = atoi(e->Attribute("value"));

					printf("PM lookup size: %d\n", p.pm_photons_nn);
				} else if (!strcmp("globalPhotons", e->Attribute("name"))) {
					p.pm_photon_shots = atoi(e->Attribute("value"));

					printf("PM photon shots: %d\n", p.pm_photon_shots);
				}
				/*else if (!strcmp("maxPasses", e->Attribute("name"))) {
				 p.vol_path_tracing_samples = atoi(e->Attribute("value"));
				 }*/
			}

			for (tinyxml2::XMLElement* e = element->FirstChildElement("float");
					e != nullptr; e = e->NextSiblingElement("float")) {
				if (!strcmp("alpha", e->Attribute("name"))) {
					p.pm_kernel_shrinkage = atof(e->Attribute("value"));

					printf("PM kernel shrinkage: %f\n", p.pm_kernel_shrinkage);
				}
			}

			printf("PM path tracing samples: %d\n", p.vol_path_tracing_samples);
		}
#else
		if (!strcmp(integratorType.c_str(), "bdpt")) {
			for (tinyxml2::XMLElement* e = element->FirstChildElement("integer"); e != nullptr; e =
					e->NextSiblingElement("integer")) {
				if (!strcmp("samples", e->Attribute("name"))) {
					p.vol_path_tracing_samples = atoi(e->Attribute("value"));

					printf("BDPT samples: %d\n", p.vol_path_tracing_samples);
				}
			}
		}
#endif
		else {
			throw std::runtime_error("Unknown integrator.");
		}
	}

	/********************** SENSOR **********************/

	//Supported: Perspective camera
	element = base->FirstChildElement("sensor");

	if (element != nullptr) {
		tinyxml2::XMLElement* sensorElement;
		std::string sensorType = element->Attribute("type");
		for (tinyxml2::XMLElement* e = element->FirstChildElement("float"); e != nullptr;
				e = e->NextSiblingElement("float")) {
			if (!strcmp("fov", e->Attribute("name"))) {
				Real fov = M_PI / 180. * std::stof(e->Attribute("value")) / 2;
				p.camera_view_plane_dist = 1. / tan(fov);

				printf("Fov: %f, View plane dist: %f\n", fov, p.camera_view_plane_dist);
			}
		}

		sensorElement = element->FirstChildElement("sampler");
		if (sensorElement != nullptr) {
			//std::string samplerType = sensorElement->Attribute("type"); //What the fuck is that

			for (tinyxml2::XMLElement* e = sensorElement->FirstChildElement("integer");
					e != nullptr; e = e->NextSiblingElement("integer")) {
				if (!strcmp("sampleCount", e->Attribute("name"))) {
					p.sqrt_spp = atoi(e->Attribute("value"));

					printf("Camera spp: %d\n", p.sqrt_spp);
				}
			}
		}

		sensorElement = element->FirstChildElement("film");
		if (sensorElement != NULL) {
			//i dont care about the type, i guess
			for (tinyxml2::XMLElement* e = sensorElement->FirstChildElement("integer");
					e != nullptr; e = e->NextSiblingElement("integer")) {
				if (!strcmp("width", e->Attribute("name"))) {
					p.width = atoi(e->Attribute("value"));

					printf("Camera width: %d\n", p.width);
				} else if (!strcmp("height", e->Attribute("name"))) {
					p.height = atoi(e->Attribute("value"));

					printf("Camera height: %d\n", p.height);
				} else if (!strcmp("time", e->Attribute("name"))) {
					p.time = atoi(e->Attribute("value"));

					printf("Camera time res: %d\n", p.time);
				}
			}

			for (tinyxml2::XMLElement* e = sensorElement->FirstChildElement("float"); e != nullptr;
					e = e->NextSiblingElement("float")) {
				if (!strcmp("timeResolution", e->Attribute("name"))) {
					p.time_resolution = atof(e->Attribute("value"));

					printf("Camera time resolution: %f\n", p.time_resolution);
				}
				if (!strcmp("timeOffset", e->Attribute("name"))) {
					p.film_offset = atof(e->Attribute("value"));

					printf("Film offset: %f\n", p.film_offset);
				}
			}

			for (tinyxml2::XMLElement* e = sensorElement->FirstChildElement("boolean");
					e != nullptr; e = e->NextSiblingElement("boolean")) {
				if (!strcmp("steadyState", e->Attribute("name"))) {
					p.transient = strcmp(e->Attribute("value"), "true");
				}
			}
		}

		sensorElement = element->FirstChildElement("transform");
		if (sensorElement != nullptr) {
			//I only support lookat for now (i mean, if you dont define a camera using lookat you are definitely a bad person)
			tinyxml2::XMLElement* transformType = sensorElement->FirstChildElement("lookat");
			if (transformType != nullptr) {
				if (transformType->Attribute("target") != nullptr) {
					std::stringstream ss(transformType->Attribute("target"));
					std::string x, y, z;
					std::getline(ss, x, ',');
					std::getline(ss, y, ',');
					std::getline(ss, z);
					p.camera_looking_at = Vector3(atof(x.c_str()), atof(y.c_str()),
							atof(z.c_str()));

					printf("Looking at: %f %f %f\n", p.camera_looking_at[0], p.camera_looking_at[1],
							p.camera_looking_at[2]);
				}
				if (transformType->Attribute("origin") != nullptr) {
					std::stringstream ss(transformType->Attribute("origin"));
					std::string x, y, z;
					std::getline(ss, x, ',');
					std::getline(ss, y, ',');
					std::getline(ss, z);
					p.camera_position = Vector3(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));

					printf("Camera position: %f %f %f\n", p.camera_position[0],
							p.camera_position[1], p.camera_position[2]);
				}
				if (transformType->Attribute("up") != nullptr) {
					std::stringstream ss(transformType->Attribute("up"));
					std::string x, y, z;
					std::getline(ss, x, ',');
					std::getline(ss, y, ',');
					std::getline(ss, z);
					p.camera_up = Vector3(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));

					printf("Up vector: %f %f %f\n", p.camera_up[0], p.camera_up[1], p.camera_up[2]);
				}
			}
		}

		//NEED TO ADD TRANSFORMS!!!!!!
		//Color spectra? blackbodies?
	}

	/*********************** EMITTERS *********************/

	//Supported: Point light source, spot light source
	for (tinyxml2::XMLElement* e = base->FirstChildElement("emitter"); e != nullptr;
			e = e->NextSiblingElement("emitter")) {
		if (!strcmp(e->Attribute("type"), "point")) {
			tinyxml2::XMLElement* lightElement = e->FirstChildElement("spectrum");
			Spectrum color(1.);
			if (lightElement != nullptr) {
				//intensity it is
				color = parse_mitsuba_spectrum(lightElement);
			}
			lightElement = e->FirstChildElement("rgb");

			if (lightElement != nullptr) {
				//intensity it is
				color = parse_mitsuba_RGB(lightElement);
			}
			lightElement = e->FirstChildElement("point");

			if (lightElement != nullptr) {
				if (!(lightElement->Attribute("x") == nullptr
						|| lightElement->Attribute("y") == nullptr
						|| lightElement->Attribute("z") == nullptr)) {
					Vector3 pos = Vector3(atof(lightElement->Attribute("x")),
							atof(lightElement->Attribute("y")), atof(lightElement->Attribute("z")));

#ifdef _POLARIZATION_
					w.add_light(new PointLightSource<3, Radiance>(&w, pos, Radiance(color, color, PolarizationFrame<3>(normalize(p.camera_up), normalize(p.camera_looking_at - p.camera_position)))));
#else
					w.add_light(new PointLightSource<3, Radiance>(&w, pos, color));

					Vector3 rgb = color.to_rgb();
					printf("Light source (point): %f %f %f (%f, %f, %f)\n", pos[0], pos[1], pos[2],
							rgb[0], rgb[1], rgb[2]);
#endif
				}
			}
		} else if (!strcmp(e->Attribute("type"), "spot")) {
			float intensity = 1;
			float cutoffAngle = 20; //In deg though
			Vector3 from(0, 0, 0);
			Vector3 to(0, 0, 1);

			tinyxml2::XMLElement* lightElement = e->FirstChildElement("spectrum");
			if (lightElement != nullptr) {
				//intensity it is
				intensity = atof(lightElement->Attribute("value"));
			}
			lightElement = e->FirstChildElement("cutoffAngle");
			if (lightElement != nullptr) {
				//intensity it is
				cutoffAngle = atof(lightElement->Attribute("value"));
			}
			lightElement = e->FirstChildElement("transform");
			if (lightElement != nullptr)
				lightElement = e->FirstChildElement("lookat");
			if (lightElement != nullptr) {
				if (lightElement->Attribute("target") != nullptr) {
					std::stringstream ss(lightElement->Attribute("target"));
					std::string x, y, z;
					std::getline(ss, x, ',');
					std::getline(ss, y, ',');
					std::getline(ss, z);
					to = Vector3(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));
				}
				if (lightElement->Attribute("origin") != nullptr) {
					std::stringstream ss(lightElement->Attribute("origin"));
					std::string x, y, z;
					std::getline(ss, x, ',');
					std::getline(ss, y, ',');
					std::getline(ss, z);
					from = Vector3(atof(x.c_str()), atof(y.c_str()), atof(z.c_str()));
				}
			}

			w.add_light(
					new SpotLightSource<3, Radiance>(&w, from, to, cutoffAngle,
							Radiance(intensity)));
			printf("Light source (spot): \n");
		}
	}
}

inline void parse_args(int argc, char* argv[], ProgramParams& p, World<DIM, Radiance>& w,
	UberMaterial3D& m_grey)
{
	for (int i = 1; i < argc;) {
		// Store current loop iter
		const int i_start = i;

		// Default material BSDF
		if (!strcmp("-lambertian-rho", argv[i])) {
			i++;
			p.grey = Spectrum(atof(argv[i++]));
			continue;
		}
#ifdef _USE_EMBREE_
		// Mesh
		if (!strcmp("-name-mesh", argv[i])) {
			i++;
			char* mesh_name = argv[i++];

			Material3D* mat = load_material(i, argc, argv, static_cast<Material3D*>(&m_grey));
			w.add_triangle_mesh(mesh_name, mat);
			continue;
		}
		// Sphere
		if (!strcmp("-sphere", argv[i])) {
			i++;
			throw std::runtime_error(
					"Unsupported in the Embree-based ray tracing. Transform it in a triangle mesh.");
		}
		// Plane
		if (!strcmp("-plane", argv[i])) {
			i++;
			throw std::runtime_error(
					"Unsupported in the Embree-based ray tracing. Transform it in a triangle mesh.");
		}
#else
		// Mesh
		if (!strcmp("-name-mesh", argv[i])) {
			i++;
			char* mesh_name = argv[i++];

			Material3D* mat = load_material(i, argc, argv, static_cast<Material3D*>(&m_grey));
			Mesh* m = new Mesh(mesh_name, mat);
			w.add_object(Object3DPointer(m));
			continue;
		}
		// Sphere
		if (!strcmp("-sphere", argv[i])) {
			i++;
			Vector3 p(atof(argv[i]), atof(argv[i+1]), atof(argv[i+2]));
			i += 3;
			Real r = atof(argv[i++]);

			Material3D* mat = load_material(i, argc, argv, static_cast<Material3D*>(&m_grey));
			w.add_object(Object3DPointer(new Sphere(p, r, mat)));
			continue;
		}
		// Plane
		if (!strcmp("-plane", argv[i])) {
			i++;
			Vector3 p(atof(argv[i]), atof(argv[i+1]), atof(argv[i+2]));
			i += 3;
			Vector3 n(atof(argv[i]), atof(argv[i+1]), atof(argv[i+2]));
			i += 3;

			Material3D *mat = load_material(i, argc, argv, static_cast<Material3D*>(&m_grey));
			w.add_object(Object3DPointer(new Plane(p, n, mat)));
			continue;
		}
#endif
		/// Film parameters
		if (!strcmp("-film-name", argv[i])) {
			i++;
			p.image_file = argv[i++];
			continue;
		}
		if (!strcmp("-film-size-x", argv[i])) {
			i++;
			p.width = atoi(argv[i++]);
			p.height = (DIM == 2) ? 1 : p.width;
			continue;
		}
		if (!strcmp("-film-size-y", argv[i])) {
			i++;
			p.height = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-film-size-t", argv[i])) {
			i++;
			p.time = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-film-aspect-ratio", argv[i])) {
			i++;
			p.aspectRatioY = atof(argv[i++]);
			continue;
		}
		if (!strcmp("-film-exposure", argv[i])) {
			i++;
			p.time_resolution = atof(argv[i++]);
			continue;
		}
		if (!strcmp("-film-scanline", argv[i])) {
			i++;
			p.scanline = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-film-offset", argv[i])) {
			i++;
			p.film_offset = atof(argv[i++]);
			continue;
		}
		if (!strcmp("-film-single-pixel", argv[i])) {
			i++;
			p.single_pixel = true;
			p.single_pixel_x = atoi(argv[i++]);
			p.single_pixel_y = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-film-reveal", argv[i])) {
			i++;
			p.film_usereveal = true;
			continue;
		}
		if (!strcmp("-film-store-depth", argv[i])) {
			i++;
			p.film_depth = true;
			continue;
		}
		if (!strcmp("-film-store-positions", argv[i])) {
			i++;
			p.film_positions = true;
			continue;
		}
		if (!strcmp("-film-store-normals", argv[i])) {
			i++;
			p.film_normals = true;
			continue;
		}

		/// Camera parameters
		if (!strcmp("-camera-position", argv[i])) {
			i++;
			if (argc > i + 2) {
				p.camera_position = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
				i += 3;
			}
			continue;
		}
		if (!strcmp("-camera-focus", argv[i])) {
			i++;
			if (argc > i + 2) {
				p.camera_looking_at = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
				i += 3;
			}
			continue;
		}
		if (!strcmp("-camera-up", argv[i])) {
			i++;
			if (argc > i + 2) {
				p.camera_up = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
				i += 3;
				p.camera_up.normalize();
			}
			continue;
		}
		if (!strcmp("-camera-view-plane", argv[i])) {
			i++;
			p.camera_view_plane_dist = atof(argv[i++]);
			continue;
		}
		if (!strcmp("-camera-spp", argv[i])) {
			i++;
			p.sqrt_spp = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-camera-fov", argv[i])) {
			i++;
			Real fov = M_PI / 180. * atof(argv[i++]) / 2;
			p.camera_view_plane_dist = 1. / tan(fov);
			continue;
		}
		if (!strcmp("-camera-view-plane-size", argv[i])) {
			i++;
			p.camera_w = atof(argv[i++]);
			p.camera_h = atof(argv[i++]);
			continue;
		}

		/// Lights
		if (!strcmp("-cosine-light-source", argv[i])) {
			i++;
			Vector3 pos = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Vector3 dir = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			dir.normalize();
			Radiance intensity(atof(argv[i++]));

			w.add_light(new CosineLightSource<3, Radiance>(&w, pos, dir, intensity));
			continue;
		}
		if (!strcmp("-hemispherical-light-source", argv[i])) {
			i++;
			Vector3 pos = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Vector3 dir = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			dir.normalize();
			Radiance intensity(atof(argv[i++]));

			w.add_light(new HemisphericalLightSource<3, Radiance>(&w, pos, dir, intensity));
			continue;
		}
		if (!strcmp("-point-light-source", argv[i])) {
			i++;
			Vector3 pos = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
#ifdef _POLARIZATION_
			PolarizationFrame<DIM> frame(normalize(p.camera_up), normalize(p.camera_looking_at - p.camera_position));
			Radiance intensity = Radiance(atof(argv[i++]), frame);
#else
			Radiance intensity(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
#endif
			w.add_light(new PointLightSource<3, Radiance>(&w, pos, intensity));
			continue;
		}
		if (!strcmp("-point-light-source-polarized", argv[i])) {
			i++;
			Vector3 pos = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Spectrum I = Spectrum(atof(argv[i++]));
			Spectrum Q = Spectrum(atof(argv[i++]));
			Spectrum U = Spectrum(atof(argv[i++]));
			Spectrum V = Spectrum(atof(argv[i++]));
#ifdef _POLARIZATION_
			PolarizationFrame<DIM> frame(normalize(p.camera_up), normalize(p.camera_looking_at - p.camera_position));
			Radiance intensity(I, Q, U, V, frame);

			w.add_light(new PointLightSource<3, Radiance>(&w, pos, intensity));
#else
			Radiance intensity(I);

			w.add_light(new PointLightSource<3, Radiance>(&w, pos, intensity));
#endif
			continue;
		}
		if (!strcmp("-spot-light-source", argv[i])) {
			i++;
			Vector3 pos = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Vector3 dir = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			dir.normalize();

			Real angle = atof(argv[i++]);
			Radiance intensity(atof(argv[i++]));

			w.add_light(new SpotLightSource<3, Radiance>(&w, pos, dir, angle, intensity));
			continue;
		}
		if (!strcmp("-rectangular-area-light-source", argv[i])) {
			i++;
			Vector3 pos(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Vector3 dir(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			Vector2 size(atof(argv[i]), atof(argv[i + 1]));
			i += 2;
			Radiance intensity(atof(argv[i++]));

			Texture* tex = 0;
			if (argc > i + 2 && !strcmp("-emission-texture", argv[i++])) {
				tex = new Texture(Imaging::load<Real>(argv[i++]));
			}

			w.add_light(
					new RectangularAreaLightSource<Radiance>(&w, pos, (dir - pos).normalized(),
							size, intensity, tex));
			continue;
		}
		if (!strcmp("-directional-light-source", argv[i])) {
			i++;
			Vector3 dir(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			dir.normalize();
			Radiance intensity(atof(argv[i++]));

#ifdef _POLARIZATION_
			w.add_light(new DirectionalLightSource<3, Radiance>(&w, dir, intensity));
#else
			w.add_light(new DirectionalLightSource<3, Radiance>(&w, dir, intensity));
#endif
			continue;
		}

		/// VPL generation
		if (!strcmp("-virtual-point-light-emitter", argv[i])) {
			i++;
			p.emitter = true;

			p.emitter_position = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			p.emitter_looking_at = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;
			p.emitter_up = Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
			i += 3;

			if (!strcmp("-ortographic", argv[i])) {
				i++;
				p.emitter_w = atof(argv[i++]);
				p.emitter_h = atof(argv[i++]);

				p.emitter_othographic = true;
			} else if (!strcmp("-perspective", argv[i])) {
				i++;
				Real fov = M_PI / 180. * atof(argv[i++]) / 2;
				p.emitter_view_plane_dist = 1. / tan(fov);
			} else {
				throw std::runtime_error("Unknown VPL emitter type");
			}

			continue;
		}
		if (!strcmp("-virtual-point-light-source", argv[i])) {
			i++;
			if (!strcmp("-image-coordinates", argv[i])) {
				i++;
				p.image_coordinates_vpl.emplace_back(
						Vector3(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2])));
				i += 3;
			} else {
				Vector3 pos(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
				i += 3;
				Vector3 dir(atof(argv[i]), atof(argv[i + 1]), atof(argv[i + 2]));
				i += 3;
				Spectrum power(atof(argv[i++]));

				w.add_light(
						new CosineLightSource<DIM, Radiance>(&w, Ray<DIM>(pos, dir.normalized()),
								power));
			}
			continue;
		}

		/// Sampling
		if (!strcmp("-transient-state", argv[i])) {
			i++;
			p.transient = true;
			continue;
		}
		if (!strcmp("-steady-state", argv[i])) {
			i++;
			p.transient = false;
			continue;
		}
		if (!strcmp("-time-sampling", argv[i])) {
			i++;
			p.time_sampling = true;
			continue;
		}
		if (!strcmp("-standard-sampling", argv[i])) {
			i++;
			p.time_sampling = false;
			continue;
		}

		/// Log file
		if (!strcmp("-log-name", argv[i])) {
			i++;
			p.log_file = argv[i++];
			continue;
		}

		/// Random numbers
		if (!strcmp("-new-seed", argv[i])) {
			i++;
			p.new_seed = atoi(argv[i++]);
			continue;
		}

		/// Participating media
		if (!strcmp("-homogeneous-medium", argv[i]) || !strcmp("-mie-medium", argv[i])) {
			w.set_medium(load_medium(i, argc, argv));
			continue;
		}
		if (!strcmp("-max-nb-bounces", argv[i])) {
			i++;
			p.max_nb_bounces = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-scattering-level", argv[i])) {
			i++;
			if (!strcmp("multiple", argv[i])) {
				p.scat_level = Scattering::MULTIPLE;
			} else if (!strcmp("single", argv[i])) {
				p.scat_level = Scattering::SINGLE;
			} else if (!strcmp("none", argv[i])) {
				p.scat_level = Scattering::NONE;
			} else {
				p.scat_level = Scattering::ALL;
			}
			i++;

			continue;
		}

		// Integrator
		if (!strcmp("-bidirectional-path-tracing", argv[i])) {
			i++;
#ifdef _USE_PM_
			printf("WARNING: BDPT is not working in this compilation. Recompile undefining _USE_PM_ in bunnykiller.h");
#else
			p.vol_path_tracing_samples = atoi(argv[i++]);
#endif
			continue;
		}
		if (!strcmp("-photon-mapping", argv[i])) {
			i++;
#ifdef _USE_PM_
			p.pm_photon_shots = atoi(argv[i++]);
			p.pm_photons_nn = atoi(argv[i++]);
			p.pm_kernel_shrinkage = atof(argv[i++]);
			p.vol_path_tracing_samples = atoi(argv[i++]);
#else
			printf(
					"WARNING: PPM is not working in this compilation. Recompile defining _USE_PM_ in bunnykiller.h");
#endif
			continue;
		}

		// Transient params
		if (!strcmp("-multibounce-streak", argv[i])) {
			i++;
			p.is_multibounce_streak = true;
			p.nb_multibounces = atoi(argv[i++]);
			continue;
		}
		if (!strcmp("-kernel-based-streak", argv[i])) {
			i++;
			p.streak_kernel_based = true;

			while (!strcmp("-r", argv[i]) || !strcmp("-nn", argv[i]) || !strcmp("-alpha", argv[i])) {
				if (!strcmp("-r", argv[i])) {
					i++;
					p.streak_kernel_radius = atof(argv[i++]);
				}

				if (!strcmp("-nn", argv[i])) {
					i++;
					p.streak_nn = atoi(argv[i++]);
				}

				if (!strcmp("-alpha", argv[i])) {
					i++;
					p.streak_alpha = atof(argv[i++]);
				}
			}

			continue;
		}
		if (!strcmp("-camera-unwarp", argv[i])) {
			i++;
			p.camera_unwarp = true;
			continue;
		}

		// Args file parsed IN ARG POSITION. It can be seen as an extension of the standard command line
		if (!strcmp("-parse-args-file", argv[i])) {
			i++;
			char* filename = argv[i++];

			std::vector<char*> args;

			char* buffer = nullptr;
			if (parse_args_from_file(filename, args, buffer)) {
				parse_args(args.size(), args.data(), p, w, m_grey);
			}
			delete buffer;

			continue;
		}

		// XML in Mitsuba-like format
		if (!strcmp("-parse-xml-file", argv[i])) {
			i++;
			char* filename = argv[i++];

			tinyxml2::XMLDocument doc;
			printf("Parsing args from XML:\n");
			int errcode = doc.LoadFile(filename);
			if (errcode == 0) {
				//Get relative path
				int lastOccur = 0;
				int ci = 0;
				while (filename[ci] != '\0') {
					if (filename[ci] == '/' || filename[ci] == '\\')
						lastOccur = ci + 1; //Want to include the / or \ in the size too
					ci++;
				}
				std::string relPath(filename, lastOccur);
				parse_args_from_mitsuba(doc, p, w, m_grey, relPath);
			} else {
				fprintf(stderr,
						"WARNING: Unable to load XML params file (Maybe it is bad formed). Error code: %d. Error msg:%s\n",
						errcode, doc.GetErrorStr1());
			}

			printf("\n");
			continue;
		}

		// Keep advancing arguments even if none matched
		if (i == i_start) {
			i++;
		}
	}

	// Add BSDF to default material
	Lambertian3D rho_grey(p.grey);
	m_grey.add_bsdf(&rho_grey);
}

int main(int argc, char* argv[])
{
	Timer timer;

	// Prints all input arguments.
	for (int i = 0; i < argc; i++) {
		fprintf(stdout, "argv[%d]: %s\n", i, argv[i]);
	}
	fprintf(stdout, "\n");

	// WORLD
	World<DIM, Radiance> w = World<DIM, Radiance>();

	// Default material when unspecified
	UberMaterial3D m_grey;

	// Read all program parameters
	ProgramParams p(DIM);
	parse_args(argc, argv, p, w, m_grey);

	bool sensor_mode = false;

	// CREATE ALL FILENAMES INCLUDING FOLDER
	if (!p.log_file) {
		p.log_file = (char*) malloc(sizeof(char) * 512);
		sprintf(p.log_file, "%s.log", p.image_file);
	}
	
	std::string log_path = Filesystem::file_path(p.log_file);
	if (!log_path.empty()) {
		Filesystem::create_directory(log_path.c_str());
	}

	std::string img_path = Filesystem::file_path(p.image_file);
	if (!img_path.empty()) {
		Filesystem::create_directory(img_path.c_str());
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	FILE* f_log = fopen(p.log_file, "w+");

	//LOG FILM DATA
	fprintf(f_log, "FILM\n");
	fprintf(f_log, "====\n");
	fprintf(f_log, "Name: %s\n", p.image_file);
	fprintf(f_log, "Size(x-y-t): %d %d %d\n", p.width, p.height, p.time);
	fprintf(f_log, "Exposure: %.10f\n", p.time_resolution);
	fprintf(f_log, "t0: %.10f\n", p.film_offset);
	fprintf(f_log, "Sqrt-spp: %d\n", p.sqrt_spp);
	fprintf(f_log, "Seed: %d\n", (unsigned int) p.new_seed);
	if (sensor_mode) {
		fprintf(f_log, "Camera-unwarp: --\n");
	} else if (p.camera_unwarp) {
		fprintf(f_log, "Camera-unwarp: TRUE\n");
	} else {
		fprintf(f_log, "Camera-unwarp: FALSE\n");
	}
	fprintf(f_log, "\n");

	//LOG MEDIUM DATA
	Medium<DIM>* medium = w.get_medium();

	if (medium != nullptr) {
		HomogeneousMedium<DIM, PhaseFunction::Isotropic<DIM>>* m = static_cast<HomogeneousMedium<
				DIM, PhaseFunction::Isotropic<DIM>>*>(medium);

		Vector3 medium_sigma_a = m->get_absorption(VectorN<DIM>(0.)).to_rgb();
		Vector3 medium_sigma_s = m->get_scattering(VectorN<DIM>(0.)).to_rgb();

		fprintf(f_log, "MEDIUM\n");
		fprintf(f_log, "======\n");
		fprintf(f_log, "Absorption: %.10f %.10f %.10f\n", medium_sigma_a[0], medium_sigma_a[1],
				medium_sigma_a[2]);
		fprintf(f_log, "Scattering: %.10f %.10f %.10f\n", medium_sigma_s[0], medium_sigma_s[1],
				medium_sigma_s[2]);
		fprintf(f_log, "PF: g = %.10f\n", p.hg_g);
		switch (p.scat_level) {
			case Scattering::ALL:
				fprintf(f_log, "Scattering level: ALL\n");
				break;
			case Scattering::SINGLE:
				fprintf(f_log, "Scattering level: SINGLE\n");
				break;
			case Scattering::MULTIPLE:
				fprintf(f_log, "Scattering level: MULTI\n");
				break;
			default:
				break;
		}
		fprintf(f_log, "\n");
	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////

	PinholeCamera pinhole_camera = PinholeCamera(p.camera_position, p.camera_looking_at,
			p.camera_up, p.camera_view_plane_dist);
	OrthographicCamera ortho_camera = OrthographicCamera(p.camera_position, p.camera_looking_at,
			p.camera_up, p.camera_w, p.camera_h);

	const Camera3D& camera =
			(p.camera_view_plane_dist == std::numeric_limits<Real>::infinity()) ?
					(const Camera3D&) ortho_camera : (const Camera3D&) pinhole_camera;

	const Camera3D* emitter = &camera;
	if (p.emitter) {
		if (p.emitter_othographic) {
			emitter = new OrthographicCamera(p.emitter_position, p.emitter_looking_at, p.emitter_up,
					p.emitter_w, p.emitter_h);
		} else {
			emitter = new PinholeCamera(p.emitter_position, p.emitter_looking_at, p.emitter_up,
					p.emitter_view_plane_dist);
		}
	}

	printf("Preparing scene to be rendered...\r");

	timer.start();

	// Construct acceleration structures and such
	w.set_background(Spectrum(0.));
	w.freeze();
	{
		Real secs = timer.get_secs();
		int hours = static_cast<int>(secs) / 3600;
		secs -= hours * 3600;
		int minutes = static_cast<int>(secs) / 60;
		secs -= minutes * 60;
		printf("Prepared scene to render: [%d:%d:%f]             \n", hours, minutes, secs);

		if (f_log) {
			fprintf(f_log, "Prepared scene to render: %d:%d:%f\n\n", hours, minutes, secs);
		}
	}

	fprintf(f_log, "SCENE\n");
	fprintf(f_log, "=====\n");

	AABB<DIM> bb = w.get_bounding_volume();
	fprintf(f_log, "BB: [%f, %f, %f] - [%f, %f, %f]\n", bb._min[0], bb._min[1], bb._min[2],
			bb._max[0], bb._max[1], bb._max[2]);
	camera.print(f_log);
	fprintf(f_log, "\n");

	// ----------------------------------------------------------------------
	// Light sources
	if (w.nb_lights() < 1) {
		for (const VectorN<DIM> vpl : p.image_coordinates_vpl) {
			Ray<DIM> r = emitter->get_ray(Vector2(vpl[0], vpl[1]));

			w.add_light(new CosineLightSource<DIM, Radiance>(&w, r, Radiance(vpl[2])));
		}

		if (w.nb_lights() < 1) {
			w.add_light(new PointLightSource<DIM, Radiance>(&w, p.camera_position, Radiance(1.f)));
		}
	}

	fprintf(f_log, "LIGHTS\n");
	fprintf(f_log, "======\n");
	if (f_log) {
		for (unsigned int i = 0; i < w.nb_lights(); ++i) {
			w.light(i)->print(f_log);
		}
	}
	fprintf(f_log, "\n");

	fflush(f_log);

	// ----------------------------------------------------------------------
	// Create Film and Render engine
	std::unique_ptr<FilmS> film = nullptr;
	StratifiedSampler sampler(p.width, p.height, p.sqrt_spp, true);
	
	RenderEngine<DIM, Radiance, RadianceAttenuation>* engine;
	
	/* Select film components */
	FilmComponents comp = FilmComponents::RADIANCE;
	if (p.film_depth) {
		comp = comp | FilmComponents::DEPTH;
	}
	if (p.film_positions) {
		comp = comp | FilmComponents::POSITIONS;
	}
	if (p.film_normals) {
		comp = comp | FilmComponents::NORMALS;
	}

	if (p.transient) {
		// If rendering in transient state...
		// ...create streak-camera film...
		if (!p.is_multibounce_streak) {
			std::unique_ptr<StreakCameraFilmS> sfilm;

			if (!p.film_usereveal) {
				BoxFilter* filter = new BoxFilter(1.);
				if (p.streak_kernel_based) {
					printf("Using a kernel-based streak film...\n");
					if (p.streak_nn > 0) {
						// Kernel bandwidth based on NN search
						sfilm = std::make_unique<KernelBasedStreakCameraFilmS>(p.width, p.height, p.time,
								p.time_resolution, filter, p.streak_nn, p.camera_unwarp,
								p.streak_alpha, p.streak_kernel_radius);
					} else {
						// Fixed initial kernel bandwidth
						sfilm = std::make_unique<KernelBasedStreakCameraFilmS>(p.width, p.height, p.time,
								p.time_resolution, filter, p.streak_kernel_radius, p.camera_unwarp,
								p.streak_alpha);
					}
				} else {
					printf("Using a streak film...\n");
					sfilm = std::make_unique<StreakCameraFilmS>(p.width, p.height, p.time, p.time_resolution,
							p.camera_unwarp, filter, comp);
				}
			}
			sfilm->set_name(p.image_file, "hdr");

			if (p.film_offset > 0.0)
				sfilm->set_offset(p.film_offset);

			film = std::move(sfilm);
		} else {
			std::unique_ptr<MultibounceStreakFilmS> mbfilm;

			if (p.streak_kernel_based) {
				printf("Using a kernel-based multi-bounce film...\n");
				if (p.streak_nn > 0) {
					// Kernel bandwitdh based on NN search
					mbfilm = std::make_unique<MultibounceStreakFilmS>(p.width, p.height, p.time, p.time_resolution,
							p.nb_multibounces, new BoxFilter(1.), p.streak_nn, p.streak_alpha,
							p.streak_kernel_radius, p.camera_unwarp);
				} else {
					// Fixed initial kernel bandwith
					mbfilm = std::make_unique<MultibounceStreakFilmS>(p.width, p.height, p.time, p.time_resolution,
							p.nb_multibounces, new BoxFilter(1.), p.streak_kernel_radius,
							p.streak_alpha, p.camera_unwarp);
				}
			} else {
				printf("Using a multi-bounce film...\n");
				mbfilm = std::make_unique<MultibounceStreakFilmS>(p.width, p.height, p.time, p.time_resolution,
						p.nb_multibounces, new BoxFilter(1.), p.camera_unwarp, comp);
			}

			mbfilm->set_name(p.image_file, "hdr");
			if (p.film_offset > 0.0)
				mbfilm->set_offset(p.film_offset);

			film = std::move(mbfilm);
		}

		// ... and the transient renderer.
		engine = new TransientRenderer<DIM, Radiance, RadianceAttenuation>(f_log, p.time_sampling);
		static_cast<TransientRenderer<DIM, Radiance, RadianceAttenuation>*>(engine)->set_sensor_mode(
				sensor_mode);
	} else {
		printf("Using a steady-state film...\n");
		// If rendering in steady-state...
		// ... create regular film...
		BoxFilter* filter = new BoxFilter(1.);

		film = std::make_unique<FilmS>(p.width, p.height, filter, comp);
		// ... and steady-state renderer.
		engine = new RenderEngine<DIM, Radiance, RadianceAttenuation>(f_log);
	}

	film->set_aspect_ratio(p.aspectRatioY); //Set film aspect ratio

	if (p.single_pixel && p.single_pixel_x > -1 && p.single_pixel_x < p.width
			&& p.single_pixel_y > -1 && p.single_pixel_y < p.height) {
		sampler.set_single_pixel(p.single_pixel_x, p.single_pixel_y);
		if (p.transient) {
			if (!p.is_multibounce_streak) {
				static_cast<StreakCameraFilmS*>(film.get())->set_streak(p.single_pixel_y);
			} else {
				static_cast<MultibounceStreakFilmS*>(film.get())->set_streak(p.single_pixel_y);
			}
		}
	}

	if (p.scanline > -1) {
		sampler.set_scanline(p.scanline);
		if (p.transient) {
			if (!p.is_multibounce_streak) {
				static_cast<StreakCameraFilmS*>(film.get())->set_streak(p.scanline);
			} else {
				static_cast<MultibounceStreakFilmS*>(film.get())->set_streak(p.scanline);
			}
		}
	}
	
	// Reseed the global RNG
	Random::StdRNG.seed(p.new_seed);

	// ----------------------------------------------------------------------
	// Create the integrator to be used
#ifdef _USE_PM_
	typedef PhotonMapping<DIM, Radiance, RadianceAttenuation> BPT;
	BPT bpt(w, p.pm_photon_shots, p.pm_photons_nn, p.pm_kernel_shrinkage, p.vol_path_tracing_samples, f_log, film);
	printf("Rendering with Progressive Photon Mapping...\n");
#else
	typedef BidirectionalPathTracing<DIM, Radiance, RadianceAttenuation> BPT;
	BPT bpt(w, p.vol_path_tracing_samples, f_log, film.get());
	printf("Rendering with Bidirectional Path Tracing...\n");
#endif

	bpt.set_max_path_size(p.max_nb_bounces);
	bpt.set_path_tracing_only(p.path_tracing);
	bpt.set_scattering_components(p.scat_level);

	if (p.transient) {
		if (p.time_sampling) {
			bpt.set_mode(BPT::SamplingMode::Transient);
		} else {
			bpt.set_mode(BPT::SamplingMode::SteadyState);
		}
	}

	fflush(f_log);

	// ----------------------------------------------------------------------
	// GO RENDERING!
	engine->render(p.image_file, w, &bpt, camera, film.get(), &sampler);

	fclose(f_log);
}
