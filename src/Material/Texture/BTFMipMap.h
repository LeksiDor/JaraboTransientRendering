#ifndef _BTF_MIPMAP_H_
#define _BTF_MIPMAP_H_

#include "Material/Material.h"
#include "BTF.h"
#include "Filter.h"
#include <vector>

/**	Class that implements the BTFs MipMap hierarchy, for each
	domain (spatial+angular). It is templatized as a function 
	of the filter used when building the hierarchy			*/
class BTFMipMap: public Material
{
	Real m_delta_solidangle_view, m_delta_solidangle_light;
	unsigned int m_nb_levels_view, m_nb_levels_light; 
	std::vector<std::vector<BTF*> > m_hierarchy;
	
	void build_hierarchy(BTF *btf);
	void clear();
public:
	BTFMipMap( Medium *default_medium = 0, Real default_n = DEFAULT_REFRACTION_INDEX )
		:Material(default_medium, default_n) {}

	BTFMipMap( Real default_n, Medium *default_medium = 0 )
		:Material(default_medium, default_n) {}

	BTFMipMap( const char *filename, Filter *filter, Medium *default_medium = 0, Real default_n = DEFAULT_REFRACTION_INDEX )
		:Material(default_medium, default_n)
	{	load(filename, filter);	}

	BTFMipMap( const char *filename, Filter *filter, Real default_n, Medium *default_medium = 0 )
		:Material(default_medium, default_n)
	{	load(filename, filter);	}

	BTFMipMap( BTF *btf, Filter *filter )
		:Material(btf->get_dafault_medium(), btf->get_default_refraction_index())
	{	build_hierarchy(btf, filter);	}

	~BTFMipMap()
	{	clear();	}
	
	void load(const char *filename, Filter *filter)
	{	BTF *btf = new BTF(filename); build_hierarchy(btf, filter);	}

	//-----------------------------------------------------------
	//----- Functions inherited from base class 'Material'---------
	virtual Vector3 f(const Intersection &it, const LightSample &ls)const
	{
		// Compute ray differentials

		// Compute solid angles from ray differentials

		// Get hierarchy level
		
		// Shade

	}
	virtual Vector3 f(const Vector3 &omega_o, const Vector3 &omega_i, const Vector3 &normal, const Vector2 &uv) const
	{
		return get_reflectance(uv, omega_i, omega_o,normal);
	}
	virtual Vector3 sample_direction(const Vector3 &omega_o, Vector3 &omega_i, const Vector3 &normal, const Vector2 &uv, Real &pdf) const
	{
		Vector3 wi;
		Sampling::cosine_weighted(wi, pdf);
		omega_i = wi.transform_matrix_to(Vector3(0,1,0), normal);

		return f(omega_o, omega_i, normal, uv);
	}
	virtual Vector3 sample_outgoing_ray(const Intersection &it, Ray &new_ray, Real &pdf) const
	{
		Vector3 wi, omega_i; 
		Sampling::cosine_weighted(wi, pdf);
		omega_i = wi.transform_matrix_to(Vector3(0,1,0), it.get_normal());
		new_ray = Ray(it.get_position(), omega_i,
			false, it.get_ray().get_level()+1, it.get_ray().get_refraction_index(), 
			it.get_ray().get_medium());

		return f(it.get_ray().get_direction(), omega_i, it.get_normal(), it.get_uv());
		
	}
}; //BTFMipMap

void BTFMipMap::build_hierarchy(BTF *btf, Filter *filter)
{
	m_hierarchy.clear();
	std::vector<BTF*> v;
	m_hierarchy.push_back(v);
	
	//Add finest level (non-mipmapped BTF)
	m_hierarchy.back().push_back(btf);	

	//While BTF can be reduced...
	BTF *nbtf = btf;
	while( nbtf->get_width() > 1 && nbtf->get_height() > 1 )
	{
		// ... create a new level in the hierarchy...
		m_hierarchy.push_back(v);

		// ... and create mipmapped version for each solid angle.
		for( int v=0; v < m_nb_levels_view; ++v )
		for( int l=0; l < m_nb_levels_light; ++l )
		{
			// Compute directions where we are evaluating light...
			Real solid_angle_v = static_cast<Real>(v)*m_delta_solidangle_view;
			Real solid_angle_l = static_cast<Real>(l)*m_delta_solidangle_view;
			// ... and use them to mipmap the BTF (and store it).
			m_hierarchy.back().push_back(nbtf->mipmap(filter, solid_angle_v, solid_angle_l));
		}
		// Finally, update the BTF to downsample
		nbtf = m_hierarchy.back().front();
	}
}

#endif //_BTF_MIPMAP_H_