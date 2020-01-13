#ifndef _DOF_CAMERA_H_
#define _DOF_CAMERA_H_

#include "bunnykiller.h"

#include "Camera/Camera.h"
#include "Utils/RandomNumbers.h"

class FiniteApertureCamera: public PinholeCamera
{
	Real r_apperture;
	Real focus_distance;
	// Gives a (x,y) position on the camera apperture
	Vector3 sample_aperture()const
	{
		float PI = 3.14159265;
		float Phi = r_apperture*sqrt(RNG::StdRandom::get_real());
		float Theta = 2*PI*RNG::StdRandom::get_real();

		float x = Phi*cos(Theta);
		float y = Phi*sin(Theta);

		return Vector3(x, y, 0.);
	}
public:
	FiniteApertureCamera(const Vector3& c, const Vector3& f, const Vector3& u,
		const Real fd, const Real _r_apperture, const Real _focus_distance) :
		PinholeCamera(c,f,u,fd), r_apperture(_r_apperture), focus_distance(_focus_distance)
	{}
	
	virtual ~FiniteApertureCamera()
	{}

	/// Calculates the ray 
	virtual Ray3D get_ray(const Vector2& p) const
	{
		//We first sample the apperture
		Vector3 point_lens = sample_aperture();
		//And translate it to camera coordinates
		point_lens = m_uv*point_lens[1]+m_vv*point_lens[0];

		//We now sample the point in the focus plane
		Vector3 pvp(m_nv);
		/*pvp *= focal_dist; 
		pvp += uv*p[0];
		pvp += vv*p[1];*/

		pvp *= focus_distance;
		pvp += m_uv*p[0]*(focus_distance/m_focal_dist);
		pvp += m_vv*p[1]*(focus_distance/m_focal_dist);
		
		

		//We now get the direction of the ray
		Vector3 dir = (pvp - point_lens).normalized();

		return Ray3D(m_position+point_lens, dir);
	}
}; //FiniteApertureCamera
#endif //_DOF_CAMERA_H_
