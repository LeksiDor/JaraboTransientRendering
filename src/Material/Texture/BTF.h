#ifndef _BTF_H_
#define _BTF_H_

#include "Material/Material.h"
#include "Texture.h"
#include "Sampling/Sampling.h"
#include "Filter.h"

/** Class that implements a Material defined by a BTF.
	It assumes that the BTF is in a uncompressed format. */
class BTF: public Material, Texture
{
	// File containing the BTF (out-of-core)	
	char m_name[1024];
	FILE *m_fp;
	int m_data_offset;

	
	// Directional information of the BTF including...
	//	... angular dimensions...
	int m_nb_v, m_nb_l;
	//	... sampling directions...
	vector<Vector3> directions_v, directions_l;
	//	... info for nearest direction sampling...
	vector<vector<Real> > rows_v, rows_l;
	vector<vector<int> > row_to_index_v, row_to_index_l;
	//	... and min sampled cosine. 
	float min_cos_v, min_cos_l;

	bool b_float;
	
	// Load function
	void load_from_BTF(const char *filename);
	
	// Data access functions
	Vector3 get_reflectance(const Vector2 &x, const Vector3 &v, const Vector3 &l, const Vector3 &n)const;
	
	// Auxiliar functions to access data
	void get_surface_offset(const Vector2 &c, std::vector<int> &offset_coordinates, std::vector<float> &weights_x)const;
	void get_closest_directions(const Vector3 &v, 
		const vector<vector<Real> > &rows, const vector<vector<int> > &row_to_index, 
		const vector<Vector3> directions, const Real min_cos,
		std::vector<int> &dir_indices, std::vector<Real> &dir_weights)const;
	
	void get_raw_data(const int x_offset, void *data)const;
	Vector3 get_data(const int x_offset, const int v_offset, const int l_offset)const;

	// Function to clear BTF data
	void clear();
public:
	BTF( )
		:Material(),m_fp(0), m_data_offset(-1), 
		 min_cos_v(1.), min_cos_l(1.),m_nb_v(0), m_nb_l(0){}
	BTF( const char *filename )
		:Material(),m_fp(0), m_data_offset(-1), 
		 min_cos_v(1.), min_cos_l(1.),m_nb_v(0), m_nb_l(0)
	{	load(filename); }
	~BTF()
	{
		clear();
	}
	void load(const char *filename);

	// Duplicate the BTF in file 'filename' and returns it
	BTF *duplicate(const char *filename)const;


	//-----------------------------------------------------------
	//----- Functions inherited from base class 'Texture'--------
	//Resize (and potentionally reshape) the texture
	template<class Filter>
	void resize(int w, int h)
	{
		//Need to be done!
	}
	
	// Access data
	virtual const Vector3 &operator()(Real x, Real y)const
	{
		return f(Vector3(0,0,1),Vector3(0,0,1), Vector3(0,0,1), Vector2(x,y));
	}
	virtual const Vector3 &operator()(int x, int y)const
	{
		return f(Vector3(0,0,1),Vector3(0,0,1), Vector3(0,0,1), 
				 Vector2(static_cast<Real>(x)*m_inv_width, static_cast<Real>(y)*m_inv_height));
	}
	
	//-----------------------------------------------------------
	//----- Functions inherited from base class 'Material'-------
	virtual Vector3 f(const Vector3 &omega_o, const Vector3 &omega_i, const Vector3 &normal, const Vector2 &uv) const
	{
		if( !m_fp )
			throw( "No BTF loaded! Thus, you cannot shade!");

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
	
	//-----------------------------------------------------------
	//----- Mipmap function -------------------------------------
	BTF* mipmap(Filter *filter, Real solid_angle_v, Real solid_angle_l)const;

	/** DEPRECATED: Shade at intersection it */
	Vector3 shade(const Intersection &it)const
	{
		if(!m_fp )
			throw( "No BTF loaded! Thus, you cannot shade!");
		
		return Vector3();
	}

}; //BTF


void BTF::clear()
{
#ifdef X64
	if( data4 || data1 )
	{
		if(data4) delete[] data4;
		if(data1) delete[] data1;
		data4 = 0; data1 = 0;
		
		delete[] vec_angle_v;
		delete[] vec_angle_l;
		
		m_data_offset = -1;
		vec_angle_v = vec_angle_l = 0; 
		m_width = m_height = m_nb_v = m_nb_l = element_size = 0;
		b_float = false;
	}
#else
	if( m_fp ) //It means one BTF has been loaded...
	{
		fclose(m_fp);
		m_fp = 0;

		m_data_offset = -1;
		m_width = m_height = m_nb_v = m_nb_l = 0;
		b_float = false;

	}
#endif
}
void BTF::load(const char *filename)
{
	clear();
	
	char *ext;
	char aux_filename[512];
	strcpy(aux_filename, filename);
	ext = strtok(aux_filename, ".");
	strcpy(m_name, aux_filename);
	ext = strtok(NULL, ".");

	if (strcmp(ext, "btf_byte") == 0)
	{
		b_float = false;
		m_size_element = SizeElement::Byte;
	} 
	else if (strcmp(ext, "btf_float") == 0)
	{
		b_float = true;
		m_size_element = SizeElement::Float;
	}
	else throw("Unknown file extension on BTF loading");

	load_from_BTF(filename);

	//print();
}

void BTF::load_from_BTF(const char *filename)
{
	if(!(m_fp = fopen(filename, "rb")))
		throw("Impossible to open file on BTF loading");

	fread(&m_data_offset, sizeof(m_data_offset), 1, m_fp);
	fread(&m_width, sizeof(m_width), 1, m_fp);
	fread(&m_height, sizeof(m_height), 1, m_fp);
	m_inv_width = 1./static_cast<Real>(m_width);
	m_inv_height = 1./static_cast<Real>(m_height);
	
	vector<Vector2> dir;

	min_cos_v = 1.;
	fread(&m_nb_v, sizeof(m_nb_v), 1, m_fp);
	for (int i = 0; i < m_nb_v; i++)
	{
		float t, p;
		fread(&t, sizeof(t), 1, m_fp);
		fread(&p, sizeof(p), 1, m_fp);

		t = t/180.0f*M_PI;
		p = p/180.0f*M_PI;

		float cos_t = cos(t);
		min_cos_v = (cos_t<min_cos_v)?cos_t:min_cos_v;

		dir.push_back(Vector2(t,p));
		directions_v.push_back(Vector3(sin(t)*cos(p),sin(t)*sin(p),cos_t));
	}
	
	if( !dir.empty() )
	{
		vector<Real> v;
		rows_v.push_back(v);
		rows_v[0].push_back(dir[0][0]);
		rows_v[0].push_back(dir[0][1]);
		vector<int> vi;
		row_to_index_v.push_back(vi);
		row_to_index_v[0].push_back(0);
		for( int i=1; i<dir.size(); ++i)
		{
			if( dir[i][0] > rows_v.back()[0] )
			{
				vector<Real> v;
				rows_v.push_back(v);
				rows_v.back().push_back(dir[i][0]);
				vector<int> vi;
				row_to_index_v.push_back(vi);
			}
			rows_v.back().push_back(dir[i][1]);
			row_to_index_v.back().push_back(i);

		}
		for( int i=0; i<rows_v.size(); ++i)
		{
			if( rows_v[i].size() == 2 )
			{
				for(int j=2; j < rows_v[i+1].size(); ++j )
				{	
					rows_v[i].push_back( rows_v[i+1][j]);
					row_to_index_v[i].push_back(row_to_index_v[i].back());
				}
			}
			rows_v[i].insert(rows_v[i].begin()+1, rows_v[i].back()-2.*M_PI);
			row_to_index_v[i].insert(row_to_index_v[i].begin(), row_to_index_v[i].back());
			rows_v[i].push_back(rows_v[i][2]+2.*M_PI);
			row_to_index_v[i].push_back(row_to_index_v[i][1]);
	}	
	}

	dir.clear();

	min_cos_l = 1.;
	fread(&m_nb_l, sizeof(m_nb_l), 1, m_fp);
	for (int i = 0; i < m_nb_l; i++)
	{
		float t, p;
		fread(&t, sizeof(t), 1, m_fp);
		fread(&p, sizeof(p), 1, m_fp);

		t = t/180.0f*M_PI;
		p = p/180.0f*M_PI;

		float cos_t = cos(t);
		min_cos_l = (cos_t<min_cos_l)?cos_t:min_cos_l;

		dir.push_back(Vector2(t,p));
		directions_l.push_back(Vector3(sin(t)*cos(p),sin(t)*sin(p),cos_t));
	}

	if( !dir.empty() )
	{
		vector<Real> v;
		rows_l.push_back(v);
		rows_l[0].push_back(dir[0][0]);
		rows_l[0].push_back(dir[0][1]);
		
		vector<int> vi;
		row_to_index_l.push_back(vi);
		row_to_index_l[0].push_back(0);
		for( int i=1; i<dir.size(); ++i)
		{
			if( dir[i][0] > rows_l.back()[0] )
			{
				vector<Real> v;
				rows_l.push_back(v);
				rows_l.back().push_back(dir[i][0]);
				vector<int> vi;
				row_to_index_l.push_back(vi);
			}
			rows_l.back().push_back(dir[i][1]);
			row_to_index_l.back().push_back(i);

		}
		for( int i=0; i<rows_l.size(); ++i)
		{
			if( rows_l[i].size() == 2 )
			{
				for(int j=2; j < rows_l[i+1].size(); ++j )
				{
					rows_l[i].push_back( rows_l[i+1][j]);
					row_to_index_l[i].push_back(row_to_index_l[i].back());
				}
			}
			rows_l[i].insert(rows_l[i].begin()+1, rows_l[i].back()-2.*M_PI);
			row_to_index_l[i].insert(row_to_index_l[i].begin(), row_to_index_l[i].back());
			rows_l[i].push_back(rows_l[i][2]+2.*M_PI);
			row_to_index_l[i].push_back(row_to_index_l[i][1]);
		}
	}
}
void BTF::get_surface_offset(const Vector2 &c, std::vector<int> &offset_coordinates, std::vector<float> &weights_x)const
{
	float c_x = m_width * c[0]; 
	float c_y = m_height * c[1];
	
	if( m_spatial_filtering == Filtering::NearestNeighbor )
	{
		int t_x = static_cast<int>(c_x+((c_x>0.)?.5:-.5))%(m_width);
		int t_y = static_cast<int>(c_y+((c_y>0.)?.5:-.5))%(m_height);

		t_x = (t_x<0)?m_width+t_x:t_x;
		t_y = (t_y<0)?m_height+t_y:t_y;

		offset_coordinates.push_back(t_y*m_width+t_x);
		weights_x.push_back(1.);
	}
	else if( m_spatial_filtering == Filtering::Bilinear )
	{

		int t_x = static_cast<int>(c_x);
		int t_y = static_cast<int>(c_y);

		float u = fabs(c_x-static_cast<float>(t_x));
		float v = fabs(c_y-static_cast<float>(t_y));

		t_x = t_x%m_width; t_y = t_y%m_height;

		t_x = (t_x<0)?m_width+t_x:t_x;
		t_y = (t_y<0)?m_height+t_y:t_y;


		
		offset_coordinates.push_back(t_y*m_width+t_x);
		weights_x.push_back((1.-u)*(1.-v));
		offset_coordinates.push_back(t_y*m_width+(t_x+1)%m_width);
		weights_x.push_back(u*(1.-v));
		offset_coordinates.push_back((t_y+1)%m_height*m_width+t_x);
		weights_x.push_back(v*(1.-u));
		offset_coordinates.push_back((t_y+1)%m_height*m_width+(t_x+1)%m_width);
		weights_x.push_back(u*v);
	}
	else
		throw( "Non supported filtering on BTFs");

}

void BTF::get_closest_directions(const Vector3 &v, 
	const vector<vector<Real> > &rows, const vector<vector<int> > &row_to_index, 
	const vector<Vector3> directions, const Real min_cos,
	std::vector<int> &dir_indices, std::vector<Real> &dir_weights)const
{
	dir_indices.clear(); dir_weights.clear();

	Real phi = atan2(v[1], v[0]); phi = phi<0.?(2.*M_PI+phi):phi;
	if( v[2]>min_cos )
	{
		Real theta = acos(v[2]);

		int row = 0;
		while( rows[row+1][0] < theta )
			++row;
		
		vector<int> indices;
		vector<Real>::const_iterator up;
		
		up = upper_bound(rows[row].begin()+1, rows[row].end(), phi);
		int upper_r = static_cast<int>(up-rows[row].begin())-1;
		int lower_r = upper_r - 1; 
		indices.push_back(row_to_index[row][lower_r]);
		indices.push_back(row_to_index[row][upper_r]);

		up = upper_bound(rows[row+1].begin()+1, rows[row+1].end(), phi);
		int upper_r1 = static_cast<int>(up-rows[row+1].begin())-1;
		int lower_r1 = upper_r1 - 1;
		indices.push_back(row_to_index[row+1][lower_r1]);
		indices.push_back(row_to_index[row+1][upper_r1]);
	
		int idx_1 = 2;
		Vector3 v_(theta, phi,0.), 
			v1(rows[row][0], rows[row][lower_r+1], 0.), 
			v2(rows[row+1][0], rows[row+1][lower_r1+1], 0.), 
			v3(rows[row+1][0], rows[row+1][upper_r1+1], 0.);

		Real cu = ( (v2[1]-v3[1])*(v_[0]-v3[0])+(v3[0]-v2[0])*(v_[1]-v3[1]) )/( (v2[1]-v3[1])*(v1[0]-v3[0])+(v3[0]-v2[0])*(v1[1]-v3[1]) );
		Real cv = ( (v3[1]-v1[1])*(v_[0]-v3[0])+(v1[0]-v3[0])*(v_[1]-v3[1]) )/( (v2[1]-v3[1])*(v1[0]-v3[0])+(v3[0]-v2[0])*(v1[1]-v3[1]) );

		if( cu < 0. || cv < 0. || (cu+cv)> 1.)
		{
			idx_1 = 1;
			v2 = Vector3(rows[row][0], rows[row][upper_r+1], 0.); 
			cu = ( (v2[1]-v3[1])*(v_[0]-v3[0])+(v3[0]-v2[0])*(v_[1]-v3[1]) )/( (v2[1]-v3[1])*(v1[0]-v3[0])+(v3[0]-v2[0])*(v1[1]-v3[1]) );
			cv = ( (v3[1]-v1[1])*(v_[0]-v3[0])+(v1[0]-v3[0])*(v_[1]-v3[1]) )/( (v2[1]-v3[1])*(v1[0]-v3[0])+(v3[0]-v2[0])*(v1[1]-v3[1]) );
		}
		
		if( cu < 0. || cv < 0. || (cu+cv)> 1.)
		{
			cu = 0.; cv = 0.;
		}
		dir_indices.push_back(indices[0]);
		dir_indices.push_back(indices[idx_1]);
		dir_indices.push_back(indices[3]);
	
		dir_weights.push_back(cu);
		dir_weights.push_back(cv);
		dir_weights.push_back(1.-cu-cv);
	}
	else
	{
		vector<Real>::const_iterator up;
		up = upper_bound(rows.back().begin()+1, rows.back().end(), phi);
		int upper = static_cast<int>(up-rows.back().begin())-1;
		int lower = upper - 1;
		upper = upper %(rows.back().size()-1);
		dir_indices.push_back(row_to_index.back()[lower]);
		dir_indices.push_back(row_to_index.back()[upper]);

		Real dist1 = (directions[dir_indices[0]]-v).length();
		Real dist2 = (directions[dir_indices[1]]-v).length();
		Real inv_sum = 1./(dist1+dist2);
		
		dir_weights.push_back(dist2*inv_sum);
		dir_weights.push_back(dist1*inv_sum);
	}
}

void BTF::get_raw_data(const int x_offset, void *data)const
{
	
#ifdef X64
	if( b_float )
		data = &(data4[__int64(x_offset*m_nb_l*m_nb_v)*3*sizeof(float)]);
	else
		data = &(data1[__int64(x_offset*m_nb_l*m_nb_v)*3*sizeof(unsigned char)]);
#else
	if( b_float )
	{
		__int64 offset = m_data_offset+__int64(x_offset*m_nb_l*m_nb_v)*3*sizeof(float);
		_fseeki64(m_fp, offset, SEEK_SET);
		fread(data, sizeof(float), 3*m_nb_l*m_nb_v, m_fp);		
	}
	else
	{
		__int64 offset = m_data_offset+__int64(x_offset*m_nb_l*m_nb_v)*3*sizeof(unsigned char);
		_fseeki64(m_fp, offset, SEEK_SET);
		fread(data, sizeof(unsigned char), 3*m_nb_l*m_nb_v, m_fp);
	}
#endif
}

Vector3 BTF::get_data(const int x_offset, const int v_offset, const int l_offset)const
{
#ifdef X64
	if( b_float )
	{
		__int64 offset = __int64(x_offset*m_nb_l*m_nb_v + v_offset*m_nb_l + l_offset)*3*sizeof(float);
		float *data = &(data4[offset]);
		return Vector3(data[0], data[1], data[2]);
	}
	else
	{
		__int64 offset = __int64(x_offset*m_nb_l*m_nb_v + v_offset*m_nb_l + l_offset)*3*sizeof(unsigned char);
		unsigned char *data = &(data1[offset]);
		return Vector3(static_cast<float>(data[0])/255., 
					 static_cast<float>(data[1])/255.,
					 static_cast<float>(data[2])/255.);
	}
#else
	if( b_float )
	{
		__int64 offset = m_data_offset+__int64(x_offset*m_nb_l*m_nb_v + v_offset*m_nb_l + l_offset)*3*sizeof(float);
		_fseeki64(m_fp, offset, SEEK_SET);
		float data[3];
		fread(data, sizeof(float), 3, m_fp);
		return Vector3(data[0], data[1], data[2]);
		
	}
	else
	{
		__int64 offset = m_data_offset+__int64(x_offset*m_nb_l*m_nb_v + v_offset*m_nb_l + l_offset)*3*sizeof(unsigned char);
		_fseeki64(m_fp, offset, SEEK_SET);
		unsigned char data[3];
		fread(data, sizeof(unsigned char), 3, m_fp);
		return Vector3(static_cast<float>(data[0])/255., 
							 static_cast<float>(data[1])/255.,
							 static_cast<float>(data[2])/255.);
	}
#endif
}
Vector3 BTF::get_reflectance(const Vector2 &tcoor, const Vector3 &view, const Vector3 &light, const Vector3 &normal)const
{
	Vector3 reflected;

	vector<int> offset_coordinates;
	vector<float> weights_x;
	
	//Get offset coordinates for the sample
	get_surface_offset(tcoor, offset_coordinates, weights_x);
	Vector3 v = (-view).transform_matrix_to_z(normal, Vector3(0,-1,0));
	
	if( v[2] < 0.)
		return Vector3();

	vector<int> idx_directions_v;
	vector<float> weights_v;
	get_closest_directions(v, rows_v, row_to_index_v, directions_v, min_cos_v, idx_directions_v, weights_v);
	//idx_directions_v.push_back(0);
	//weights_v.push_back(1.);

	//Get light directions
	Vector3 l = (-light).transform_matrix_to_z(normal, Vector3(0,-1,0));
	if( l[2] < 0.)
		return Vector3();

	vector<int> idx_directions_l;
	vector<float> weights_l;
	get_closest_directions(l, rows_l, row_to_index_l, directions_l, min_cos_l, idx_directions_l, weights_l);
	
	//return Vector3(idx_directions_l[0],idx_directions_l[1],(idx_directions_l.size()>2)?idx_directions_l[2]:0.);
	//idx_directions_l.push_back(0);
	//weights_l.push_back(1.);

	//Iterate over all texels that contribute...
	for( int sx = 0; sx < offset_coordinates.size(); ++sx )
		for( int sv = 0; sv < idx_directions_v.size(); ++sv )
			for( int sl = 0; sl < idx_directions_l.size(); ++sl )
				reflected += get_data(offset_coordinates[sx], idx_directions_v[sv], idx_directions_l[sl])
					*weights_x[sx]*weights_v[sv]*weights_l[sl];
				
	return reflected*((l[2] < min_cos_l)?(l[2]/min_cos_l):1.);
}





BTF* BTF::mipmap(Filter *filter, Real solid_angle_v, Real solid_angle_l)const
{
	return 0;
}
#endif //_BTF_H_