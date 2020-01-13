/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef _MESH_H_
#define _MESH_H_

#include "Geometry/Object.h"
#include "Geometry/Primitives/3D/Triangle.h"
#include "LinearAlgebra/MatrixN.h"
#include "Material/Material.h"
#include "RayTracing/Intersection.h"
#include "RayTracing/Ray.h"
#include <vector>
#include <string>
#include <string>
#include <fstream>
#include <stdio.h>
#include <iostream>

class Mesh: public Object3D
{
	std::vector<Triangle> triangles;
	AABB3D bb;
public:
 
	Mesh(const std::vector<Triangle> &_triangles, Material3D *_m)
		:Object3D(_m), triangles(_triangles)
	{
		for (unsigned int i=0; i<triangles.size(); ++i)
			bb.add(triangles[i].get_bounding_box());
	}
	Mesh(const std::string &name_file, Material3D *_m)
		:Object3D(_m)
	{
		load(name_file);
		for (unsigned int i=0; i<triangles.size(); ++i)
			bb.add(triangles[i].get_bounding_box());

	}
	virtual void transform(const Matrix4x4& m)
	{
	}
  
	/// Intersect the object with a ray.
	bool intersect(Ray3D& r, Intersection3D &it, float max_distance ) const
	{
		Real t;
		bool result (false);
		
		if( bb.intersect(r,t) && t<r.get_parameter())
		{
			for( unsigned int i=0; i<triangles.size(); ++i)
			{
				result = result | triangles[i].intersect(r,it,max_distance);
			}
		}
		return result;
	}
	bool intersect(const Ray3D& r, float max_distance ) const
	{
		bool result (false);
		Real t;
		if( bb.intersect(r,t) && t<r.get_parameter())
			for( unsigned int i=0; i<triangles.size() && !result; ++i)
			{
				result = result | triangles[i].intersect(r,max_distance);
			}
		
		return result;
	}

	AABB3D get_bounding_box()const
	{
		return bb;
	}

	Vector3 get_center()const
	{
		return bb.get_center();
	}
	
	void push_primitives(std::vector<Object3D*> &objects)
	{
		for(unsigned int i=0; i<triangles.size(); ++i)
			objects.push_back(&triangles[i]);
	}
	void push_primitives(std::vector<Triangle> &objects)
	{
		for(unsigned int i=0; i<triangles.size(); ++i)
			objects.push_back(triangles[i]);
	}
	Real get_intersection_cost()const
	{
		if( triangles.size() )
			return triangles.size()*triangles[0].get_intersection_cost();
		else
			return 0;
	}
	
	//Loads the mesh from an obj file
	void load(const std::string &name_file );
}; // Mesh

void Mesh::load(const std::string &name_file )
{
	std::ifstream fin(name_file.c_str());
	if(!fin.is_open()) printf("Error: failed to open model file.\n");
	
	char c;
	Real x, y, z;
	int Va, Vb, Vc;
	int Na, Nb, Nc;
	int Ta, Tb, Tc;
	std::string identifier;
	
	std::vector<Vector3> vertices;
	std::vector<Vector3> normals;
	std::vector<Vector2> texcoords;

	for(fin >> identifier; !fin.eof(); fin >> identifier)
	{
		//Store vertices
		if(!identifier.compare("v"))
		{
			fin >> x >> y >> z;
			vertices.push_back(Vector3(x, y, z));
		}
		
		//Store normals
		else if(!identifier.compare("vn"))
		{ 
			fin >> x >> y >> z;
			normals.push_back(Vector3(x, y, z));
		}
		
		//Store texcoords
		else if(!identifier.compare("vt"))
		{
			//by now, we ignore texcoords
			fin >> x >> y;
			texcoords.push_back(Vector2(x, y));
		}
		
		//Store triangles
		else if(!identifier.compare("f"))
		{
			//If only vertices - f v1 v2 v3
			if(normals.empty() && texcoords.empty())
			{
				fin >> Va >> Vb >> Vc;
				triangles.push_back(Triangle(vertices[--Va], vertices[--Vb], vertices[--Vc], mat));
			}

			//If vertices and normals - f v1//vn1 v2//vn2 v3//vn3
			else if(!normals.empty() && texcoords.empty())
			{
				fin >> Va >> c >> c >> Na
					>> Vb >> c >> c >> Nb
					>> Vc >> c >> c >> Nc;
				triangles.push_back(Triangle(vertices[--Va], normals[--Na], vertices[--Vb], normals[--Nb], vertices[--Vc], normals[--Nc], mat));
			}

			//If vertices and texcoords and normals - f v1/vt1/vn1 v2/vt2/vn2 v3/vt3/vn3
			else if(!normals.empty() && !texcoords.empty())
			{
				fin >> Va >> c >> Ta >> c >> Na
					>> Vb >> c >> Tb >> c >> Nb
					>> Vc >> c >> Tc >> c >> Nc;
				triangles.push_back(Triangle(vertices[--Va], normals[--Na], texcoords[--Ta], vertices[--Vb], normals[--Nb], texcoords[--Tb], vertices[--Vc], normals[--Nc], texcoords[--Tc], mat));
				//triangles.push_back(Triangle(vertices[--Va], normals[--Na], vertices[--Vb], normals[--Nb], vertices[--Vc], normals[--Nc], mat));

			}
			
			else
			{
				printf("Error: unsupported mesh format.\n");
			}
		}
	}
	fin.close();
}


#endif //_MESH_H_