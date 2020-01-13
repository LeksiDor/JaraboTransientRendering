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

/* 
	Created by:			Julio Marco 
	Date of creation:	March 2013
	Description:		SVG parser to extract circles, rectangles, polylines, cameras
						and lights from a SVG file, and to add them to a World object. 
 */
#include "LinearAlgebra/Vector2.h"
#include "RayTracing/World.h"
#include "Camera/LambertianCamera.h"
#include "Geometry/Primitives/2D/Rectangle.h"
#include "Geometry/Primitives/2D/Circle.h"
#include "Geometry/Primitives/2D/Line.h"
#include "LightSource/LightSource.h"
#include <stdio.h>
#include <iostream>
#include <string>

#define LSIZE 8192
#define MAXPOINTS 64
#define MAXOBJECTS 200
#define ENDPATH 10
#define ENDRECT 20
#define CONTINUE 1
#define IGNORE 0

#define DELTA_DIST 0.001

namespace LoaderSVG {
	class SVGParser
	{
		std::string name_file;
		Real scale;
		FILE *f;
		UberMaterial2D *m_color[5];
	private:
		typedef enum {NONE=0, CIRCLE, RECT, POLYLINE, CAMERA, LIGHT} e_object;
		typedef enum {BLACK=0, RED, GREEN, BLUE, WHITE} e_color;
		
		typedef struct{
			float x,y;
		} t_xy;
		
		typedef struct{
			e_object type; 
			char ident[256];
			char label[16]; //for the color
			//Position for circle (center), rectangle (upper-left corner), camera and light
			t_xy pos;
			e_color color;
			bool isLambCam;
			//CIRCLE
			float radius;
			
			//POLYLINE
			t_xy points[MAXPOINTS];
			int npoints;
			
			//RECTANGLE
			t_xy dim; //width (x) and height (y)
		} t_object;
		
		bool in_rect;
		bool in_path;
		
		t_object curr_object;
		t_object vobjects[MAXOBJECTS];
		
		int nobjects;
		
		//AUXILIARY FUNCTIONS
		int parsepoints (char *str, t_xy *vp)
		{ /* parses a "m 87.438332,578.25214 0,-450.79318 L 40.718 -5.98 ..." format string,
		   sets the absolute points xy coords in vp and returns the number of points */
			
			char mode = 'm', *token;
			int np = 0;
			
			token = strtok(str, " ,");
			while (token != NULL){
				if (isalpha(token[0])) //first character is always alphanumeric
					mode = token[0];
				else{
					switch (mode){
						case 'L': //absolute
						case 'M': //absolute
							vp[np].x = atof(token); //x coordinate
							vp[np].y = atof(strtok(NULL, " ,")); //y coordinate (should exist!!)
							np++;
							break;
						case 'l': //relative
						case 'm': //relative
							if (np == 0) //first point is absolute
							{
								vp[np].x = atof(token); //x coordinate
								vp[np].y = atof(strtok(NULL, " ,")); //y coordinate (should exist!!)
								np++;
							}
							else
							{
								vp[np].x = atof(token) + vp[np-1].x; //x coordinate
								vp[np].y = atof(strtok(NULL, " ,"))  + vp[np-1].y; //y coordinate (should exist!!)
								np++;
							}
							break;
						case 'H': //absolute
							vp[np].x = atof(token); //x coordinate
							vp[np].y = vp[np-1].y; //y coordinate (same as previous)
							np++;
							break;
						case 'h': //relative
							vp[np].x = atof(token) + vp[np-1].x; //x coordinate
							vp[np].y = vp[np-1].y; //y coordinate (same as previous)
							np++;
							break;
						case 'V': //absolute
							vp[np].y = atof(token); //x coordinate
							vp[np].x = vp[np-1].x; //x coordinate (same as previous)
							np++;
							break;
						case 'v': //relative
							vp[np].y = atof(token) + vp[np-1].y; //x coordinate
							vp[np].x = vp[np-1].x; //x coordinate (same as previous)
							np++;
							break;
					}
				}
				token = strtok(NULL, " ,");
			}
			return np;
		}
		
		int parseline(char *l)
		{
			char str[LSIZE];
			
			if (in_rect) //<rect
			{
				curr_object.type = RECT;
				if (strstr(l, "width="))
					sscanf(l, "%*[ \t]width=\"%f\"", &curr_object.dim.x);
				else if (strstr(l, "height="))
					sscanf(l, "%*[ \t]height=\"%f\"", &curr_object.dim.y);
				else if (strstr(l, "x="))
					sscanf(l, "%*[ \t]x=\"%f\"", &curr_object.pos.x);
				else if (strstr(l, "y="))
					sscanf(l, "%*[ \t]y=\"%f\"", &curr_object.pos.y);
				else if (strstr(l, "inkscape:label="))
				{
					sscanf(l, "%*[ \t]inkscape:label=\"%[-,. 0-9A-Za-z]\"", curr_object.label);
					if (strcmp("black", curr_object.label) == 0)
						curr_object.color = BLACK;
					else if (strcmp("red", curr_object.label) == 0)
						curr_object.color = RED;
					else if (strcmp("green", curr_object.label) == 0)
						curr_object.color = GREEN;
					else if (strcmp("blue", curr_object.label) == 0)
						curr_object.color = BLUE;
					else curr_object.color = WHITE;
				}
				if (strstr(l, "/>")) return ENDRECT;
				else return CONTINUE;
			}
			else if (in_path) //<path
			{
				if (strstr(l, "sodipodi:type")) //Defining type of path (star, spiral or arc)
				{
					if (strstr(l,"star"))
						curr_object.type = LIGHT;
					//else if (strstr(l,"spiral"))
					//	curr_object.type = CAMERA;
					else if (strstr(l,"arc"))
						curr_object.type = CIRCLE;
				}
				else if (strstr(l, "style=") && curr_object.type == NONE) //THERE'S NO TYPE DEFINED => POLYLINE
					curr_object.type = POLYLINE;
				else //TYPE ALREADY DEFINED
				{
					switch(curr_object.type){
						case CIRCLE:
							if (strstr(l, "sodipodi:cx"))
								sscanf(l, "%*[ \t]sodipodi:cx=\"%f\"", &curr_object.pos.x);
							else if (strstr(l, "sodipodi:cy"))
								sscanf(l, "%*[ \t]sodipodi:cy=\"%f\"", &curr_object.pos.y);
							else if (strstr(l, "sodipodi:rx"))
								sscanf(l, "%*[ \t]sodipodi:rx=\"%f\"", &curr_object.radius);
							else if (strstr(l, "id="))
							{
								sscanf(l, "%*[ \t]id=\"%[-,. 0-9A-Za-z]\"", curr_object.ident);
								curr_object.isLambCam = false;
								if (strcmp("camera", curr_object.ident) == 0)
								{	
									curr_object.isLambCam = true;
									curr_object.color = WHITE;
								}
							}
							else if (strstr(l, "inkscape:label="))
							{
								sscanf(l, "%*[ \t]inkscape:label=\"%[-,. 0-9A-Za-z]\"", curr_object.label);
								if (strcmp("black", curr_object.label) == 0)
									curr_object.color = BLACK;
								else if (strcmp("red", curr_object.label) == 0)
									curr_object.color = RED;
								else if (strcmp("green", curr_object.label) == 0)
									curr_object.color = GREEN;
								else if (strcmp("blue", curr_object.label) == 0)
									curr_object.color = BLUE;
								else curr_object.color = WHITE;
							}
							break;
						case POLYLINE:
							if (strstr(l, "id="))
							{
								sscanf(l, "%*[ \t]id=\"%[-,. 0-9A-Za-z]\"", curr_object.ident);
								curr_object.isLambCam = false;
								if (strcmp("camera", curr_object.ident) == 0)
								{	
									curr_object.isLambCam = true;
									curr_object.color = WHITE;
								}
							}
							else if (strstr(l, "inkscape:label="))
							{
								sscanf(l, "%*[ \t]inkscape:label=\"%[-,. 0-9A-Za-z]\"", curr_object.label);
								if (strcmp("black", curr_object.label) == 0)
									curr_object.color = BLACK;
								else if (strcmp("red", curr_object.label) == 0)
									curr_object.color = RED;
								else if (strcmp("green", curr_object.label) == 0)
									curr_object.color = GREEN;
								else if (strcmp("blue", curr_object.label) == 0)
									curr_object.color = BLUE;
								else curr_object.color = WHITE;
							}
							else if (sscanf(l, "%*[ \t]d=\"%[-,. 0-9A-Za-z]\"", str))
								curr_object.npoints = parsepoints(str, curr_object.points);						
							break;
						case CAMERA:
						case LIGHT:
							if (strstr(l, "sodipodi:cx"))
								sscanf(l, "%*[ \t]sodipodi:cx=\"%f\"", &curr_object.pos.x);
							else if (strstr(l, "sodipodi:cy"))
								sscanf(l, "%*[ \t]sodipodi:cy=\"%f\"", &curr_object.pos.y);
							break;
					}
				}
				if (strstr(l, "/>")) return ENDPATH;
				else return CONTINUE;
			}
			else if (strstr(l, "<path")) {
				in_path = true;
				curr_object.color =	WHITE;
				return CONTINUE;
			}
			else if (strstr(l, "<rect")) {
				in_rect = true;
				curr_object.color =	WHITE;
				return CONTINUE;
			}
			else return IGNORE;
		}
		
		void parsefile()
		/* Parses SVG file and adds the objects to a t_object 1D array */
		{
			char line[LSIZE];
			int ret;
			do{
				fgets(line, LSIZE, f);
				ret = parseline(line);
				switch (ret){
					case ENDPATH:
						vobjects[nobjects++] = curr_object;
						memset(&curr_object, 0, sizeof(t_object));
						curr_object.type = NONE;
						in_path = false; 
						break;
					case ENDRECT:
						vobjects[nobjects++] = curr_object;
						memset(&curr_object, 0, sizeof(t_object));
						curr_object.type = NONE;
						in_rect = false;
						break;
					case CONTINUE: 
						break;
					case IGNORE: 
						break;
				}
			} while(!feof(f));
		}
		
		
		void addtoworld(World<2>* w, LambertianCamera2D& lcam)
		/* Inserts the parsed objects into w */
		{
			std::vector< Object2DPointer > objects;
			//printf("nobjects = %d\n", nobjects);
			for (int i=0; i < nobjects; i++){
				t_object *o = &(vobjects[i]);
				switch (o->type){
					case CIRCLE:
					{	
						//printf("Circle (%f %f - %f)\n", o->pos.x, o->pos.y, o->radius);
						Circle *c = new Circle(scale*Vector2(o->pos.x, o->pos.y), scale*o->radius, m_color[o->color]);
						//ADD CIRCLE
						objects.push_back(c);
						
						if (o->isLambCam)
						{
							//printf("\tCircle as a Lambertian camera\n");
							lcam = LambertianCamera2D(c, DELTA_DIST);
						}
						break;
					}
					case RECT:
					{
#ifndef WIN32
						//printf("Rectangle...\n");
						Rectangle *r = new Rectangle(scale*Vector2(o->pos.x, o->pos.y), 
												 scale*Vector2(o->pos.x+o->dim.x, o->pos.y+o->dim.y), m_color[o->color]);
						//ADD RECTANGLE
						objects.push_back(r);
#endif
						break;

					}
					case POLYLINE:
					{
						for (int j=0; j < o->npoints-1; j++) 
						{
							//printf("Line (%f %f - %f %f)\n", o->points[j].x, o->points[j].y, o->points[j+1].x, o->points[j+1].y);
							Line *l = new Line(scale*Vector2(o->points[j].x, o->points[j].y),
											   scale*Vector2(o->points[j+1].x, o->points[j+1].y), m_color[o->color]);
							//ADD LINE
							objects.push_back(l);
							//If it's a lambertian camera, create the camera associated to it
							if (o->isLambCam)
							{
								//printf("\tLine as a Lambertian camera\n");
								lcam = LambertianCamera2D(l, DELTA_DIST);
							}
						}
						break;
					}
					case CAMERA:
					{	
						//ADD CAMERA
						break;
					}
					case LIGHT:
					{	
						//printf("Light (%f %f)\n", o->pos.x, o->pos.y);
						LightSource<2>* ls = new PointLightSource<2>(w, scale*Vector2(o->pos.x, o->pos.y), Spectrum(1.0));
						w->add_light(ls);
						break;
					}
					default:

						break;
				}
			}
			w->add_objects(objects);
		}
		
	public:
		SVGParser(const std::string &_name_file, float _scale): name_file(_name_file), scale(_scale)
		{
			//printf("svgparser initialized");
			in_rect = false;
			in_path = false;
			nobjects = 0;
			Lambertian2D *rho_white = new Lambertian2D(Spectrum(1.0));
			Lambertian2D *rho_red = new Lambertian2D(Spectrum(1.0, 0.0, 0.0));
			Lambertian2D *rho_green = new Lambertian2D(Spectrum(0.0, 1.0, 0.0));
			Lambertian2D *rho_blue = new Lambertian2D(Spectrum(0.0, 0.0, 1.0));
			Lambertian2D *rho_black = new Lambertian2D(Spectrum(0.0));
			
			for (int i = 0; i < 5; i++) 
				m_color[i] = new UberMaterial2D();
			
			m_color[BLACK]->add_bsdf(rho_black);
			m_color[RED]->add_bsdf(rho_red);
			m_color[GREEN]->add_bsdf(rho_green);
			m_color[BLUE]->add_bsdf(rho_blue);
			m_color[WHITE]->add_bsdf(rho_white);
		}
		
		const int loadSVG(World<2>* w, LambertianCamera2D& lcam)
		{
			curr_object.type = NONE;
			
			f = fopen(name_file.c_str(), "r");
			
			if (!f) return -1;
			
			parsefile();
			fclose(f);
			addtoworld(w, lcam);
			return 1;
		}
	};
	//namespace functions
	void loadfile (const std::string &name_file, World<2> *w, LambertianCamera2D& lcam, float scale = 1.)
	{
		SVGParser parser(name_file, scale);
		if (parser.loadSVG(w, lcam) == -1)
			fprintf(stdout, "ERROR: \"LoaderSVG::load(...)\" failed opening SVG file \"%s\".\n", name_file.c_str());
	}
}

