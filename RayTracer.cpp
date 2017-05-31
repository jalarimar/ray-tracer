/*========================================================================
* COSC 363  Computer Graphics (2017)
* Ray tracer
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include <math.h>
#include "Sphere.h"
#include "Cylinder.h"
#include "Cone.h"
#include "Plane.h"
#include "SceneObject.h"
#include "Ray.h"
#include "PerlinNoise.h"
#include "TextureBMP.h"
#include <GL/glut.h>
using namespace std;

const float WIDTH = 20.0;  
const float HEIGHT = 20.0;
const float EDIST = 40.0;
const int NUMDIV = 500;
const int MAX_STEPS = 5;
const float XMIN = -WIDTH * 0.5;
const float XMAX =  WIDTH * 0.5;
const float YMIN = -HEIGHT * 0.5;
const float YMAX =  HEIGHT * 0.5;
const float PI = 3.14159;

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene
TextureBMP texturePurple;
TextureBMP textureBlue;
TextureBMP textureGreen;
PerlinNoise perlin;

glm::vec3 textureBackground(Ray ray) {
	float texcoords = (ray.xpt.x + 40) / 80;
	float texcoordt = (ray.xpt.y + 20) / 80;
	return textureBlue.getColorAt(texcoords, texcoordt);
}

glm::vec3 textureFloor(Ray ray) {
	float texcoords = (ray.xpt.x + 40) / 100;
	float texcoordt = (ray.xpt.z + 200 + 40) / 200;
	return textureBlue.getColorAt(texcoords, 1 - texcoordt);
}

glm::vec3 textureLeftside(Ray ray) {
	float texcoords = (ray.xpt.z + 200) / 200;
	float texcoordt = (ray.xpt.y + 20 + 10) / 140;
	return textureBlue.getColorAt(texcoords, 1 - texcoordt);
}

glm::vec3 textureRightside(Ray ray) {
	float texcoords = (ray.xpt.z + 200) / 200;
	float texcoordt = (ray.xpt.y + 20) / 80;
	return textureBlue.getColorAt(1 - texcoords, texcoordt);
}

glm::vec3 textureTopside(Ray ray) {
	float texcoords = (ray.xpt.x + 40) / 100;
	float texcoordt = (ray.xpt.z + 240) / 200;
	return textureBlue.getColorAt(texcoords, 1 - texcoordt);
}

glm::vec3 textureBehindside(Ray ray) {
	float texcoords = (ray.xpt.x + 80) / 160;
	float texcoordt = (ray.xpt.y + 40) / 160;
	return textureBlue.getColorAt(texcoords, texcoordt);
}

glm::vec3 texture2D(Ray ray) {
	float norx = (ray.xpt.x + 18) / 8; // normalise x
	float nory = (ray.xpt.y + 10) / 8; // normalise y
	float norz = (ray.xpt.z + 115) / 8; // normalise z
	float texcoords = 0.5 + atan2(norz, norx) / (2 * PI);
	float texcoordt = 0.5 - asin(nory) / PI;
	return texturePurple.getColorAt(texcoords, texcoordt);
}

glm::vec3 texture2DGreen(Ray ray) {
	float norx = (ray.xpt.x + 5) / 15; // normalise x
	float nory = (ray.xpt.y + 5) / 15; // normalise y
	float norz = (ray.xpt.z + 145) / 15; // normalise z
	float texcoords = 0.5 + atan2(norz, norx) / (2 * PI);
	float texcoordt = 0.5 - asin(nory) / PI;
	return textureGreen.getColorAt(texcoords, texcoordt);
}

glm::vec3 textureProcedural(Ray ray) {
	float norx = (ray.xpt.x - 2.5) / 16; // normalise x
	float nory = (ray.xpt.y - 5) / 20; // normalise y
	double noise = perlin.noise(nory * 5, norx * 5, 0);
	return glm::vec3(1, (1 + sin(30 * nory + 20 * noise)) / 2.5, 0); 
	// decreasing 20 increases smokiness, decreasing 2.5 decreases blending
}

glm::vec3 textureProceduralBox(Ray ray) {
	float norx = (ray.xpt.x + 20) / 10; // normalise x
	float nory = (ray.xpt.y + 10) / 20; // normalise y
	float norz = (ray.xpt.z + 110) / 20; // normalise z	
	double noise = perlin.noise(nory * 5, norx * 5, norz);
	return glm::vec3(1 + sin(10 * nory + 20 * noise), (1 + sin(30 * nory + 20 * noise)) / 2.5, 1 + sin(30 * norx + 20 * noise));
}


//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step) // step starts at 1
{
	glm::vec3 backgroundCol(0);
	glm::vec3 ambientCol(0.1);
	glm::vec3 light(10, 40, -3);

    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour

    glm::vec3 col = sceneObjects[ray.xindex]->getColor(); // object's colour
    
	glm::vec3 normalVector = sceneObjects[ray.xindex]->normal(ray.xpt);
	glm::vec3 lightVector = light - ray.xpt;
	glm::vec3 unitLightVector = glm::normalize(lightVector); // normalize
	float lDotn = glm::dot(unitLightVector, normalVector); // dot product    
    glm::vec3 reflVector = glm::reflect(-unitLightVector, normalVector);
    float rDotn = glm::dot(reflVector, normalVector);
    
    glm::vec3 specularCol;
    if (rDotn < 0) {
		specularCol = glm::vec3(0);
	} else {
		float term = pow(rDotn, 40); // (r.n)^f where f is Phong's constant (shininess)
		specularCol = glm::vec3(term * glm::vec3(1, 1, 1)); // light colour
	}

	if (ray.xindex == 0) {
		return textureBackground(ray);
	}
	if (ray.xindex == 1) {
		return textureFloor(ray);
	}
	if (ray.xindex == 2) {
		return textureLeftside(ray);
	}
	if (ray.xindex == 3) {
		return textureRightside(ray);
	}
	if (ray.xindex == 4) {
		return textureTopside(ray);
	}
	if (ray.xindex == 5) {
		return textureBehindside(ray);
	}
	
	Ray shadow(ray.xpt, unitLightVector);
	shadow.closestPt(sceneObjects);
	float lightDist = glm::length(lightVector);
	   
    if (lDotn <= 0 || (shadow.xindex > -1 && shadow.xdist < lightDist)) { // facing away from view or shadow ray intersects object and object is not behind light source
		return ambientCol * col; // behind the sphere is in shadow (only ambient)
	} else {
		glm::vec3 colorSum(0);
	
		if (ray.xindex == 9) {
			return 0.4f * (ambientCol * col + lDotn * col) + 0.6f * texture2D(ray);
		}
		if (ray.xindex == 6) {
			col = (0.6f * col) + (0.4f * texture2DGreen(ray));
		}

		if (ray.xindex == 8) {
			col = textureProcedural(ray);
			return (ambientCol * col + lDotn * col); // no specular
		}

		if (ray.xindex == 14 || ray.xindex == 15 || ray.xindex == 17) {
			col = textureProceduralBox(ray);
			return (ambientCol * col + lDotn * col); // no specular
		}

		if (ray.xindex == 6 && step < MAX_STEPS) {
			// reflection
			glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
			Ray reflectedRay(ray.xpt, reflectedDir);
			glm::vec3 reflectedCol = trace(reflectedRay, step + 1);
			colorSum = colorSum + (0.9f * reflectedCol);
		}
		
		// refraction
		if (ray.xindex == 7 && step < MAX_STEPS) {
			float eta = 1 / 1.02;
			glm::vec3 refractedDir1 = glm::refract(ray.dir, normalVector, eta);
			
			Ray refractedRay(ray.xpt, refractedDir1);
			refractedRay.closestPt(sceneObjects);
			if (refractedRay.xindex == -1) return backgroundCol;
			glm::vec3 normal = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
			glm::vec3 refractedDir2 = glm::refract(refractedDir1, -normal, 1.0f/eta);

			Ray outputRay(refractedRay.xpt, refractedDir2);
			outputRay.closestPt(sceneObjects);
			if (outputRay.xindex == -1) return backgroundCol;
			glm::vec3 refrCol = trace(outputRay, step+1); // has to be trace
			return (0.4f * lDotn * col) + refrCol;
		}
		// transparent
		if ((ray.xindex == 11) && step < MAX_STEPS) {
			float eta = 1.0;
			glm::vec3 refractedDir1 = glm::refract(ray.dir, normalVector, eta);

			Ray refractedRay(ray.xpt, refractedDir1);
			refractedRay.closestPt(sceneObjects);
			if (refractedRay.xindex == -1) return backgroundCol;
			glm::vec3 normal = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
			glm::vec3 refractedDir2 = glm::refract(refractedDir1, -normal, 1.0f / eta);

			Ray outputRay(refractedRay.xpt, refractedDir2);
			outputRay.closestPt(sceneObjects);
			if (outputRay.xindex == -1) return backgroundCol;
			glm::vec3 refrCol = trace(outputRay, step + 1); // has to be trace
			return ambientCol * col + lDotn * col + refrCol;
		}
		
		return (ambientCol * col + lDotn * col + specularCol + colorSum); // ambient + diffuse + specular + reflection
	}
	
}

//---The main display module -----------------------------------------------------------
// In a ray tracing application, it just displays the ray traced image by drawing
// each cell as a quad.
//---------------------------------------------------------------------------------------
void display()
{
	float xp, yp;  //grid point
	float cellX = (XMAX-XMIN)/NUMDIV;  //cell width
	float cellY = (YMAX-YMIN)/NUMDIV;  //cell height

	glm::vec3 eye(0., 0., 0.);  //The eye position (source of primary rays) is the origin

	glClear(GL_COLOR_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	glBegin(GL_QUADS);  //Each cell is a quad.

	for(int i = 0; i < NUMDIV; i++)  	//For each grid point xp, yp
	{
		xp = XMIN + i*cellX;
		for(int j = 0; j < NUMDIV; j++)
		{
			yp = YMIN + j*cellY;

		    glm::vec3 dir(xp+0.5*cellX, yp+0.5*cellY, -EDIST);	//direction of the primary ray

		    Ray ray = Ray(eye, dir);		//Create a ray originating from the camera in the direction 'dir'
			ray.normalize();				//Normalize the direction of the ray to a unit vector
		    glm::vec3 col = trace (ray, 1); //Trace the primary ray and get the colour value

			glColor3f(col.r, col.g, col.b);
			glVertex2f(xp, yp);				//Draw each cell with its color value
			glVertex2f(xp+cellX, yp);
			glVertex2f(xp+cellX, yp+cellY);
			glVertex2f(xp, yp+cellY);
        }
    }

    glEnd();
    glFlush();
}




//---This function initializes the scene ------------------------------------------- 
//   Specifically, it creates scene objects (spheres, planes, cones, cylinders etc)
//     and add them to the list of scene objects.
//   It also initializes the OpenGL orthographc projection matrix for drawing the
//     the ray traced image.
//----------------------------------------------------------------------------------
void initialize()
{
    glMatrixMode(GL_PROJECTION);
    gluOrtho2D(XMIN, XMAX, YMIN, YMAX);
    glClearColor(0, 0, 0, 1);

	texturePurple = TextureBMP("pink-planet.bmp");
	textureBlue = TextureBMP("plain-blue-space.bmp");
	textureGreen = TextureBMP("green-planet.bmp");

	Plane *background = new Plane(glm::vec3(-40, -20, -200), // back left
		glm::vec3(40, -20, -200), // back right
		glm::vec3(40, 60, -200),  // top right
		glm::vec3(-40, 60, -200), // top left
		glm::vec3(0.5, 0.5, 0)); // default colour

	Plane *floor = new Plane(glm::vec3(-40, -20, -40), // a - front left
		glm::vec3(40, -20, -40), // b - front right
		glm::vec3(40, -20, -200), // c - back right
		glm::vec3(-40, -20, -200), // d - back left
		glm::vec3(0.5, 0, 0.5)); // colour

	Plane *leftside = new Plane(glm::vec3(-40, -20, -40), // front bottom
		glm::vec3(-40, -20, -200), // back bottom
		glm::vec3(-40, 60, -200),  // top back
		glm::vec3(-40, 60, -40), // top front
		glm::vec3(0.5, 0, 0)); // default colour

	Plane *rightside = new Plane(glm::vec3(40, -20, -40), // front bottom
		glm::vec3(40, 60, -40), // top front
		glm::vec3(40, 60, -200),  // top back
		glm::vec3(40, -20, -200), // back bottom
		glm::vec3(0, 0.5, 0)); // default colour
	
	Plane *topside = new Plane(glm::vec3(-40, 60, -40), // front bottom
		glm::vec3(-40, 60, -200), // top front
		glm::vec3(40, 60, -200),  // top back
		glm::vec3(40, 60, -40), // back bottom
		glm::vec3(1, 1, 1)); // default colour

	Plane *behindside = new Plane(glm::vec3(-80, -40, 5),
		glm::vec3(-80, 80, 5),
		glm::vec3(80, 80, 5),
		glm::vec3(80, -40, 5),
		glm::vec3(1, 1, 1)); // default colour	

	//-- Create a pointer to a sphere object: x, y, z, radius, color
	Sphere *sphereBlue = new Sphere(glm::vec3(-5.0, -5.0, -145.0), 15.0, glm::vec3(0, 0, 0.8));
	Sphere *sphereRed = new Sphere(glm::vec3(5.0, 2, -125.0), 2.3, glm::vec3(0.5, 0.5, 0.5));
	Cylinder *cylinderRed = new Cylinder(glm::vec3(17, -10, -100), 4, 10, glm::vec3(1, 0, 1));
	Sphere *sphereGreen = new Sphere(glm::vec3(15.0, 5, -185.0), 20.0, glm::vec3(1, 0, 0));
	Sphere *sphereGrey = new Sphere(glm::vec3(-18, -10, -115.0), 8.0, glm::vec3(0.6, 0, 0.6));
	Cone *coneGrey = new Cone(glm::vec3(17, 0, -100), 4, 4, glm::vec3(0.9, 0, 0.9));
	Sphere *sphereOrange = new Sphere(glm::vec3(10, -15, -85.0), 4.0, glm::vec3(1, 0.6, 0));

	// box
	float min_x = 0; // -20;
	float max_x = 10; //-10;
	float min_y = 0; //10;
	float max_y = 5; //15;
	float min_z = 0; //-110;
	float max_z = -30; //-140;
	glm::vec3 boxcol(0.5, 0.5, 0.2);

	
	// translate
	glm::mat4 t(1, 0, 0, -20,
				0, 1, 0, 10,
				0, 0, 1, -110,
				0, 0, 0, 1);
	
	// shear
	glm::mat4 s(1, -0.7, -0.1, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1);

	glm::vec3 flb(min_x, min_y, min_z);
	glm::vec3 flt(min_x, max_y, min_z);
	glm::vec3 frb(max_x, min_y, min_z);
	glm::vec3 frt(max_x, max_y, min_z);
	glm::vec3 blt(min_x, max_y, max_z);
	glm::vec3 blb(min_x, min_y, max_z);
	glm::vec3 brt(max_x, max_y, max_z);
	glm::vec3 brb(max_x, min_y, max_z);
	
	glm::vec4 flb4 = glm::vec4(glm::vec3(flb), 1.0) * s * t;
	glm::vec4 flt4 = glm::vec4(glm::vec3(flt), 1.0) * s * t;
	glm::vec4 frb4 = glm::vec4(glm::vec3(frb), 1.0) * s * t;
	glm::vec4 frt4 = glm::vec4(glm::vec3(frt), 1.0) * s * t;
	glm::vec4 blt4 = glm::vec4(glm::vec3(blt), 1.0) * s * t;
	glm::vec4 blb4 = glm::vec4(glm::vec3(blb), 1.0) * s * t;
	glm::vec4 brt4 = glm::vec4(glm::vec3(brt), 1.0) * s * t;
	glm::vec4 brb4 = glm::vec4(glm::vec3(brb), 1.0) * s * t;
	
	flb = glm::vec3(flb4.x / flb4.w, flb4.y / flb4.w, flb4.z / flb4.w); // normalise
	flt = glm::vec3(flt4.x / flt4.w, flt4.y / flt4.w, flt4.z / flt4.w);
	frb = glm::vec3(frb4.x / frb4.w, frb4.y / frb4.w, frb4.z / frb4.w);
	frt = glm::vec3(frt4.x / frt4.w, frt4.y / frt4.w, frt4.z / frt4.w);
	blb = glm::vec3(blb4.x / blb4.w, blb4.y / blb4.w, blb4.z / blb4.w);
	blt = glm::vec3(blt4.x / blt4.w, blt4.y / blt4.w, blt4.z / blt4.w);
	brb = glm::vec3(brb4.x / brb4.w, brb4.y / brb4.w, brb4.z / brb4.w);
	brt = glm::vec3(brt4.x / brt4.w, brt4.y / brt4.w, brt4.z / brt4.w);
	
	Plane *left = new Plane(blb, flb, flt, blt, boxcol); 
	Plane *right = new Plane(brt, frt, frb, brb, boxcol); 
	Plane *bottom = new Plane(flb, blb, brb, frb, boxcol);
	Plane *top = new Plane(flt, frt, brt, blt, boxcol);
	Plane *front = new Plane(flb, frb, frt, flt, boxcol); 
	Plane *back = new Plane(blb, blt, brt, brb, boxcol);
	
	//--Add the above to the list of scene objects.
	sceneObjects.push_back(background); //0 star texture
	sceneObjects.push_back(floor); // star texture
	sceneObjects.push_back(leftside); // star texture
	sceneObjects.push_back(rightside); // star texture
	sceneObjects.push_back(topside); // star texture
	sceneObjects.push_back(behindside); // star texture

	sceneObjects.push_back(sphereBlue); //6 reflective sphere
	sceneObjects.push_back(cylinderRed); //7 refractive cylinder - near procedural
	sceneObjects.push_back(sphereGreen); //8 procedural sphere
	sceneObjects.push_back(sphereGrey); //9 2d texture sphere - near reflective
	sceneObjects.push_back(coneGrey); //10 cone - on top of cylinder
	sceneObjects.push_back(sphereOrange); //11 transparent sphere
	sceneObjects.push_back(sphereRed); //12 shadow caster
	
	sceneObjects.push_back(left); //13 sheared cube - shadow caster, near reflective
	sceneObjects.push_back(right); // sheared cube
	sceneObjects.push_back(bottom); // sheared cube
	sceneObjects.push_back(top); // sheared cube
	sceneObjects.push_back(front); // sheared cube
	sceneObjects.push_back(back); // sheared cube
	
}



int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB );
    glutInitWindowSize(600, 600);
    glutInitWindowPosition(20, 20);
    glutCreateWindow("Raytracer");

    glutDisplayFunc(display);
    initialize();

    glutMainLoop();
    return 0;
}
