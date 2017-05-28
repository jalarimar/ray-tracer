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

vector<SceneObject*> sceneObjects;  //A global list containing pointers to objects in the scene
TextureBMP texturePurple;
TextureBMP textureBlue;


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
    
	/** Remember to change indices when uncomment
	// background 
	if (ray.xindex == 0) { // INDEC
		float texcoords = (ray.xpt.x + 40) / 80;
		float texcoordt = (ray.xpt.y + 20) / 80;
		return textureBlue.getColorAt(texcoords, texcoordt);
	}
	// floor
	if (ray.xindex == 1) {
		float texcoords = (ray.xpt.x + 40) / 100;
		float texcoordt = (ray.xpt.z + 200 + 40) / 200;
		return textureBlue.getColorAt(texcoords, 1 - texcoordt);
	}
	// leftside
	if (ray.xindex == 2) { // INDEC
		float texcoords = (ray.xpt.z + 200) / 200;
		float texcoordt = (ray.xpt.y + 20 + 10) / 140;
		return textureBlue.getColorAt(texcoords, 1 - texcoordt);
	}
	// rightside
	if (ray.xindex == 3) { // INDEC
		float texcoords = (ray.xpt.z + 200) / 200;
		float texcoordt = (ray.xpt.y + 20) / 80;
		return textureBlue.getColorAt(1 - texcoords, texcoordt);
	}
	
	// topside
	if (ray.xindex == 4) { // INDEC
		float texcoords = (ray.xpt.x + 40) / 100;
		float texcoordt = (ray.xpt.z + 240) / 200;
		return textureBlue.getColorAt(texcoords, 1 - texcoordt);
	}
	// behindside
	if (ray.xindex == 5) { // INDEC
		float texcoords = (ray.xpt.x + 80) / 160;
		float texcoordt = (ray.xpt.y + 40) / 160;
		return textureBlue.getColorAt(texcoords, texcoordt);
	}
	*/

	// sphere texture
	//if (ray.xindex == 3) { //INDEC
	//	float pi = 3.14;
	//	glm::vec3 n = sceneObjects[ray.xindex]->normal(ray.xpt);
	//	float norx = (ray.xpt.x + 15) / 8;
	//	float nory = (ray.xpt.y + 10) / 8;
	//	float norz = (ray.xpt.z + 115) / 8;
	//	float texcoords = 0.5 + atan2(norz, norx) / (2 * pi);
	//	float texcoordt = 0.5 - asin(nory) / pi;
	//	return texturePurple.getColorAt(1 - (texcoords/2+texcoords/4), (texcoordt/2+texcoordt/4));
	//	
	//	/*float radius = 16;
	//	float t = acos(norz / 4) / pi;
	//	if (nory >= 0) {
	//		float s = acos(norx / (radius * sin(pi*t))) / (2*pi);
	//		return texturePurple.getColorAt(s, t);
	//	} else {
	//		float s = (pi + acos(norx / (radius * sin(pi*t)))) / (2*pi);
	//		return texturePurple.getColorAt(s, t);
	//	}*/
	//}
	/*
	// procedural pattern
	if (ray.xindex == 1) { // INDEC
		float norx = (ray.xpt.x - 2.5) / 16; // red 4.6
		float nory = (ray.xpt.y - 5) / 20; // green
		PerlinNoise perlin;
		double noise = perlin.noise(nory * 5, norx * 5, 0);
		//col = glm::vec3(1, ((1 + sin(nory + noise)) / 2), 0);
		//double random = ((double)rand() / (double)RAND_MAX);
		col = glm::vec3(1, (1 + sin(nory + 20 * noise)) / 2.5, 0); // decreasing 20 increases smokiness, decreasing 2.5 decreases blending
	}
	*/
	
    //--- mine ---
    glm::vec3 normalVector = sceneObjects[ray.xindex] -> normal(ray.xpt);
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
	
	Ray shadow(ray.xpt, unitLightVector);
	shadow.closestPt(sceneObjects);
	float lightDist = glm::length(lightVector);
	   
    if (lDotn <= 0 || (shadow.xindex > -1 && shadow.xdist < lightDist)) { // facing away from view or shadow ray intersects object and object is not behind light source
		return ambientCol * col; // behind the sphere is in shadow (only ambient)
	} else {
		
		glm::vec3 colorSum(0);
		// reflective
		/*if (ray.xindex == 3 && step < MAX_STEPS) {
			glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
			Ray reflectedRay(ray.xpt, reflectedDir);
			glm::vec3 reflectedCol = trace(reflectedRay, step + 1);
			colorSum = colorSum + (0.9f * reflectedCol);
		}*/
		/*
		if (ray.xindex == 2) { // procedural
			// no specular, no reflection INDEC
			return (ambientCol * col + lDotn * col);
		}*/
		
		// refraction
		if (ray.xindex == 0 && step < MAX_STEPS) {
			float eta = 1/1.05;
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
			return refrCol;
		}
		// transparent
		/*if (ray.xindex == 3 && step < MAX_STEPS) {
			float eta = 1.0;
			glm::vec3 refractedDir1 = glm::refract(ray.dir, normalVector, eta);

			Ray refractedRay(ray.xpt, refractedDir1);
			refractedRay.closestPt(sceneObjects);
			glm::vec3 normal = sceneObjects[refractedRay.xindex]->normal(refractedRay.xpt);
			glm::vec3 refractedDir2 = glm::refract(refractedDir1, -normal, 1.0f / eta);

			Ray outputRay(refractedRay.xpt, refractedDir2);
			outputRay.closestPt(sceneObjects);
			if (outputRay.xindex == -1) return lDotn * col;
			glm::vec3 refrCol = trace(outputRay, step + 1);
			return ambientCol * col + lDotn * col + refrCol;
		}*/
		
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

	// TODO BEFORE SUBMISSION find a way of loading based on relative path
	texturePurple = TextureBMP("..\\purple-space.bmp");
	textureBlue = TextureBMP("C:\\Users\\Jay\\Documents\\Engineering2017\\363\\Assignment2\\ray-tracer\\plain-blue-space.bmp");

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
	Sphere *sphereBlue = new Sphere(glm::vec3(-5.0, -5.0, -150.0), 15.0, glm::vec3(0, 0, 0.8));
	//Sphere *sphereRed = new Sphere(glm::vec3(5.0, 2, -130.0), 2.3, glm::vec3(1, 0, 0));
	Cylinder *cylinderRed = new Cylinder(glm::vec3(15, -30, -100), 4, 15, glm::vec3(1, 0, 0));
	Sphere *sphereGreen = new Sphere(glm::vec3(15.0, -16, -185.0), 25.0, glm::vec3(0, 1, 0));
	Sphere *sphereGrey = new Sphere(glm::vec3(-15, -10, -115.0), 8.0, glm::vec3(0.7, 0.7, 0.7));
	Cone *coneGrey = new Cone(glm::vec3(15, -15, -100), 4, 5, glm::vec3(0.7, 0.7, 0.7));
	Sphere *sphereOrange = new Sphere(glm::vec3(15, -10, -115.0), 8.0, glm::vec3(0.7, 0.7, 0.7));

	
	// box
	float min_x = -10;//20;
	float max_x = 30; // 20 is out of range??
	float min_y = -15;
	float max_y = -12; // -10 is out of range??
	float min_z = -110;
	float max_z = -120; // wtf
	glm::vec3 boxcol(1, 0.7, 0);

	/*
	// translate using mat4 before shearing?
	glm::mat4 t(1, 0, 0, min_x,
				0, 1, 0, min_y,
				0, 0, 1, min_z,
				0, 0, 0, 1);

	glm::mat4 m(1, 2, 0, 0,
				0, 1, 0, 0,
				0, 0, 1, 0,
				0, 0, 0, 1);
	*/
	glm::vec3 flb(min_x, min_y, min_z);
	glm::vec3 flt(min_x, max_y, min_z);
	glm::vec3 frb(max_x, min_y, min_z);
	glm::vec3 frt(max_x, max_y, min_z);
	glm::vec3 blt(min_x, max_y, max_z);
	glm::vec3 blb(min_x, min_y, max_z);
	glm::vec3 brt(max_x, max_y, max_z);
	glm::vec3 brb(max_x, min_y, max_z);
	/*
	glm::vec4 flb4 = m * glm::vec4(glm::vec3(flb), 1.0);
	glm::vec4 flt4 = m * glm::vec4(glm::vec3(flt), 1.0);
	glm::vec4 frb4 = m * glm::vec4(glm::vec3(frb), 1.0);
	glm::vec4 frt4 = m * glm::vec4(glm::vec3(frt), 1.0);
	glm::vec4 blb4 = m * glm::vec4(glm::vec3(blt), 1.0);
	glm::vec4 blt4 = m * glm::vec4(glm::vec3(blb), 1.0);
	glm::vec4 brb4 = m * glm::vec4(glm::vec3(brt), 1.0);
	glm::vec4 brt4 = m * glm::vec4(glm::vec3(brb), 1.0);
	
	flb = glm::vec3(flb4.x / flb4.w, flb4.y / flb4.w, flb4.z / flb4.w); // normalise
	flt = glm::vec3(flt4.x / flt4.w, flt4.y / flt4.w, flt4.z / flt4.w);
	frb = glm::vec3(frb4.x / frb4.w, frb4.y / frb4.w, frb4.z / frb4.w);
	frt = glm::vec3(frt4.x / frt4.w, frt4.y / frt4.w, frt4.z / frt4.w);
	blb = glm::vec3(blb4.x / blb4.w, blb4.y / blb4.w, blb4.z / blb4.w);
	blt = glm::vec3(blt4.x / blt4.w, blt4.y / blt4.w, blt4.z / blt4.w);
	brb = glm::vec3(brb4.x / brb4.w, brb4.y / brb4.w, brb4.z / brb4.w);
	brt = glm::vec3(brt4.x / brt4.w, brt4.y / brt4.w, brt4.z / brt4.w);
	*/
	//flb = glm::normalize(glm::vec4(glm::vec3(flb), 1.0));
	/*
	flb = shear * flb;
	flt = shear * flt;
	frb = shear * frb;
	frt = shear * frt;
	blb = shear * blb;
	blt = shear * blt;
	brb = shear * brb;
	brt = shear * brt;
	*/
	
	Plane *left = new Plane(blb, flb, flt, blt, boxcol); 
	Plane *right = new Plane(brt, frt, frb, brb, boxcol); 
	Plane *bottom = new Plane(flb, blb, brb, frb, boxcol);
	Plane *top = new Plane(flt, frt, brt, blt, boxcol);
	Plane *front = new Plane(flb, frb, frt, flt, boxcol); 
	Plane *back = new Plane(blb, blt, brt, brb, boxcol);
	
	//--Add the above to the list of scene objects.
	//sceneObjects.push_back(background); // star texture
	//sceneObjects.push_back(floor); // star texture
	//sceneObjects.push_back(leftside); // star texture
	//sceneObjects.push_back(rightside); // star texture
	//sceneObjects.push_back(topside); // star texture
	//sceneObjects.push_back(behindside); // star texture

	sceneObjects.push_back(sphereBlue); // reflective sphere
	sceneObjects.push_back(cylinderRed); // refractive cylinder - near procedural
	sceneObjects.push_back(sphereGreen); // procedural sphere
	sceneObjects.push_back(sphereGrey); // 2d texture sphere - near reflective
	sceneObjects.push_back(coneGrey); // cone - on top of cylinder
	sceneObjects.push_back(sphereOrange); // transparent sphere - near reflective
	
	//sceneObjects.push_back(left); // sheared cube - shadow caster, near reflective
	//sceneObjects.push_back(right); // sheared cube
	//sceneObjects.push_back(bottom); // sheared cube
	//sceneObjects.push_back(top); // sheared cube
	//sceneObjects.push_back(front); // sheared cube
	//sceneObjects.push_back(back); // sheared cube
	
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
