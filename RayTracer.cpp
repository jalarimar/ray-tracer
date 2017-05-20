/*========================================================================
* COSC 363  Computer Graphics (2017)
* Ray tracer
*=========================================================================
*/
#include <iostream>
#include <cmath>
#include <vector>
#include <glm/glm.hpp>
#include "Sphere.h"
#include "Plane.h"
#include "SceneObject.h"
#include "Ray.h"
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
TextureBMP texture;

//---The most important function in a ray tracer! ---------------------------------- 
//   Computes the colour value obtained by tracing a ray and finding its 
//     closest point of intersection with objects in the scene.
//----------------------------------------------------------------------------------
glm::vec3 trace(Ray ray, int step)
{
	glm::vec3 backgroundCol(0);
	glm::vec3 ambientCol(0.1);
	glm::vec3 light(10, 40, -3);

    ray.closestPt(sceneObjects);		//Compute the closest point of intersetion of objects with the ray

    if(ray.xindex == -1) return backgroundCol;      //If there is no intersection return background colour

    glm::vec3 col = sceneObjects[ray.xindex]->getColor(); // object's colour
    
	// background
	if (ray.xindex == 0) {
		float texcoords = (ray.xpt.x + 40) / 80;
		float texcoordt = (ray.xpt.y + 20) / 80;
		return texture.getColorAt(texcoords, texcoordt);
	}
	// floor
	if (ray.xindex == 1) {
		float texcoords = (ray.xpt.x + 40) / 100;
		float texcoordt = (ray.xpt.z + 200 + 40) / 200;
		return texture.getColorAt(texcoords, 1 - texcoordt);
	}
	// leftside
	if (ray.xindex == 2) {
		float texcoords = (ray.xpt.z + 200) / 200;
		float texcoordt = (ray.xpt.y + 20 + 10) / 140;
		return texture.getColorAt(texcoords, 1 - texcoordt);
	}
	// rightside
	if (ray.xindex == 3) {
		float texcoords = (ray.xpt.z + 200) / 200;
		float texcoordt = (ray.xpt.y + 20) / 80;
		return texture.getColorAt(1 - texcoords, texcoordt);
	}
	
	// topside
	if (ray.xindex == 4) {
		float texcoords = (ray.xpt.x + 40) / 100;
		float texcoordt = (ray.xpt.z + 240) / 200;
		return texture.getColorAt(texcoords, 1 - texcoordt);
	}
	// behindside
	if (ray.xindex == 5) {
		float texcoords = (ray.xpt.x + 80) / 160;
		float texcoordt = (ray.xpt.y + 40) / 160;
		return texture.getColorAt(texcoords, texcoordt);
	}

    /// mine
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
	
	// Lab 8
	Ray shadow(ray.xpt, unitLightVector);
	shadow.closestPt(sceneObjects);
	float lightDist = glm::length(lightVector);
	   
    if (lDotn <= 0 || (shadow.xindex > -1 && shadow.xdist < lightDist)) { // facing away from view or shadow ray intersects object and object is not behind light source
		return ambientCol * col; // behind the sphere is in shadow (only ambient)
	} else {
		
		// Lab 8
		glm::vec3 colorSum(0);
		// makes the blue sphere reflective
		if (ray.xindex == 6 && step < MAX_STEPS) {
			glm::vec3 reflectedDir = glm::reflect(ray.dir, normalVector);
			Ray reflectedRay(ray.xpt, reflectedDir);
			glm::vec3 reflectedCol = trace(reflectedRay, step + 1);
			colorSum = colorSum + (0.9f * reflectedCol);
		}
		
		return (ambientCol * col + lDotn * col + specularCol + colorSum); // ambient + diffuse + specular
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
	texture = TextureBMP("C:\\Users\\Jay\\Documents\\Engineering2017\\363\\Assignment2\\ray-tracer\\plain-blue-space.bmp");

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
	Sphere *sphereRed = new Sphere(glm::vec3(5.0, 2, -130.0), 2.3, glm::vec3(1, 0, 0));
	Sphere *sphereGreen = new Sphere(glm::vec3(15.0, 10, -85.0), 5.0, glm::vec3(0, 1, 0));
	Sphere *sphereGrey = new Sphere(glm::vec3(-10, -10, -65.0), 4.0, glm::vec3(0.7, 0.7, 0.7));
	
	// box
	float min_x = 20;
	float max_x = 30;
	float min_y = -20;
	float max_y = -17;
	float min_z = -110;
	float max_z = -140;
	glm::vec3 boxcol(0.7, 0.7, 0.7);
	
	glm::vec3 flb(min_x, min_y, min_z);
	glm::vec3 flt(min_x, max_y, min_z);
	glm::vec3 frb(max_x, min_y, min_z);
	glm::vec3 frt(max_x, max_y, min_z);
	glm::vec3 blt(min_x, max_y, max_z);
	glm::vec3 blb(min_x, min_y, max_z);
	glm::vec3 brt(max_x, max_y, max_z);
	glm::vec3 brb(max_x, min_y, max_z);
	
	Plane *left = new Plane(blb, flb, flt, blt, boxcol); 
	Plane *right = new Plane(brt, frt, frb, brb, boxcol); 
	Plane *bottom = new Plane(flb, blb, brb, frb, boxcol);
	Plane *top = new Plane(flt, frt, brt, blt, boxcol);
	Plane *front = new Plane(flb, frb, frt, flt, boxcol); 
	Plane *back = new Plane(blb, blt, brt, brb, boxcol);
	

	//--Add the above to the list of scene objects.
	sceneObjects.push_back(background); // 0
	sceneObjects.push_back(floor);
	sceneObjects.push_back(leftside);
	sceneObjects.push_back(rightside);
	sceneObjects.push_back(topside);
	sceneObjects.push_back(behindside);

	sceneObjects.push_back(sphereBlue); // 6
	sceneObjects.push_back(sphereRed);
	//sceneObjects.push_back(sphereGreen);
	//sceneObjects.push_back(sphereGrey);
	
	sceneObjects.push_back(left); // 10
	sceneObjects.push_back(right);
	sceneObjects.push_back(bottom);
	sceneObjects.push_back(top);
	sceneObjects.push_back(front);
	sceneObjects.push_back(back);
	
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
