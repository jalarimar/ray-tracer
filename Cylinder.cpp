/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cylinder class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cylinder.h"
#include "Ray.h"
#include <math.h>

const double PI = 3.14159;

/**
* Cylinder's constructor
*/
Cylinder::Cylinder(glm::vec3 centre, float rad, float hight, glm::vec3 col)
{
	center = centre;
	radius = rad;
	height = hight;
	color = col;
}

/**
* Cylinder's intersection method.  The input is a ray (pos, dir). 
*/
float Cylinder::intersect_internal(glm::vec3 posn, glm::vec3 dir)
{
	float xc = center.x;
	float yc = center.y;
	float zc = center.z;

	float x0 = posn.x;
	float y0 = posn.y;
	float z0 = posn.z;

	float dx = dir.x;
	float dy = dir.y;
	float dz = dir.z;

	// ray-cylinder intersect: t^2(dx^2+dz^2) + 2t(dx(x0-xc) + dz(z0-zc)) + (x0-xc)^2 + (z0-zc)^2 - r^2 = 0
	// quadratic formula: ax^2 + bx + c = 0 ==> x = (-b +- sqrt(b^2 - 4ac)) / 2a
	float a = dx*dx + dz*dz;
	float b = 2 * (dx*(x0 - xc) + dz*(z0 - zc));
	float c = (x0 - xc)*(x0 - xc) + (z0 - zc)*(z0 - zc) - (radius * radius);
	float discrim = b*b - 4*a*c;

	if (fabs(discrim) < 0.001) return -1.0;
	if (discrim < 0.0) return -1.0;

	float t1 = (-b - sqrt(discrim)) / (2*a);
	float t2 = -b + sqrt(discrim) / (2*a);
	float t = -5.0;
	float tf = -5.0;

	if (fabs(t2) < 0.001) t2 = -1.0;
	if (fabs(t1) < 0.001)
	{
		if (t2 > 0) t = t2;
		else t1 = -1.0;
	}

	if (t == -5.0) {
		t = (t1 < t2) ? t1 : t2;
		tf = (t1 < t2) ? t2 : t1;
	}
	// return t;
	if (t < 0.0) {
		return t;
	}

	float y = y0 + dy*t;
	if ((y - yc) < 0.0 || (y - yc) > height) {
		float y2 = y0 + dy*tf;
		if (!((y2 - yc) < 0.0 || (y2 - yc) > height) && fabs(tf) < 0.001) {
			return tf; // valid second point
		}
		return -1.0; //invalid
	}
	return t;
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the cylinder.
*/
glm::vec3 Cylinder::normal_internal(glm::vec3 point)
{
	float x = point.x;
	float xc = center.x;
	float z = point.z;
	float zc = center.z;
    glm::vec3 n((x-xc)/radius, 0, (z-zc)/radius);
    return n;
}
