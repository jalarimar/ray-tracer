/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cone class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#include "Cone.h"
#include "Ray.h"
#include <math.h>

const double PI = 3.14159;

/**
* Cone's constructor
*/
Cone::Cone(glm::vec3 centre, float rad, float hight, glm::vec3 col)
{
	center = centre;
	radius = rad;
	height = hight;
	color = col;
}

/**
* Cone's intersection method.  The input is a ray (pos, dir). 
*/
float Cone::intersect(glm::vec3 posn, glm::vec3 dir)
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

	// ray-cone
	// (x-xc)2 + (z-zc)2 = (R/h)2 * (h-y+yc)2
	// (x-xc)2 + (z-zc)2 - (R/h)2 * (h-y+yc)2 = 0
	// ((x0+dx*t)-xc)2 + ((z0+dz*t)-zc)2 - (R/h)2 * (h-(y0+dy*t)+yc)2 = 0
	// (x0+dx*t)2-2(xc(x0+dx*t))+xc2 + (z0+dz*t)2-2(zc(z0+dz*t))+zc2 - (R/h)2 * (h-y0-dy*t+yc)2 = 0
	// quadratic formula: ax^2 + bx + c = 0 ==> x = (-b +- sqrt(b^2 - 4ac)) / 2a
	float a = dx*dx + dz*dz - (radius/height)*(radius/height)*dy*dy;
	float b = 2 * (dx*(x0 - xc) + dz*(z0 - zc) + (radius/height)*(radius/height)*dy*(height-y0+yc));
	float c = (x0 - xc)*(x0 - xc) + (z0 - zc)*(z0 - zc) - (radius/height)*(radius/height)*(height-y0+yc);
	float discrim = b*b - 4*a*c;

	if (fabs(discrim) < 0.001) return -1.0;
	if (discrim < 0.0) return -1.0;

	float t1 = (-b - sqrt(discrim)) / (2*a);
	float t2 = -b + sqrt(discrim) / (2*a);

	if (fabs(t2) < 0.001) t2 = -1.0;
	if (fabs(t1) < 0.001)
	{
		if (t2 > 0) return t2;
		else t1 = -1.0;
	}
	return (t1 < t2) ? t1 : t2; // the smallest
}

/**
* Returns the unit normal vector at a given point.
* Assumption: The input point p lies on the cone.
*/
glm::vec3 Cone::normal(glm::vec3 point)
{
	float x = point.x;
	float xc = center.x;
	float z = point.z;
	float zc = center.z;

	float alpha = atan(x-xc / z-zc);
	float theta = atan(radius / height);

    glm::vec3 n = glm::vec3(sin(alpha) * cos(theta), sin(theta), cos(alpha) * cos(theta));
    return n;
}
