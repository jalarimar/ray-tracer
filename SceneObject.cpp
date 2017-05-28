/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The object class
*  This is a generic type for storing objects in the scene
*  Sphere, Plane etc. must be defined as subclasses of Object.
*  Being an abstract class, this class cannot be instantiated.
-------------------------------------------------------------*/

#include "SceneObject.h"

glm::vec3 SceneObject::getColor()
{
	return color;
}

void SceneObject::setColor(glm::vec3 col)
{
	color = col;
}

glm::mat4 SceneObject::getTransform()
{
	return transform;
}

void SceneObject::setTransform(glm::mat4 trans)
{
	transformInv = glm::inverse(trans);
	transformInvTrans = glm::transpose(glm::inverse(trans));
	transform = trans;
}

glm::vec3 transformProj(glm::mat4 matx, glm::vec3 pos) {
	glm::vec4 p4 = matx * glm::vec4(pos, 1.0);
	return glm::vec3(p4.x / p4.w); // TODO
}

float SceneObject::intersect(glm::vec3 pos, glm::vec3 dir) {
	return intersect_internal(transformProj(transformInv, pos), dir);
}

glm::vec3 SceneObject::normal(glm::vec3 pos) {
	return normal_internal(pos);
}