/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The object class
*  This is a generic type for storing objects in the scene.
*  Being an abstract class, this class cannot be instantiated.
*  Sphere, Plane etc, must be defined as subclasses of Object
*      and provide implementations for the virtual functions
*      intersect()  and normal().
-------------------------------------------------------------*/

#ifndef H_SOBJECT
#define H_SOBJECT
#include <glm/glm.hpp>


class SceneObject 
{
protected:
	glm::vec3 color;
	glm::mat4 transform;
	glm::mat4 transformInv;
	glm::mat4 transformInvTrans;

	virtual float intersect_internal(glm::vec3 pos, glm::vec3 dir) = 0;
	virtual glm::vec3 normal_internal(glm::vec3 pos) = 0;
public:
	SceneObject() {}
	float intersect(glm::vec3 pos, glm::vec3 dir);
	glm::vec3 normal(glm::vec3 pos);

	virtual ~SceneObject() {}
	glm::vec3 getColor();
	void setColor(glm::vec3 col);
	glm::mat4 getTransform();
	void setTransform(glm::mat4 col);
};

#endif
