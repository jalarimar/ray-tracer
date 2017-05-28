/*----------------------------------------------------------
* COSC363  Ray Tracer
*
*  The cone class
*  This is a subclass of Object, and hence implements the
*  methods intersect() and normal().
-------------------------------------------------------------*/

#ifndef H_CONE
#define H_CONE
#include <glm/glm.hpp>
#include "SceneObject.h"

/**
 * Defines a simple Cone located at 'center' 
 * with the specified radius
 */
class Cone : public SceneObject
{

private:
    glm::vec3 center;
    float radius;
	float height;

public:	
	Cone(glm::vec3 center, float radius, float height, glm::vec3 col);

	float intersect_internal(glm::vec3 posn, glm::vec3 dir);

	glm::vec3 normal_internal(glm::vec3 p);

};

#endif //!H_CONE
