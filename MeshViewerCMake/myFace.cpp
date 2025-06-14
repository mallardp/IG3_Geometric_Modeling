#include "myFace.h"
#include "myVector3D.h"
#include "myHalfedge.h"
#include "myVertex.h"
#include <GL/glew.h>

myFace::myFace(void)
{
	adjacent_halfedge = NULL;
	normal = new myVector3D(1.0, 1.0, 1.0);
}

myFace::~myFace(void)
{
	if (normal) delete normal;
}

int myFace::countEdges(){
	myHalfedge *h = adjacent_halfedge;
	int count = 0;
	if (h == NULL) return 0;

	do {
		count++;
		h = h->next;
	} while (h != adjacent_halfedge);

	return count;
}

void myFace::computeNormal()
{
	myHalfedge *h = adjacent_halfedge;
	myVector3D v1, v2;
	v1.dX = h->source->point->X - h->next->source->point->X;
	v1.dY = h->source->point->Y - h->next->source->point->Y;
	v1.dZ = h->source->point->Z - h->next->source->point->Z;

	h = h->next;
	v2.dX = h->source->point->X - h->next->source->point->X;
	v2.dY = h->source->point->Y - h->next->source->point->Y;
	v2.dZ = h->source->point->Z - h->next->source->point->Z;

	normal->crossproduct(v1, v2);
	normal->normalize();
}
