#include "myHalfedge.h"
#include "myVertex.h"

#include <cmath>

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL;  
}

void myHalfedge::copy(myHalfedge *ie)
{
/**** TODO ****/
}

myHalfedge::~myHalfedge(void)
{
	
}

double myHalfedge::length() {
	double vX = source->point->X - twin->source->point->X;
	double vY = source->point->Y - twin->source->point->Y;
	double vZ = source->point->Z - twin->source->point->Z;
	 
	return sqrt(vX*vX + vY*vY + vZ*vZ);
}
