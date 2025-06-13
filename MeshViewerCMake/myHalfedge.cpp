#include "myHalfedge.h"
#include "myVertex.h"
#include <iostream>
#include <ostream>
#include <cmath>

myHalfedge::myHalfedge(void)
{
	source = NULL; 
	adjacent_face = NULL; 
	next = NULL;  
	prev = NULL;  
	twin = NULL; 
	double length = -1.0;
}

void myHalfedge::copy(myHalfedge *ie)
{
/**** TODO ****/
}

myHalfedge::~myHalfedge(void)
{
	
}

void myHalfedge::calculateLength() {
	if (!source || !source->point || !twin || !twin->source || !twin->source->point) {
        std::cout << "Error: Null pointer in myHalfedge::calculateLength()" << std::endl;
        length = 0.0;
		return;
    }
	double vX = source->point->X - twin->source->point->X;
	double vY = source->point->Y - twin->source->point->Y;
	double vZ = source->point->Z - twin->source->point->Z;

	this->length = sqrt(vX*vX + vY*vY + vZ*vZ);
}
