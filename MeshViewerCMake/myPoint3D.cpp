#include "myPoint3D.h"
#include "myVector3D.h"
#include <iostream>
#include <cmath>

myPoint3D::myPoint3D()
{
	X = Y = Z = 0.0;
}

myPoint3D::myPoint3D(double x, double y, double z)
{
	X = x;
	Y = y;
	Z = z;
}

myPoint3D myPoint3D::operator+(myVector3D  v1)
{
	return myPoint3D(X + v1.dX, Y + v1.dY, Z + v1.dZ);
}

myPoint3D myPoint3D::operator+(myPoint3D  v1)
{
	return myPoint3D(X + v1.X, Y + v1.Y, Z + v1.Z);
}
myPoint3D & myPoint3D::operator+=(myVector3D  v1)
{
	X += v1.dX;
	Y += v1.dY;
	Z += v1.dZ;
	return *this;
}

myPoint3D & myPoint3D::operator+=(myPoint3D  v1)
{
	X += v1.X;
	Y += v1.Y;
	Z += v1.Z;
	return *this;
}

myPoint3D & myPoint3D::operator/=(double d)
{
	if (d != 0.0)
	{
		X /= d;
		Y /= d;
		Z /= d;
	}
	return *this;
}

myPoint3D & myPoint3D::operator*=(double d)
{
	X *= d;
	Y *= d;
	Z *= d;
	return *this;
}

myPoint3D myPoint3D::operator/(double d)
{
	return myPoint3D(X / d, Y / d, Z / d);
}

myPoint3D myPoint3D::operator*(double d)
{
	return myPoint3D(X * d, Y * d, Z * d);
}

myVector3D myPoint3D::operator-(myPoint3D p1)
{
	return myVector3D(X - p1.X, Y - p1.Y, Z - p1.Z);
}

double myPoint3D::dist(myPoint3D p1)
{
	return sqrt((p1.X - X)*(p1.X - X) + (p1.Y - Y)*(p1.Y - Y) + (p1.Z - Z)*(p1.Z - Z));
}

void myPoint3D::rotate(myVector3D & lp, double theta)
{
	myVector3D tmp(X, Y, Z);
	tmp.rotate(lp, theta);
	X = tmp.dX; Y = tmp.dY; Z = tmp.dZ;
}

void myPoint3D::print(char *s)
{
	std::cout << s << X << ", " << Y << ", " << Z << "\n";
}

double myPoint3D::dist(myPoint3D *p1, myPoint3D *p2)
{
	//distance between current point, and the segment defined by p1,p2.
	/**** TODO ****/

	myVector3D p1p2 = myVector3D(p2->X - p1->X, p2->Y - p1->Y, p2->Z - p1->Z);
	myVector3D cp1 = myVector3D(p1->X - this->X, p1->Y - this->Y, p1->Z - this->Z);

	double dot = abs(p1p2 * cp1); // produit scalaire
	double t = dot / (p1p2.length() * p1p2.length()); // scalar projection
	t = std::max(0.0, std::min(1.0, t));
    // Calcul de la projection sur le segment
    myPoint3D proj = myPoint3D(p1->X + t * p1p2.dX, p1->Y + t * p1p2.dY, p1->Z + t * p1p2.dZ);

	double distP1 = this->dist(*p1); // la projection tombe avant P1
	double distP2 = this->dist(*p2); // la projection tombe après P2
	double distPROJ = this->dist(proj);

	return std::min(std::min(distP1, distP2), distPROJ); // on prend la plus petite distance
}

double myPoint3D::dist(myPoint3D *p1, myPoint3D *p2, myPoint3D *p3)
{
	//distance  between current point, and the triangled defined by p1,p2,p3.
	/**** TODO ****/
	return 0.0;
}

void myPoint3D::circumcenter(myPoint3D *p1, myPoint3D *p2, myPoint3D *p3, myPoint3D *p4)
{
	double xba, yba, zba, xca, yca, zca, xda, yda, zda;
	double balength, calength, dalength;
	double xcrosscd, ycrosscd, zcrosscd;
	double xcrossdb, ycrossdb, zcrossdb;
	double xcrossbc, ycrossbc, zcrossbc;
	double denominator;
	double xcirca, ycirca, zcirca;

	/* Use coordinates relative to point `a' of the tetrahedron. */
	xba = p2->X - p1->X;
	yba = p2->Y - p1->Y;
	zba = p2->Z - p1->Z;
	xca = p3->X - p1->X;
	yca = p3->Y - p1->Y;
	zca = p3->Z - p1->Z;
	xda = p4->X - p1->X;
	yda = p4->Y - p1->Y;
	zda = p4->Z - p1->Z;
	/* Squares of lengths of the edges incident to `a'. */
	balength = xba * xba + yba * yba + zba * zba;
	calength = xca * xca + yca * yca + zca * zca;
	dalength = xda * xda + yda * yda + zda * zda;
	/* Cross products of these edges. */
	xcrosscd = yca * zda - yda * zca;
	ycrosscd = zca * xda - zda * xca;
	zcrosscd = xca * yda - xda * yca;
	xcrossdb = yda * zba - yba * zda;
	ycrossdb = zda * xba - zba * xda;
	zcrossdb = xda * yba - xba * yda;
	xcrossbc = yba * zca - yca * zba;
	ycrossbc = zba * xca - zca * xba;
	zcrossbc = xba * yca - xca * yba;

	denominator = 0.5 / (xba * xcrosscd + yba * ycrosscd + zba * zcrosscd);

	/* Calculate offset (from `a') of circumcenter. */
	xcirca = (balength * xcrosscd + calength * xcrossdb + dalength * xcrossbc) *
		denominator;
	ycirca = (balength * ycrosscd + calength * ycrossdb + dalength * ycrossbc) *
		denominator;
	zcirca = (balength * zcrosscd + calength * zcrossdb + dalength * zcrossbc) *
		denominator;

	X = xcirca + p1->X;
	Y = ycirca + p1->Y;
	Z = zcirca + p1->Z;
}