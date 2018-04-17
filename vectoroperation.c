/***************************************************
* vector operations
* calculate the lattice surface area of the Bravais lattice
* spanned by vectors r1, r2, r3
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

//vmd unic cell parameters
// {a b c alpha beta gamma}
// alpha - bc beta - ca gamma - ab
// a should along x-axis
void BasisVector(VECTOR a,VECTOR b, VECTOR c, int n)
{
 Lx[n] = Norm(a);
 Ly[n] = Norm(b);
 Lz[n] = Norm(c);

 Anglez[n] = VectorAngle(a,b);
 Anglex[n] = VectorAngle(b,c);
 Angley[n] = VectorAngle(c,a);

 GetVertice(n);

 return;
}

void GetVertice(int n)
{
 vertice[n][0].x = 0.0;
 vertice[n][0].y = 0.0;
 vertice[n][0].z = 0.0;

 vertice[n][1] = Add(vertice[n][0],bravais[n][0]);
 vertice[n][2] = Add(vertice[n][0],bravais[n][1]);
 vertice[n][3] = Add(vertice[n][0],bravais[n][2]);

 vertice[n][4] = Add(vertice[n][1],vertice[n][3]);
 vertice[n][5] = Add(vertice[n][1],vertice[n][2]);
 vertice[n][6] = Add(vertice[n][3],vertice[n][2]);

 vertice[n][7] = Add(vertice[n][6],vertice[n][1]);

 return;
}

double VectorAngle(VECTOR r1,VECTOR r2) // return in degrees
{
 double theta12;

 theta12 = DotProduct(r1,r2)/(Norm(r1)*Norm(r2));

 theta12 = acos(theta12)/M_PI*180;

 return theta12;
}

VECTOR Equal(VECTOR r)
{
 VECTOR rr;

 rr.x = r.x;
 rr.y = r.y;
 rr.z = r.z;

 return rr;
}

double Norm(VECTOR r)
{
 double rnorm;

 rnorm = sqrt(SQR(r.x) + SQR(r.y) + SQR(r.z));

 return rnorm;
}

VECTOR ScalarProduct(double lambda, VECTOR r)
{
 VECTOR rr;

 rr.x = lambda*r.x;
 rr.y = lambda*r.y;
 rr.z = lambda*r.z;

 return rr;
}

double DotProduct(VECTOR r1, VECTOR r2)
{
 double rr;

 rr = r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;

 return rr;
}

VECTOR CrossProduct(VECTOR r1,VECTOR r2)
{
 VECTOR r1r2;

 r1r2.x = r1.y*r2.z - r1.z*r2.y;
 r1r2.y = r1.z*r2.x - r1.x*r2.z;
 r1r2.z = r1.x*r2.y - r1.y*r2.x;

 return r1r2;
}

VECTOR Add(VECTOR r1,VECTOR r2)
{
 VECTOR r1r2;

 r1r2.x = r1.x + r2.x;
 r1r2.y = r1.y + r2.y;
 r1r2.z = r1.z + r2.z;

 return r1r2;
}

VECTOR Subtract(VECTOR r1,VECTOR r2)
{
 VECTOR r1r2;

 r1r2.x = r1.x - r2.x;
 r1r2.y = r1.y - r2.y;
 r1r2.z = r1.z - r2.z;

 return r1r2;
}


double CrossNorm(VECTOR r1,VECTOR r2)
{
 VECTOR r1r2;
 double r1r2norm;

 r1r2.x = r1.y*r2.z - r1.z*r2.y;
 r1r2.y = r1.z*r2.x - r1.x*r2.z;
 r1r2.z = r1.x*r2.y - r1.y*r2.x;

 r1r2norm = SQR(r1r2.x) + SQR(r1r2.y) + SQR(r1r2.z);

 r1r2norm = sqrt(r1r2norm);

 return r1r2norm;
}

double LatticeArea(VECTOR r1,VECTOR r2,VECTOR r3)
{ 
 double area;

 area  = CrossNorm(r1,r2) + CrossNorm(r2,r3) + CrossNorm(r3,r1);

 return area;
}

double LatticeVolume(VECTOR r1,VECTOR r2,VECTOR r3)
{ 
 double v;

 v = DotProduct(r1,CrossProduct(r2,r3));
 v = fabs(v);

 return v;
}


