/***************************************************
* rotate the coordinates such that Bravais vector 0
* is along x-axis for plot purpose
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"

//rotate the coordinates individual n
//with respect to unit vector rot
// counterclockwise by angle theta
void Rotation(int n,VECTOR rot,double theta)
{
 int i,d;
 double x0,y0,z0;

 A.xx = cos(theta) + SQR(rot.x)*(1.0-cos(theta));
 A.xy = rot.x*rot.y*(1.0-cos(theta)) - rot.z*sin(theta);
 A.xz = rot.x*rot.z*(1.0-cos(theta)) + rot.y*sin(theta);
 A.yx = rot.y*rot.x*(1.0-cos(theta)) + rot.z*sin(theta);
 A.yy = cos(theta) + SQR(rot.y)*(1.0-cos(theta));
 A.yz = rot.y*rot.z*(1.0-cos(theta)) - rot.x*sin(theta);
 A.zx = rot.z*rot.x*(1.0-cos(theta)) - rot.y*sin(theta);
 A.zy = rot.z*rot.y*(1.0-cos(theta)) + rot.x*sin(theta);
 A.zz = cos(theta) + SQR(rot.z)*(1.0-cos(theta));

 for(d=0;d<Dimension;d++)
 {
  x0 = bravais[n][d].x;
  y0 = bravais[n][d].y;
  z0 = bravais[n][d].z;
  bravais[n][d].x = A.xx*x0 + A.xy*y0 + A.xz*z0;
  bravais[n][d].y = A.yx*x0 + A.yy*y0 + A.yz*z0;
  bravais[n][d].z = A.zx*x0 + A.zy*y0 + A.zz*z0;  
  
  x0 = projection[n][d].x;
  y0 = projection[n][d].y;
  z0 = projection[n][d].z;
  projection[n][d].x = A.xx*x0 + A.xy*y0 + A.xz*z0;
  projection[n][d].y = A.yx*x0 + A.yy*y0 + A.yz*z0;
  projection[n][d].z = A.zx*x0 + A.zy*y0 + A.zz*z0;  
 }

 for(i=0;i<NumberOfParticles;i++)
 {
  x0 = position[n][i].x;
  y0 = position[n][i].y;
  z0 = position[n][i].z;
  position[n][i].x = A.xx*x0 + A.xy*y0 + A.xz*z0;
  position[n][i].y = A.yx*x0 + A.yy*y0 + A.yz*z0;
  position[n][i].z = A.zx*x0 + A.zy*y0 + A.zz*z0;  
  
  if(i>0)
  displacement[n][i-1] = Subtract(position[n][i],position[n][i-1]);
 }

 return;
}

// rotate the coordinates such that Bravais vector 0 is along x-axis
// and Bravais vector 1 is in the xy plane
void RotationXY(int n)
{
 VECTOR rot1,rot2;
 double theta1,theta2;
 double norm;
 VECTOR xaxis,yaxis,yzproj;
 int i;

 if(bravais[n][0].y != 0. || bravais[n][0].z != 0.)
 {
 xaxis.x = 1.0;
 xaxis.y = 0.;
 xaxis.z = 0.;
 
 theta1 = VectorAngle(xaxis,bravais[n][0]);
 theta1 = theta1/180.*M_PI;

 rot1.x = 0.;
 rot1.y = bravais[n][0].z;
 rot1.z = -bravais[n][0].y;
 norm = sqrt(SQR(rot1.y) + SQR(rot1.z));
 rot1.y = rot1.y/norm; // nan if bravais0 is already along x
 rot1.z = rot1.z/norm;

 // align x-axis
 Rotation(n,rot1,theta1);
 }

 if(bravais[n][1].z != 0.) // if not in the xy plane
 {
 yaxis.x = 0.;
 yaxis.y = 1.;
 yaxis.z = 0.;

 yzproj.x = 0.;
 yzproj.y = bravais[n][1].y;
 yzproj.z = bravais[n][1].z;
 
 theta2 = VectorAngle(yaxis,yzproj);
 theta2 = theta2/180.*M_PI;

 if(bravais[n][1].z > 0.)
 rot2.x = -1.0;
 else
 rot2.x = 1.0;

 rot2.y = 0.;
 rot2.z = 0.;

 Rotation(n,rot2,theta2);
 }

 if(bravais[n][2].z < 0.) // if not right-hand frame
 {
 bravais[n][2] = ScalarProduct(-1.,bravais[n][2]); 
 
 projection[n][2] = Equal(bravais[n][2]);
 for(i=ito[2];i<ifrom[2];i++)
 projection[n][2] = Subtract(projection[n][2],displacement[n][i]);
 }

 BasisVector(bravais[n][0],bravais[n][1],bravais[n][2],n);

 return;
}

//only rotate vector r
void RotationVector(VECTOR rot,double theta,VECTOR *r)
{
 double x0,y0,z0;

 A.xx = cos(theta) + SQR(rot.x)*(1.0-cos(theta));
 A.xy = rot.x*rot.y*(1.0-cos(theta)) - rot.z*sin(theta);
 A.xz = rot.x*rot.z*(1.0-cos(theta)) + rot.y*sin(theta);
 A.yx = rot.y*rot.x*(1.0-cos(theta)) + rot.z*sin(theta);
 A.yy = cos(theta) + SQR(rot.y)*(1.0-cos(theta));
 A.yz = rot.y*rot.z*(1.0-cos(theta)) - rot.x*sin(theta);
 A.zx = rot.z*rot.x*(1.0-cos(theta)) - rot.y*sin(theta);
 A.zy = rot.z*rot.y*(1.0-cos(theta)) + rot.x*sin(theta);
 A.zz = cos(theta) + SQR(rot.z)*(1.0-cos(theta));

  x0 = (*r).x;
  y0 = (*r).y;
  z0 = (*r).z;

  (*r).x = A.xx*x0 + A.xy*y0 + A.xz*z0;
  (*r).y = A.yx*x0 + A.yy*y0 + A.yz*z0;
  (*r).z = A.zx*x0 + A.zy*y0 + A.zz*z0;  

 return;
}

