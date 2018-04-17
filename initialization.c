/***************************************************
* initialize particle parameters
****************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void Initialization(int job)
{ 
 FILE *fp;
 char filename[20];
 int n,i,j,d;
 int nx,ny,nz;
 int overlap;
 double rij2;

 FILE *fpid;

 double xmaxp,ymaxp,zmaxp; // positive
 double xmaxn,ymaxn,zmaxn; // negative
 VECTOR ba,bb,bc; // the second and third Bravais lattice vector
 double phib,phic;

 double phi,sintheta,costheta,dr;
 double dx,dy,dz;

 if(IDtype == 0)
 {
 for(i=0;i<NumberOfParticles;i++)
 {
  if(i<NA)
  {
  sigma[i] = sigmaA;
  identity[i] = 1;
  }
  else if(i >= NA && i< NA+NB)
  {
  sigma[i] = sigmaB;
  identity[i] = 2;
  }
  else
  {
  sigma[i] = sigmaC;
  identity[i] = 3;
  }
 }
 }
 else
 {
  fpid = fopen("inputid","r");
  for(i=0;i<NumberOfParticles;i++)
  {
   fscanf(fpid,"%lf",&identity[i]);
   if(identity[i] == 1)
    sigma[i] = sigmaA;
   else if(identity[i] == 2)
    sigma[i] = sigmaB;
   else
    sigma[i] = sigmaC;
  }
  fclose(fpid);
 }

if(InitialType == 0) // from input file 
{
for(n=0;n<PopulationSize;n++)
{
sprintf(filename,"genetics_%d",n);
fp = fopen(filename,"r");
for(i=0;i<NumberOfParticles+Dimension;i++)
{
if(i<NumberOfParticles)
fscanf(fp,"%lf\t%lf\t%lf\n",&position[n][i].x,&position[n][i].y,&position[n][i].z);
else
fscanf(fp,"%lf\t%lf\t%lf\n",&projection[n][i-NumberOfParticles].x,&projection[n][i-NumberOfParticles].y,&projection[n][i-NumberOfParticles].z);
}
fclose(fp);
}
}// end input from file


if(InitialType == 2)
{
for(n=0;n<PopulationSize;n++)
{
 position[n][0].x = 0.;
 position[n][0].y = 0.;
 position[n][0].z = 0.;
 
 position[n][1].x = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 position[n][1].y = 0.;
 position[n][1].z = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 
 position[n][2].x = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 position[n][2].y = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 position[n][2].z = 0.;
 
 position[n][3].x = 0.;
 position[n][3].y = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 position[n][3].z = sigmaA/sqrt(2.)*(1.0 + RandomNumber()*0.1);
 
 bravais[n][0].x = sqrt(2.)*sigmaA*(1.2 + RandomNumber()*0.2);
 bravais[n][0].y = 0.;
 bravais[n][0].z = 0.;
 
 bravais[n][1].x = 0. + RandomNumber()*0.1;
 bravais[n][1].y = sqrt(2.)*sigmaA*(1.2 + RandomNumber()*0.2);
 bravais[n][1].z = 0.;
 
 bravais[n][2].x = 0. + RandomNumber()*0.1;
 bravais[n][2].y = 0. + RandomNumber()*0.1;
 bravais[n][2].z = sqrt(2.)*sigmaA*(1.2 + RandomNumber()*0.2);

 BasisVector(bravais[n][0],bravais[n][1],bravais[n][2],n);
}
}//end if FCC initial


if(InitialType == 1)
{
for(n=0;n<PopulationSize;n++)
{
 position[n][0].x = 0.;
 position[n][0].y = 0.;
 position[n][0].z = 0.;

do
{
 xmaxp = 0.;
 ymaxp = 0.;
 zmaxp = 0.;
 
 xmaxn = 0.;
 ymaxn = 0.;
 zmaxn = 0.;
  
 phib = M_PI/1.*RandomNumber();
 phic = M_PI/1.*RandomNumber();
 
// phib = M_PI/2. - M_PI/16.*RandomNumber();
// phic = M_PI/2. - M_PI/16.*RandomNumber();

//phib = M_PI/2. - M_PI/10.*RandomNumber();
// phic = phib/2;

 costheta = RandomNumber(); // theta in [0, pi/2]
 //costheta = 1/sqrt(2.)*RandomNumber(); // theta in [0, pi/4]
 //costheta = 1 - (1.-1/sqrt(2.))*RandomNumber(); // theta in [0, pi/4]
// costheta = 1 - (0.05)*RandomNumber(); // theta in [0, pi/4]
 //costheta = 0.18*RandomNumber(); // theta in [0, pi/4]

 sintheta = sqrt(1.-SQR(costheta));

 ba.x = 1.0;
 ba.y = 0.;
 ba.z = 0.;

 bb.x = cos(phib);
 bb.y = sin(phib);
 bb.z = 0.;

 bc.x = sintheta*cos(phic);
 bc.y = sintheta*sin(phic);
 bc.z = costheta;

 for(i=1;i<NumberOfParticles;i++)
 {
  do
  {
  phi = 2.*M_PI*RandomNumber();
  costheta = 2.0*(RandomNumber()-0.5);
  sintheta = sqrt(1.-SQR(costheta));
  dr = BoxMuller(0.6*(sigma[i-1]+sigma[i]),0.1*(sigma[i-1]+sigma[i]));

  dx = dr*sintheta*cos(phi);
  dy = dr*sintheta*sin(phi);
  dz = dr*costheta;

  position[n][i].x = position[n][i-1].x + dx;
  position[n][i].y = position[n][i-1].y + dy;
  position[n][i].z = position[n][i-1].z + dz;

  overlap = 0;

 // if(position[n][i].x <0. || position[n][i].y<0. || position[n][i].z <0.) 
 //  overlap = 1;
 // else
  for(j=0;j<i;j++)
  {
   rij2 = Distance(n,i,j);
   if(rij2 <= SQR(SigmaIJ(identity[i],identity[j])))
   {
   overlap = 1;
   break;
   }
  }
  }
  while(overlap == 1);

if(DotProduct(position[n][i],ba) >= 0.) 
{
  if(DotProduct(position[n][i],ba) > xmaxp) xmaxp = DotProduct(position[n][i],ba); 
}
else
  if(fabs(DotProduct(position[n][i],ba)) > xmaxn) xmaxn = fabs(DotProduct(position[n][i],ba)); 

if(DotProduct(position[n][i],bb) >= 0.) 
{
  if(DotProduct(position[n][i],bb) > ymaxp) ymaxp = DotProduct(position[n][i],bb); 
}
else
  if(fabs(DotProduct(position[n][i],bb)) > ymaxn) ymaxn = fabs(DotProduct(position[n][i],bb)); 

if(DotProduct(position[n][i],bc) >= 0.) 
{
  if(DotProduct(position[n][i],bc) > zmaxp) zmaxp = DotProduct(position[n][i],bc); 
}
else
  if(fabs(DotProduct(position[n][i],bc)) > zmaxn) zmaxn = fabs(DotProduct(position[n][i],bc)); 
 
 }// end loop over particles

 bravais[n][0] = ScalarProduct(xmaxp+xmaxn+1.0*sigmaA+RandomNumber()*sigmaA*0.1,ba);
 bravais[n][1] = ScalarProduct(ymaxp+ymaxn+1.0*sigmaA+RandomNumber()*sigmaA*0.1,bb);
 bravais[n][2] = ScalarProduct(zmaxp+zmaxn+1.0*sigmaA+RandomNumber()*sigmaA*0.1,bc);

}
while(OverlapLattice(n) == 1);
 
// get a b c alpha beta gamma
 BasisVector(bravais[n][0],bravais[n][1],bravais[n][2],n);

}//end loop over individuals
} //end if random initial


for(n=0;n<PopulationSize;n++)
 for(i=0;i<NumberOfParticles-1;i++)
  displacement[n][i] = Subtract(position[n][i+1],position[n][i]); //Bi is vector from i to i+1

//ito < ifrom 
if(NumberOfParticles >= 6)
{
ifrom[0] = NumberOfParticles-1;
ito[0] = 0;

ifrom[1] = NumberOfParticles-2;
ito[1] = 1;

ifrom[2] = NumberOfParticles-3;
ito[2] = 2;
}
else if (NumberOfParticles == 5)
{
ifrom[0] = 4;
ito[0] = 0;

ifrom[1] = 3;
ito[1] = 1;

ifrom[2] = 2;
ito[2] = 2;
}
else if(NumberOfParticles == 4)
{
ifrom[0] = 3;
ito[0] = 0;

ifrom[1] = 2;
ito[1] = 1;

ifrom[2] = 2;
ito[2] = 2;
}
else if (NumberOfParticles == 3)
{
ifrom[0] = 0;
ito[0] = 0;

ifrom[1] = 1;
ito[1] = 1;

ifrom[2] = 2;
ito[2] = 2;
}
else if (NumberOfParticles == 2)
{
ifrom[0] = 0;
ito[0] = 0;

ifrom[1] = 1;
ito[1] = 1;

ifrom[2] = 1;
ito[2] = 0;
}
else if (NumberOfParticles == 1)
{
ifrom[0] = 0;
ito[0] = 0;

ifrom[1] = 0;
ito[1] = 0;

ifrom[2] = 0;
ito[2] = 0;
}

printf("projection vector (La,Lb,Lc): ");
for(d=0;d<Dimension;d++)
printf(" %d->%d ",ifrom[d],ito[d]);
printf("\n");

if(InitialType == 0) // projection -> bravais
{
for(n=0;n<PopulationSize;n++)
{
for(d=0;d<Dimension;d++)
{
  bravais[n][d] = Equal(projection[n][d]);
  for(i=ito[d];i<ifrom[d];i++)
  bravais[n][d] = Add(bravais[n][d],displacement[n][i]);
}
 BasisVector(bravais[n][0],bravais[n][1],bravais[n][2],n);
}
}
else // bravais -> projection
{
for(n=0;n<PopulationSize;n++)
{
for(d=0;d<Dimension;d++)
{
 projection[n][d] = Equal(bravais[n][d]);
 for(i=ito[d];i<ifrom[d];i++)
 projection[n][d] = Subtract(projection[n][d],displacement[n][i]);
}
}
}

 
return;
}
