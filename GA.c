/*****************************************************************
* Genetic Algorithm (GA)
* find the optimum packing of 3D ternary hard sphere crystals
* Kai Zhang, Yale University, 2013
*****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

//int main(void)
int main(int argc, char *argv[])
{
 //int step;
 int m,n,i,d,ns;
 char filename[20];
 FILE *fp[MAX_PopulationSize];
 FILE *fps[MAX_PopulationSize];
 FILE *fpvert[MAX_PopulationSize];
 FILE *fpgene[MAX_PopulationSize];
 FILE *fpgenefinal[MAX_PopulationSize];

 double Count;

 int ngen;
 double minfit,maxfit;
 int fittestn;
 int p1,p2,k1,k2;
 VECTOR B1[MAX_NumberOfParticles],B2[MAX_NumberOfParticles],Bchild[MAX_NumberOfParticles]; // displacement vectors
 VECTOR L1[MAX_Dimensions],L2[MAX_Dimensions],Lchild[MAX_Dimensions]; // projection vectors
 double norm;
 int colinear;
 int betterfit;
 int clonecheck;
 VECTOR shiftk2;
 VECTOR hydrai,hydraj;
 double dotprod;
 int ntrials;

 double dE;
 double tol;

 double dfit;
 dE =  1.0e-4;
 tol = 1.0e-10;

 sscanf(argv[1],"%d",&JobIndex);
 
 printf("**************** 3D Hard Sphere Cystal Genetic Algorithm ********************");
 printf("\n");

 ReadInput();

 if(randomseed == 0.)  randomseed = (double)(JobIndex);
 printf("randomseed = %lf\n",randomseed);
 InitializeRandomNumberGenerator(randomseed); //time(0L)

 Initialization(JobIndex);
 fflush(stdout);

 for(n=0;n<PopulationSize;n++)
 {
  sprintf(filename,"movie_%d_%d.xyz",JobIndex,n);
  fp[n]=fopen(filename,"a+");
  sprintf(filename,"vertice_%d_%d.dat",JobIndex,n);
  fpvert[n]=fopen(filename,"a+");
  sprintf(filename,"genetics_%d_%d.dat",JobIndex,n);
  fpgene[n]=fopen(filename,"a+");
  sprintf(filename,"genetics_%d",n);
  fpgenefinal[n]=fopen(filename,"w");
  Writemovie(fp[n],n);
  Writegene(fpgene[n],n);
  Writegene(fpgenefinal[n],n);
  WriteVertice(fpvert[n],n);
  fclose(fp[n]);
  fclose(fpgene[n]);
  fclose(fpgenefinal[n]);
  fclose(fpvert[n]);
  
  if(InitialType != 0)
  {
   LatticeReduction(n);
   Minimize(n);
  }
  else
  {
  ForceLattice(n); // force[i]
  fitness[n] = -Upotential;
  }
 }//end loop n
 
 printf("start GA ......\n");
 printf("\n");

for(ngen=0;ngen<NGeneration;ngen++)
{
 for(n=0;n<PopulationSize;n++)
 {
  Store(n,n+PopulationSize);
  RotationXY(n);
  sprintf(filename,"movie_%d_%d.xyz",JobIndex,n);
  fp[n]=fopen(filename,"a+");
  sprintf(filename,"vertice_%d_%d.dat",JobIndex,n);
  fpvert[n]=fopen(filename,"a+");
  sprintf(filename,"genetics_%d_%d.dat",JobIndex,n);
  fpgene[n]=fopen(filename,"a+");
  sprintf(filename,"genetics_%d",n);
  fpgenefinal[n]=fopen(filename,"w");
  Writemovie(fp[n],n);
  Writegene(fpgene[n],n);
  Writegene(fpgenefinal[n],n);
  WriteVertice(fpvert[n],n);
  fclose(fp[n]);
  fclose(fpgene[n]);
  fclose(fpgenefinal[n]);
  fclose(fpvert[n]);
  Store(n+PopulationSize,n);

  volume[n] = LatticeVolume(bravais[n][0],bravais[n][1],bravais[n][2]);
  packfrac[n] = M_PI/6.*(NA*CUBIC(sigmaA)+NB*CUBIC(sigmaB)+NC*CUBIC(sigmaC))/volume[n];
  printf("n = %d fitness = %lf volume = %lf packfrac = %lf\n",n,fitness[n],volume[n],packfrac[n]);
 }//end loop over population

 fflush(stdout);

 minfit = 100000.;
 maxfit = -100000.;
 for(n=0;n<PopulationSize;n++)
 {
  if(fitness[n] < minfit)
   minfit = fitness[n];
  if(fitness[n] > maxfit)
   maxfit = fitness[n];
 }

 printf("generation = %d minfit = %lf maxfit = %lf\t",ngen,minfit,maxfit);

 for(n=0;n<PopulationSize;n++)
 {
  if(ngen > 9) fitnesshistory[10][n] = fitnesshistory[9][n];
  if(ngen > 8) fitnesshistory[9][n] = fitnesshistory[8][n];
  if(ngen > 7) fitnesshistory[8][n] = fitnesshistory[7][n];
  if(ngen > 6) fitnesshistory[7][n] = fitnesshistory[6][n];
  if(ngen > 5) fitnesshistory[6][n] = fitnesshistory[5][n];
  if(ngen > 4) fitnesshistory[5][n] = fitnesshistory[4][n];
  if(ngen > 3) fitnesshistory[4][n] = fitnesshistory[3][n];
  if(ngen > 2) fitnesshistory[3][n] = fitnesshistory[2][n];
  if(ngen > 1) fitnesshistory[2][n] = fitnesshistory[1][n];
  if(ngen > 0) fitnesshistory[1][n] = fitnesshistory[0][n];
  fitnesshistory[0][n] = fitness[n];
 }
 dfit = 0.;
 for(n=0;n<PopulationSize;n++)
  dfit += fitness[n] - fitnesshistory[10][n];
 printf("dfit = %lf\n",dfit);
 if(dfit < dE)
 //if(maxfit-minfit < dE)
 {
  printf("converged\n");
  printf("generation = %d minfit = %lf maxfit = %lf dfit = %lf\n",ngen,minfit,maxfit,dfit);
  break;
 // return 0;
 }
 
fflush(stdout);

for(ns=0;ns<OffspringNumber;ns++)
{
ntrials=0;
do
{
 ntrials++;
 // pick up two parents
 p1 = (int)(RandomNumber()*PopulationSize);
 do
 p2 = (int)(RandomNumber()*PopulationSize);
 while(p2 == p1);
 
 // prepare parents displacement vectors, with mutation
 for(i=0;i<NumberOfParticles-1;i++)
 {
  B1[i] = Subtract(position[p1][i+1],position[p1][i]); //Bi is vector from i to i+1
  if(RandomNumber() < mutationrate) // mutation of vector length
  {
 // norm = BoxMuller(0.6*(sigma[i+1]+sigma[i]),0.1*(sigma[i+1]+sigma[i]));
  norm = BoxMuller(1.0,0.1);
  B1[i] = ScalarProduct(norm*Norm(B1[i]),B1[i]);
  }

  B2[i] = Subtract(position[p2][i+1],position[p2][i]); //Bi is vector from i to i+1
  if(RandomNumber() < mutationrate) // mutation of vector length
  {
 // norm = BoxMuller(0.6*(sigma[i+1]+sigma[i]),0.1*(sigma[i+1]+sigma[i]));
  norm = BoxMuller(1.0,0.1);
  B2[i] = ScalarProduct(norm*Norm(B2[i]),B2[i]);
  }
 }// end loop over particles
 
 //crossover of projection vectors
 do
 {
 k1 = (int)(RandomNumber()*(Dimension+1)); // 0, 1, 2, 3
 //k1 = (int)(RandomNumber()*(Dimension)); // 0, 1, 2
 if(k1 == 3)
 {
 Lchild[0] = Equal(projection[p2][0]);
 Lchild[1] = Equal(projection[p2][1]);
 Lchild[2] = Equal(projection[p2][2]);
 }
 else if(k1 == 0)
 {
 Lchild[0] = Equal(projection[p1][0]);
 Lchild[1] = Equal(projection[p2][1]);
 Lchild[2] = Equal(projection[p2][2]);
 }
 else if (k1 == 1)
 {
 Lchild[0] = Equal(projection[p1][0]);
 Lchild[1] = Equal(projection[p1][1]);
 Lchild[2] = Equal(projection[p2][2]);
 }
 else
 {
 Lchild[0] = Equal(projection[p1][0]);
 Lchild[1] = Equal(projection[p1][1]);
 Lchild[2] = Equal(projection[p1][2]);
 }
 // test colinear
 colinear = 0;
 if(CrossNorm(Lchild[0],Lchild[1])<tol || CrossNorm(Lchild[1],Lchild[2])<tol || CrossNorm(Lchild[2],Lchild[0])<tol )
 colinear = 1;
 }
 while(colinear == 1);
  
 /*
//mutation of projection vector
 if(RandomNumber() < mutationrate) // mutation of vector length
 {
 norm = Norm(Lchild[0]);
 Lchild[0] = ScalarProduct(BoxMuller(norm,0.1*norm)/norm,Lchild[0]);
 } 
 if(RandomNumber() < mutationrate) // mutation of vector length
 {
 norm = Norm(Lchild[1]);
 Lchild[1] = ScalarProduct(BoxMuller(norm,0.1*norm)/norm,Lchild[1]);
 }
 if(RandomNumber() < mutationrate) // mutation of vector length
 {
 norm = Norm(Lchild[2]);
 Lchild[2] = ScalarProduct(BoxMuller(norm,0.1*norm)/norm,Lchild[2]);
 }
*/
 
//generate new child
/********************************************/
 //crossover of displacement vectors
 if(NumberOfParticles > 3)
 k2 = (int)(RandomNumber()*(NumberOfParticles-2)); // 0, 1, 2, ..., N-3
 else
 k2 = 0; 
 for(i=0;i<=k2;i++)
 Bchild[i] = Equal(B1[i]);
 for(i=k2+1;i<NumberOfParticles-1;i++)
 Bchild[i] = Equal(B2[i]);

 //projection vector
 for(d=0;d<Dimension;d++)
  projection[ns+PopulationSize][d]= Equal(Lchild[d]);
 
 //position and displacement vector
 position[ns+PopulationSize][0].x = 0.;
 position[ns+PopulationSize][0].y = 0.;
 position[ns+PopulationSize][0].z = 0.;
 for(i=0;i<NumberOfParticles-1;i++)
 {
  displacement[ns+PopulationSize][i] = Equal(Bchild[i]);
  position[ns+PopulationSize][i+1] = Add(position[ns+PopulationSize][i],displacement[ns+PopulationSize][i]);
 } 
 
 //bravais
 for(d=0;d<Dimension;d++)
 {
  bravais[ns+PopulationSize][d] = Equal(projection[ns+PopulationSize][d]);
  for(i=ito[d];i<ifrom[d];i++)
  bravais[ns+PopulationSize][d] = Add(bravais[ns+PopulationSize][d],displacement[ns+PopulationSize][i]);
 }
/********************************************/
//if(OverlapLattice(ns+PopulationSize) == 1)
if(0> 1)
{
 betterfit = 0;
}
else
{

//compress the longest lattice vector
if(CompressSwitch == 1) Compress(ns+PopulationSize);

/********************************************/
if(RotationSwitch == 1)
 MinimizeEuler((int)(PopulationSize+ns));

if(RotationSwitch == 2)
{
 //shift the joint point k2+1 to origin for Euler rotation
 shiftk2 = Equal(position[ns+PopulationSize][k2+1]);
 for(i=0;i<NumberOfParticles;i++)
 position[ns+PopulationSize][i] = Subtract(position[ns+PopulationSize][i],shiftk2);
 MinimizeEulerTwo((int)(PopulationSize+ns),k2);
 //ForceLattice(PopulationSize+ns);
// fitness[PopulationSize+ns] = -Upotential;
}
/********************************************/

 LatticeReduction((int)(PopulationSize+ns));
 Minimize((int)(PopulationSize+ns)); // offspring info is temporarily stored in array[PopulationSize+i]

 volume[ns+PopulationSize] = LatticeVolume(bravais[ns+PopulationSize][0],bravais[ns+PopulationSize][1],bravais[ns+PopulationSize][2]);
 packfrac[ns+PopulationSize] = M_PI/6.*(NA*CUBIC(sigmaA)+NB*CUBIC(sigmaB)+NC*CUBIC(sigmaC))/volume[ns+PopulationSize];
//check better fitness
/*************************************************/
 betterfit = 1;
 if(packfrac[ns+PopulationSize] > 0.8)
  betterfit = 0;

/*
 if(betterfit == 1)
 {
 for(i=0;i<NumberOfParticles-1;i++)
 {
  if(Norm(displacement[PopulationSize+ns][i]) > distancemax)
  {
  betterfit = 0;
  break;
  }
 }
 }
 if(betterfit == 1)
 {
 for(d=0;d<Dimension;d++)
 {
  if(Norm(bravais[PopulationSize+ns][d]) > NumberOfParticles*sigmaA)
  {
  betterfit = 0;
  break;
  }
 }
 }
*/
 if(betterfit == 1)
 {
 /*
 for(n=0;n<PopulationSize;n++) // compare with parents
 {
  if(fitness[n] > fitness[PopulationSize+ns])
  {
  betterfit = 0;
  break;
  }
 }
*/
  if(minfit > fitness[PopulationSize+ns])
  betterfit = 0;
 }
/*************************************************/

//check diversity
/*************************************************/
if(betterfit == 1 && similaritythreshold < 1.0)
{
 clonecheck = 1;
 for(n=0;n<PopulationSize+ns;n++) 
 {
  dotprod = 0.;
  for(i=0;i<NumberOfParticles-1;i++)
  {
   hydrai = ScalarProduct(1./Norm(displacement[ns+PopulationSize][i]),displacement[ns+PopulationSize][i]);
   hydraj = ScalarProduct(1./Norm(displacement[n][i]),displacement[n][i]);
   dotprod += DotProduct(hydrai,hydraj);
  }
  for(d=0;d<Dimension;d++)
  {
   hydrai = ScalarProduct(1./Norm(projection[ns+PopulationSize][d]),projection[ns+PopulationSize][d]);
   hydraj = ScalarProduct(1./Norm(projection[n][d]),projection[n][d]);
   dotprod += DotProduct(hydrai,hydraj);
  }
  dotprod = dotprod/(NumberOfParticles-1+Dimension);
  
  if(dotprod >similaritythreshold)
  {
  clonecheck = 0; 
  if(n>=PopulationSize && fitness[n] < fitness[PopulationSize+ns])
  {
 //  printf("offsping %d with %lf -> ",n,fitness[n]);
   Store(PopulationSize+ns,n);
   fitness[n] = fitness[PopulationSize+ns];
 //  printf("%d with %lf \n",PopulationSize+ns,fitness[n]);
  }
  break;
  }
 }//end loop n
}
/*************************************************/

}//end else no overlap

}// end loop do if not better fit
while(betterfit == 0 || clonecheck == 0);

printf("ns = %d fitness = %lf\tntrials = %d\n",ns,fitness[PopulationSize+ns],ntrials);
fflush(stdout);

}//end loop over offsprings
 
printf("\n");

//update generation
/*************************************************/
for(m=2*PopulationSize;m<3*PopulationSize;m++)
{
 maxfit = -1000000000.;
 for(n=0;n<PopulationSize+OffspringNumber;n++)
 {
  if(fitness[n] > maxfit)
  {
   maxfit = fitness[n];
   fittestn = n;
  }
 }
 Store(fittestn,m);
 fitness[m] = fitness[fittestn];
 fitness[fittestn] = -10000000000.;
}

for(n=0;n<PopulationSize;n++)
{
 Store(n+2*PopulationSize,n);
 fitness[n] = fitness[n+2*PopulationSize];
}
/*************************************************/

 
//update generation
/*
for(m=PopulationSize+OffspringNumber;m<2*PopulationSize;m++)
{
 maxfit = -10000;
 for(n=0;n<PopulationSize;n++)
 {
  if(fitness[n] > maxfit)
  {
   maxfit = fitness[n];
   fittestn = n;
  }
 }
 Store(fittestn,m);
 fitness[m] = fitness[fittestn];
 fitness[fittestn] = -10000;
}//end loop m

for(n=0;n<PopulationSize;n++)
{
 Store(n+PopulationSize,n);
 fitness[n] = fitness[n+PopulationSize];
}
*/

} // end loop over generations

 printf("\n");

 printf("****************************** the end *******************************");


 return 0;
}
