#include "system.h"

double randomseed;
int JobIndex;
int Dimension;

int InitialType,IDtype;

int NGeneration;
int OffspringNumber;
int PopulationSize; // number of individuals
int NumberOfParticles; //number of particles per unit cell
int NA,NB,NC; // NA + NB + NC = N
double fA,fB,fC;// fraction

double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
double sigmaA,sigmaB,sigmaC,sigmaAB,sigmaBC,sigmaCA;

double rho,rhoA,rhoB,rhoC; // number density rho = N/V
double packingfraction; // phi = pi/6 rho in 3D
double V,L;  // Volume
int NumberOfLatticeSites; 
int NCells;
double mutationrate,similaritythreshold,distancemax,compressrate;
int RotationSwitch,CompressSwitch;

int PotentialType;
double n_rep,n_att,epsilonfactor; // replulsion and attraction index
double epsilonA,epsilonB,epsilonC,epsilonAB,epsilonBC,epsilonCA;
double coreA,coreB,coreC,coreAB,coreBC,coreCA;

double rc;// rcutoff
double rcA,rcB,rcC,rcAB,rcBC,rcCA;
double ucA,ucB,ucC,ucAB,ucBC,ucCA; // u(r=rc) at cutoff
double ducdrA,ducdrB,ducdrC,ducdrAB,ducdrBC,ducdrCA; // du(r=rc)/dr at cutoff

MATRIX A;// rotation matrix

VECTOR position[3*MAX_PopulationSize][MAX_NumberOfParticles]; // position of particles
VECTOR displacement[3*MAX_PopulationSize][MAX_NumberOfParticles]; // displacement vectors
VECTOR projection[3*MAX_PopulationSize][MAX_Dimensions]; // lattice vectors

VECTOR bravais[3*MAX_PopulationSize][MAX_Dimensions]; // bravais lattice vectors
double Lx[3*MAX_PopulationSize],Ly[3*MAX_PopulationSize],Lz[3*MAX_PopulationSize],Anglex[3*MAX_PopulationSize],Angley[3*MAX_PopulationSize],Anglez[3*MAX_PopulationSize];

double fitness[3*MAX_PopulationSize]; // energy
double fitnesshistory[100][3*MAX_PopulationSize]; // -energy Congen generations before
double volume[3*MAX_PopulationSize]; // volume of unit cell
double packfrac[3*MAX_PopulationSize]; // packing fraction

int ifrom[MAX_Dimensions];
int ito[MAX_Dimensions];

double ministep; // energy minimization step size
int MaxIter;// max iterations
int Congen; // number of generations beyond with no improvement then converge
double converge; // convergence criterion
VECTOR force[MAX_NumberOfParticles]; // force on particle i
VECTOR hstep[MAX_NumberOfParticles]; // displacement in conjugate gradient
double Upotential,Ftotal;

int globaln,globalk2;
VECTOR rotaxis;

VECTOR vertice[3*MAX_PopulationSize][8]; // vertices of unit cell box
