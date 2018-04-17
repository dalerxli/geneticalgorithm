/***************************************
 * write movie snapshot in .xyz format
 *
 * number of particles
 * blank
 * name x y z
 * C 0.0000 0.0000 0.0000
 * S 0.0000 1.0000 0.2345
 *
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"

void Writemovie(FILE *FilePtr,int n)
{
  int i;
  
  fprintf(FilePtr,"%d\n",NumberOfParticles);
//  fprintf(FilePtr,"\n");
  fprintf(FilePtr,"%lf %lf %lf %lf %lf %lf\n",Lx[n],Ly[n],Lz[n],Anglex[n],Angley[n],Anglez[n]);

  for(i=0;i<NumberOfParticles;i++)
  {
     //if(i%100 == 0)	
     //if(i==500 || i==1000 || i==1500 ||i==2000)	
     if(identity[i] == 1)	
      fprintf(FilePtr,"%s\t","N");
     else if(identity[i] == 2)	
      fprintf(FilePtr,"%s\t","S");
     else
      fprintf(FilePtr,"%s\t","O");
     

     fprintf(FilePtr,"%lf\t",position[n][i].x);
     fprintf(FilePtr,"%lf\t",position[n][i].y);
     fprintf(FilePtr,"%lf\n",position[n][i].z);
  }

return;
}


void WriteVertice(FILE *FilePtr,int n)
{
 fprintf(FilePtr,"n = %d N = %d\n",n,NumberOfParticles);

 fprintf(FilePtr,"set vert(0) {%lf %lf %lf}\n",vertice[n][0].x,vertice[n][0].y,vertice[n][0].z);
 fprintf(FilePtr,"set vert(1) {%lf %lf %lf}\n",vertice[n][1].x,vertice[n][1].y,vertice[n][1].z);
 fprintf(FilePtr,"set vert(2) {%lf %lf %lf}\n",vertice[n][2].x,vertice[n][2].y,vertice[n][2].z);
 fprintf(FilePtr,"set vert(3) {%lf %lf %lf}\n",vertice[n][3].x,vertice[n][3].y,vertice[n][3].z);
 fprintf(FilePtr,"set vert(4) {%lf %lf %lf}\n",vertice[n][4].x,vertice[n][4].y,vertice[n][4].z);
 fprintf(FilePtr,"set vert(5) {%lf %lf %lf}\n",vertice[n][5].x,vertice[n][5].y,vertice[n][5].z);
 fprintf(FilePtr,"set vert(6) {%lf %lf %lf}\n",vertice[n][6].x,vertice[n][6].y,vertice[n][6].z);
 fprintf(FilePtr,"set vert(7) {%lf %lf %lf}\n",vertice[n][7].x,vertice[n][7].y,vertice[n][7].z);

 return;
}

void Writegene(FILE *FilePtr,int n)
{
 int i,d;
 fprintf(FilePtr,"%d %d\n",NumberOfParticles,Dimension);
 for(i=0;i<NumberOfParticles;i++)
 {
 fprintf(FilePtr,"%lf\t",position[n][i].x);
 fprintf(FilePtr,"%lf\t",position[n][i].y);
 fprintf(FilePtr,"%lf\n",position[n][i].z);
 }
 for(d=0;d<Dimension;d++)
 {
 fprintf(FilePtr,"%lf\t",projection[n][d].x);
 fprintf(FilePtr,"%lf\t",projection[n][d].y);
 fprintf(FilePtr,"%lf\n",projection[n][d].z);
 }
 //fprintf(FilePtr,"\n");

 return;
}
