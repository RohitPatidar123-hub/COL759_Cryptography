#include "mpi.h"
#include <stdio.h>


int main(int argc, char *argv[]) {
int  numtasks, rank, len, rc;
char hostname[MPI_MAX_PROCESSOR_NAME];

// initialize MPI
MPI_Init(&argc,&argv);

// get number of tasks
MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

// get my rank
MPI_Comm_rank(MPI_COMM_WORLD,&rank);
if(rank==0)
   {
    for(int i=1;i<numtasks;i++)   //send integer to each thread
       {
         MPI_Send(&i,1,MPI_INT,i,0,MPI_COMM_WORLD);
       }
    for(int i=1;i<numtasks;i++)   //send integer to each thread
       {
         int recv;
         MPI_Recv(&recv,1,MPI_INT,i,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
         //std::cout<<"Received "<<recv <<"from rank "<< i<<'\n';
         printf("Recived %d from rank %d\n",recv,i);
       }
   }else {
              int rec;
              MPI_Recv(&rec,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
              rec=rec*rec;
              MPI_Send(&rec,1,MPI_INT,0,0,MPI_COMM_WORLD);
         }

// this one is obvious
MPI_Get_processor_name(hostname, &len);
//printf ("Number of tasks= %d My rank= %d Running on %s\n", numtasks,rank,hostname);


   
MPI_Finalize();
}