/*
 ====================================================================================================
 Name        : RobertaGesumariaMC.c
 Description : calcolo dell'approssimazione del velore di Pi greco utilizzando il metodo Monte Carlo.
 ====================================================================================================
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define SEED 35791246
#define NUM 2

int main(int argc, char* argv[]){
	
	double start, end;
	int  my_rank; 	/* rank of process */
	int  p;       	/* number of processes */
	int mast = 0; 	/* rank master */
	int source;   	/* rank of sender */

	int arrR[NUM];		//array receve per la scatter
	int niter;		//numeri di iterazioni totali
	int totCount = 0;	//variabile  per la somma dei count ricevuti dai vari processi
	int count=0;		//variabili per il calcolo di pi greco
	double x,y,z,pi; 	//variabili per il calcolo di pi greco
	
	if(argc !=2){
		printf("Utilizzo: %s <numero iterazioni totali>\n", argv[0]);
		return -1;
	}
		
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	int arrS[p][NUM];	//array send per la scatter
	
	
	
	niter= atoi(argv[1]);
	
	/*creazione di un array di numeri random dal quale i processi acquisiranno i numeri per il calcolo di pi greco*/
	double *numCas;
	numCas = (double*) malloc(sizeof(double)*(NUM*niter));
	srand(SEED);
	for(int i=0; i<(NUM*niter); i++){
		numCas[i]=(double)rand()/RAND_MAX;
	}

	start = MPI_Wtime();

	if(my_rank == 0){
		int resto, quoz, inizio, disp, itps=0;
		
		/*divisione delle iterazioni per processo (compreso il master)*/
		quoz = niter/p;
		resto= niter%p;
		disp = resto * NUM;
		
		for(source=0; source < p; source++){
			if(resto != 0){
				itps = quoz + (resto!=0);
				inizio = source * itps * NUM;
				resto--;
			}else{
				itps = quoz;
				inizio = source * itps * NUM + disp;
			}
			arrS[source][0] = itps;
			arrS[source][1] = inizio;
		}
	}
	
	//invio delle iterazioni da compiere per ogni processo
	MPI_Scatter(arrS, NUM, MPI_INT, arrR, NUM, MPI_INT, mast, MPI_COMM_WORLD);
	
	//ogni processo acquisirà i valori per poter compiere il lavoro assegnato
	int it = arrR[0];
	int inc = arrR[1];
	for(int i=0; i<it; i++) {
		x = numCas[inc++];
		y = numCas[inc++];
		z = x*x+y*y;
		if(z<=1) count++;
	}
	
	/*attraverso la reduce verranno sommati i vari count ottenuti dai processiper poi inviarli al master che procederà con il calcolo
	finale di pi greco*/
	MPI_Reduce(&count, &totCount, 1, MPI_INT, MPI_SUM, mast,MPI_COMM_WORLD);
	
	if(my_rank == 0){
	//calcolo di pi greco da parte del master
		pi=(double)totCount/niter*4;
		printf("# of trials= %d , estimate of pi is %f \n",niter,pi);
		
		end = MPI_Wtime();
		printf("*** tempo impiegato %f***\n",(end-start));
	}

	MPI_Finalize();
	return 0;
}


