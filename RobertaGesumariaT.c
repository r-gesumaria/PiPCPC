/*
 ============================================================================================
 Name        : RobertaGesumariaT.c
 Description : calcolo dell'approssimazione di Pi greco utilizzando la regola del trapeziode.
 ============================================================================================
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#define N 2000000
#define d 1E-7
#define d2 1E-14
#define num 2

int main(int argc, char* argv[]){
	double start, end;
	int  my_rank; /* rank of process */
	int  p;       /* number of processes */
	int source;   /* rank of sender */
	int mast = 0;     /* rank master */
	
	double result=0.0,x2=0.0; 	//variabili per il calcolo di pi greco
	double resultTot = 0.0; 	/*variabile del result totale ottenuto dalla somma rei result parziali
					dei vari processi*/
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD, &p); 
	
	int arrS[p][num], arrR[num];	//array send e receve utilizzati nella scatter
	start = MPI_Wtime();
	
	if(my_rank == 0){
		int quoz, resto, fine , disp, inizio=0;
	
		//divisione delle iterazioni per tutti i processi (compreso il master)
		quoz = N/p;
		resto =(int) N % p;
		disp= resto;

		for(source=0; source < p; source++){
			if(resto != 0){
				fine = (quoz + (resto!=0)) * (source+1);
				inizio = source * (quoz + (resto!=0));
				resto--;
			}else{
				fine = quoz*(source+1)+disp;
				inizio = source * quoz + disp;
			}
			arrS[source][0] = inizio;
			arrS[source][1] = fine;
		}
	}

	//invio delle variabili utili ai processi per compire il lavoro
	MPI_Scatter(arrS, num, MPI_INT, arrR, num, MPI_INT, mast, MPI_COMM_WORLD);
	
	//ogni processo acquisirà i valori utili per poter compiere il lavoro assegnatovi
	for (int i = arrR[0]; i < arrR[1]; i++){
		x2=d2*i*i;
		result+=1.0/(1.0+x2);
	}

	/*attraverso la reduce verranno sommati i result parziali ottenuti dai processi per poi inviarli al master che procederà con il
	 calcolo finale di pi greco*/
	MPI_Reduce(&result, &resultTot, 1, MPI_DOUBLE, MPI_SUM, mast,MPI_COMM_WORLD);

	//calcolo di Pi greco da parte del master
	if(my_rank == 0){
		double pi=0.0;	
		pi=4*d*resultTot;
		printf("PI = %lf \n", pi);
		end= MPI_Wtime();
		printf("*** tempo impiegato: %f\n",(end-start));
	}

	MPI_Finalize(); 
	return 0;
}
