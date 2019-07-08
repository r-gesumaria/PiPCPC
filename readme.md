# Pi greco
## 
I due programmi che verranno documentati sono stati sviluppati in C utilizzando il protocollo MPI e calcolano in parallelo un'approssimazione del valore di π.
Avremo un primo programma che utilizza il metodo di Monte Carlo e un secondo che utilizza la regola del trapezoide per ottenere l'approssimazione del valore π.

## Esecuzione

**Prima di eseguire il programma**
Posizionarsi nella directory contenente il pacchetto.tar.gz del progetto ed estrarlo attraverso il comando:

	tar -xzvf RobertaGesumariaPiGreco.tar.gz
		
Una volta estratto ci si deve spostare all'interno della directory root del progetto attraverso il comando:

	cd RobertaGesumariaPiGreco

**Per eseguire il programma:**
*Nomi dei due programmi:*
- *RobertaGesumariaMC* - programma che utilizza il metodo Monte Carlo;
- *RobertaGesumariaT* - programma che utilizza la regola del trapezoide.

Compilare il programma .c con il seguente comando:

	mpicc nomeProgramma.c -o nomeProgramma

**Esecuzione dei programmi**
Eseguire il programma che utilizza il metodo Monte Carlo attraverso il comando:


    mpirun -np numProcessori RobertaGesumariaMC numeroIterazioniTotali

mentre per eseguire il programma che utilizza la regola del trapezoide usare:

	mpirun -np numProcessori RobertaGesumariaT

## METODO MONTE CARLO
Il metodo "Monte Carlo" è una strategia di risoluzione di problemi che utilizza la statistica: se la probabilità di un certo evento è P possiamo simulare in maniera random questo evento e ottenere P facendo (numero di volte in cui il nostro evento è avvenuto)/(simulazioni totali).
Per applicare questa strategia al calcolo del π, dato un cerchio di raggio 1, esso può essere inscritto in un quadrato di raggio 2. Prendiamo in considerazione solo 1/4 del quadrato. Così facendo l'area del quadrato sarà 1 e l'area dello spicchio di cerchio sarà di π /4.
Se generiamo N numeri random all'interno del quadrato il numero di punti che cadono nel cerchio M diviso il numero totale di numeri generati N dovrà approssimare appunto l'area del cerchio e quindi π /4.
In sostanza otterremo che π = 4 * M / N.

### Input
Per poter eseguire il programma sarà necessario indicare il numero di iterazioni che si vorranno svolgere per poter approssimare al meglio il valore di π, ovviamente maggiori saranno il numero di iterazioni, migliore sarà l'approssimazione.

### Scelte progettuali
Per approcciare il problema sono state effettuate le seguenti scelte progettuali:

* ##### Matrice dei numeri casuali
Come suddetto c'è bisogno di un insieme di numeri casuali.
Per poter far sì che ogni processo abbia un set di numeri differenti da usare, si è scelto di creare un array di numeri random che sia il doppio (NUM=2) del numero di iterazioni date in input poiché ad ogni processo, per ogni ciclo, occorreranno due numeri casuali.
```
	double *numCas;
	numCas = (double*) malloc(sizeof(double)*(NUM*niter));
	srand(SEED);
	for(int i=0; i<(NUM*niter); i++){
		numCas[i]=(double)rand()/RAND_MAX;
	}
```

* ##### Divisione ed invio del carico di lavoro
La divisione del carico di lavoro viene distribuita tra tutti i processi, compreso il processo master. Nel caso in cui sia presente resto nella divisione del carico, questo sara distribuito tra i processi.
Come carico di lavoro si intende il numero di iterazioni che ogni processo dovrà svolgere (*itp*) ed il punto di inizio da cui dovrà leggere i valori all'interno dell'array di numeri casuali (*inizio*).
Queste informazioni verranno memorizzate nell'array di invio (*arrS*) utilizzato dalla routine di comunicazione collettiva *MPI_Scatter* utilizzata per inviare il lavoro ad ogni singolo processo.
Il carico di lavoro viene suddiviso dal master.
```
		int resto, quoz, inizio, disp, itps=0;
		quoz = niter/p;
		resto= niter%p;
		disp = resto * NUM;

		for(source=0; source < p; source++){
			if(resto != 0){
				itps = quoz + (resto!=0);
				inizio = source * itp * NUM;
				resto--;
			}else{
				itps = quoz;
				inizio = source * itp * NUM + disp;
			}
			arrS[source][0] = itp;
			arrS[source][1] = inizio;
		}

	MPI_Scatter(arrS, NUM, MPI_INT, arrR, NUM, MPI_INT, mast,,MPI_COMM_WORLD);
```

* ##### Elaborazione e restituzione del lavoro
Ogni processo riceverà il carico di lavoro che gli spetta attraverso l'array di recezione della *Scatter* (*arrR*).
Dopo aver eseguito il calcolo del valore z e aver valutato se incrementare il contatore (*count*), che indicherà la quantità di numeri che sono all'interno dello spicchio di cerchio preso in considerazione, si procederà alla sommatoria dei singoli contatori risultati dall'esecuzione dei vari processi e all'invio del totale al processo master, attraverso la routine *MPI_Reduce*.
```
	int it = arrR[0];
	int inc = arrR[1];
	for(int i=0; i<it; i++) {
		x = numCas[inc++];
		y = numCas[inc++];
		z = x*x+y*y;
		if(z<=1) count++;
	}
	MPI_Reduce(&count, &totCount, 1, MPI_INT, MPI_SUM, mast,MPI_COMM_WORLD);
```

* ##### Elaborazione finale del processo master
In seguito alla ricezione del risultato totale (*totCount*), risultatante dall'elaborazione dei contatori locali di tutti i processi, il master procede con il calcolo dell'approssimazione di π (*pi*).
```
	pi=(double)totCount/niter*4;
	printf("# of trials= %d , estimate of pi is %f \n",niter,pi);
```

## REGOLA DEL TRAPEZOIDE
Il valore di π è calcolato attraverso una formula di quadratura trapezoidale. E' possibile approssimare π/4 con la somma Sn delle aree dei trapezi di altezza 1/n inscritti in un quarto di circonferenza di raggio unitario.
Infatti, la funzione che definisce un arco di circonferenza unitaria nel I quadrante è:
**f(x) = √(1-x^2)**
Per cui integrando numericamente f(x) in [0,1] si ottiene che:
 **∫√(1-x^2)dx ≃ π/4**

### Input
Il problema non ha bisogno di alcun input esterno. L'approssimazione del π verrà effettuata su un numero di iterazioni fisse pari a 1E7.

### Scelte progettuali
Per approcciare il problema sono state effettuate le seguenti scelte progettuali:

* ##### Divisione del carico di lavoro
La divisione del carico di lavoro viene distribuita tra tutti i processi, compreso il processo master. Nel caso in cui sia presente resto nella divisione del carico, questo sara distribuito tra i processi.
Come carico di lavoro si intende il punto di inizio e fine di ogni ciclo di iterazioni che ogni processo dovrà eseguire. Le seguenti informazioni sono memorizzate nell'array di invio (*arrS*) utilizzato dalla routine di comunicazione collettiva *MPI_Scatter* che provvederà all'invio delle informazioni ai vari processi.
Il carico di lavoro viene suddiviso dal master.
```
		int quoz, resto, fine, disp, inizio=0;
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
	MPI_Scatter(arrS, num, MPI_INT, arrR, num, MPI_INT, mast, MPI_COMM_WORLD);
```

* ##### Elaborazione e restituzione del lavoro
Ogni processo riceverà il carico di lavoro che gli spetta attraverso l'array di recezione della *Scatter* (*arrR*).
Dopo aver eseguito la parte di lavoro dei singoli processi, si procederà alla sommatoria dei singoli risultati parziali ottenuti dall'esecuzione dei vari processi e all'invio del totale al processo master, attraverso la routine *MPI_Reduce*.
```
	for (int i = arrR[0]; i < arrR[1]; i++){
		x2=d2*i*i;
		result+=1.0/(1.0+x2);
	}
	MPI_Reduce(&result, &resultTot, 1, MPI_DOUBLE, MPI_SUM, mast, MPI_COMM_WORLD);
```

* ##### Elaborazione finale del processo master
In seguito alla ricezione dei risultati dell'elaborazioni dei processi (resultTot), il master procede con il calcolo dell'approssimazione di π (*pi*).
```
		double pi=0.0;
		pi=4*d*resultTot;
		printf("PI = %lf \n", pi);
```

## MPI_Wtime()
In entrambi i progetti il tempo di esecuzione è stato calcolato utilizzando la funzione *MPI_Wtime()* per ottenere il tempo di inizio (*start*) e fine(*end*).

##Benchmark
Verranno mostrati i risultati dei test di entrambi i programmi effettuati utilizzando delle istanze di tipo m4.large (2 core) di Amazon Web Services con StarCluster-Ubuntu_12.04-x86_64-hvm - ami-52a0c53b.

### Strong scaling
Per testare la strong scaling sono stati effettuati dei test in cui si è tenuta fissa la dimensione del problema mentre si è aumentato il numero di elementi di elaborazione.
L'obiettivo di questo test è quello di trovare un "sweet spot" che consenta il completamento del calcolo in un lasso di tempo ragionevole e che non sprechi troppi cicli a causa del sovraccarico parallelo.

**Calcolo**
Se la quantità di tempo per completare un'unità di lavoro con 1 elemento di elaborazione è indicata con *t1* e la quantità di tempo per completare la stessa unità di lavoro con gli elementi di elaborazione N è *tN*, la strong scalability è data come:

**t1/(N x tN) x 100%**

###Weak scaling
Per testare la weak scaling sono stati effettuati dei test in cui la dimensione del problema (carico di lavoro) assegnata a ciascun elemento di elaborazione è rimasta costante e sono stati utilizzati elementi aggiuntivi per risolvere un problema totale più grande.
In caso di weak scalability, il ridimensionamento lineare viene ottenuto se il tempo di esecuzione rimane costante mentre il carico di lavoro viene aumentato in proporzione diretta al numero di processori.

**Calcolo**
Se la quantità di tempo per completare un'unità di lavoro con 1 elemento di elaborazione è *t1* e la quantità di tempo per completare N delle stesse unità di lavoro con gli elementi di elaborazione N è *tN*,la weak scalability è data come:

**(t1  x tN) x 100%**

## CONFRONTO DEI DUE METODI
**Valore di π = 3,14159 26.**
#### Strong scalability

* ##### Metodo Monte Carlo
Approssimazione di π = 3.141580
|Core MPI| Tempo (s) | Percentuale |
|--------|-----------|-------------|
|1		 | 2.477 	 |	-		   |
|2       | 1.258 	 |	98.44%     |
|4       | 0.625 	 |	99.11%	   |
|6       | 0.434 	 |	95.14%	   |
|8       | 0.323  	 |	95.75%	   |
|10      | 0.366 	 |	67.64%	   |
|12      | 0.351 	 |	58.83%	   |
|14      | 0.343 	 |	51.59%	   |
|16      | 0.270	 |	57.27%	   |
![](https://github.com/Edilio1995/ASimpleJacobiPcPc/blob/master/MontecarloStrong.JPG?raw=true)

Dalla tabella, risultante dai test effettuati, è possibile notare come lo sweet spot ideale sarebbe a 4 core, ma si ha un buon speedup fino a 16 core nonostante l'overhead di comunicazione tra i processi.

* ##### Regola del trapezoide
Approssimazione di π = 3.141593
|Core MPI| Tempo (s) | Percentuale |
|--------|-----------|-------------|
|1		 | 0.130 	 |	-		   |
|2       | 0.066 	 |	98.12%     |
|4       | 0.035 	 |	93.61%	   |
|6       | 0.025 	 |	87.14%	   |
|8       | 0.023  	 |	70.77%	   |
|10      | 0.021 	 |	60.92%	   |
|12      | 0.015 	 |	70.41%	   |
|14      | 0.018 	 |	52.60%	   |
|16      | 0.018	 |	42.70%	   |
![](https://github.com/Edilio1995/ASimpleJacobiPcPc/blob/master/trapezoideStrong.JPG?raw=true)

Dai test effettuati è possibile notare come lo sweet spot ideale sarebbe a 2 processori, ma si continuano ad avere ottime prestazioni fino a 14 core. Con 16 core si può notare un calo delle prestazioni dovute all'overhead di comunicazione tra i processi.

##### Confronto
1. **Approssimazione.**
	- L'approssimazione ricavata con il metodo di Monte Carlo (3.141580) si avvicina al valore di pi greco (3,14159 26) ma risulta essere più corretta quella ottenuta con la regola del trapezoide (3.141593).
2. **Tempo di esecuzione.**
	- I tempi di esecuzione ottenuti con la regola del trapezoide sono nettamente inferiori a quelli ottenuti dall'esecuzione del metodo Monte Carlo.
![](https://github.com/Edilio1995/ASimpleJacobiPcPc/blob/master/ConfrontoTempiEs.JPG?raw=true)

3. **Scalabilità (SpeedUp)**
	- Il programma che utilizza il metodo Monte Carlo scala piuttosto bene all'aumentare dei core mantenendo un buon speedup nonostante l'overhead di comunicazione. Altrettando vale per il programma che utilizza la regola del trapezoide anche se con l'utilizzo di 16 core ha un peggioramento delle prestazioni, dovuto appunto all'overhead di comunicazione necessario tra i processori.

#### Weak scalability

* ##### Metodo Monte Carlo
|Num. iterazioni| Core MPI | Tempo (s)| Percentuale |    π     |
|---------------|----------|----------|-------------|----------|
|15.000.000    	|1         |   0,155  |	-	        | 3.141436 |
|30.000.000	   	|2         |   0,157  |	98,69%      | 3.141492 |
|60.000.000	   	|4         |   0,159  |	98,06%	    | 3.141670 |
|90.000.000	   	|6         |   0,160  |	97.34%	    | 3.141655 |
|120.000.000	|8         |   0,162  |	96.02%	    | 3.141657 |
|150.000.000	|10        |   0,233  |	66.57%	    | 3.141697 |
|180.000.000	|12        |   0,243  |	64.02%	    | 3.141581 |
|210.000.000	|14        |   0,277  |	56.21%	    | 3.141629 |
|240.000.000   	|16        |   0,275  |	56.53%	    | 3.141580 |
![](https://github.com/Edilio1995/ASimpleJacobiPcPc/blob/master/MontecarloWeak.JPG?raw=true)
Dai test effettuati è possibile notare come il programma si comporta bene fine ad 8 core, per avere un peggioramento delle prestazioni con più di 10 processori, dovuto all'l'overhead di comunicazione necessario tra i processori, così come per la strong scaling.

* ##### Regola del trapezoide
Non potendo variare il numero di iterazioni da effettuare per ottenere una giusta approssimazione di π, non sarà possibile effettuare questo tipo di test sul programma che utilizza la regola del trapezoide.
