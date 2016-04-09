#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

//variabili globali

    //Gauss
    int n, i, j, k, ind=0; //dimensione dell'input e contatori
	double alfa, rmax, temp, ck, s=0;
	//alfa=variabile di appoggio che contiene il valore di matrice[i][k]
	//rmax=elemento massimo della colonna
	//temp=variabile di appoggio
	//ck=coefficiente che moltiplica le righe per la diagonalizzazione
	//variabile per il calcolo delle soluzioni all'indietro

	//AcquisisciMatrice
	fstream y;//variabile di tipo fstream
	const int massimo=1003;
	double matrice[massimo][massimo];

	int SO,q;//serve per selezionare la sintassi per i comandi della shell



//  F U N Z I O N I--------------------------------

void pausa(){
    if(SO)system("pause");
    else{
        cout<<"\nQualsiasi lettera per continuare...\n";
        getchar();
        getchar();
    }

}


void stampa(){//ok
    cout<<endl;

    //stampa
        for(int i=1;i<=n;i++){//scorrimento delle righe
            for(int j=1;j<=n+1;j++){//scorrimento delle colonne
                cout<<matrice[i][j]<<"\t";
            }
            cout<<endl;
        }
        cout<<endl<<endl;

        //pausa
        pausa();
}

void acquisisciMatrice(){//ok
    char o[15];//creo una stringa di 15 caratteri

    cout<<"\nFile disponibili...\n\n";
    if(SO)system("dir");
    else system("ls");
    cout<<"________________________________"<<endl;
    cout<<"\nInserisci il nome del file contente la matrice da analizzare \n[specifica l'estensione .txt]\n";
    cin>>o;
    y.open(o,ios::in);//nel file di testo vanno inseriti in testa due caratteri
                            //indicanti le dimensioni della matrice

    //ACQUISIZIONE
    if(y.good()){
        //recupera il numero delle equazioni
        y>>n;

        for(int i=1;i<=n;i++){//scorrimento delle righe
            for(int j=1;j<=n+1;j++){//scorrimento delle colonne
                y>>matrice[i][j];//assegnamento
            }
        }


    }
    y.close();
    cout<<"Acquisizione avvenuta...";

    pausa();
    cout<<endl<<endl;

    return;
}


void Gauss(){//ok
    double soluzione[n];
    for(k=1;k<=n-1;k++){
	    rmax=abs(matrice[k][k]);
	    ind=k; //la variabile ind indica la riga sulla quale si sta operando
	    for(i=k+1;i<=n;i++){
	    	alfa=abs(matrice[i][k]);
	    	if(alfa>rmax){
	    		rmax=alfa;
	    		ind=i; //la variabile ind va ad indicare la riga che contiene il massimo
			}
		}

	    if(rmax==0){//ok
			cout<<"\n Il sistema non ammette soluzione\n\n";
            //pausa
            pausa();
		}

		if(ind!=k){ //Se ind e k avessero lo stesso valore si scambierebbe la riga con se stessa     ok

			for(j=k;j<=n+1;j++){
				temp=matrice[ind][i];
				matrice[ind][i]=matrice[k][j];
				matrice[k][j]=temp;
			}
		}

		//Diagonalizzazione della matrice dei coefficienti
		for(i=k+1;i<=n;i++){//ok
			ck=matrice[i][k]/matrice[k][k];
			for(j=k/*+1*/;j<=n+1;j++) matrice[i][j]-=ck*matrice[k][j];
		}

    }//chiusura for k

    //-------------------------------------------------
    //Calcolo del vettore delle soluzioni
	if(matrice[n][n]==0){//ok
		cout<<"\n Il sistema e' indeterminato\n\n";
        //pausa
        pausa();
	}
	soluzione[n]=matrice[n][n+1]/matrice[n][n];
	for(i=n-1;i>=1;i--){//ok
		for(j=i+1;j<=n;j++){
			s+=matrice[i][j]*soluzione[j];
		}
		soluzione[i]=(matrice[i][n+1]-s)/matrice[i][i];
	}
    //-------------------------------------------------
	//Stampa del vettore soluzione
	cout<<"\n La soluzione e':\n\n";
	for(i=1;i<=n;i++) cout<<soluzione[i]<<endl<<endl;
    //-------------------------------------------------

    //pausa
    pausa();
}

void Jacobi(){//ok
    cout<<"\nJacobi\n\n";
    //variabili
    int  i=1, j=1;
    double eps, Rmax, dif;

    for (i=1; i<=n; i++){
            double s = 0;
            for (j=1; j<=n; j++){
                s += abs(matrice[i][j]);
            }
            s -= abs(matrice[i][i]);
            if (abs(matrice[i][i])<s){
                cout <<"\nLa matrice non e' a diagonale dominante.\n\n";
                pausa();
                return ;
            }
        }

    cout << "\nInserire il coefficiente di tolleranza\n\n";
    cin >> eps;

    double tentativo[n], tentativo2[n];
    for(i=1;i<=n;i++){
        tentativo[i]=matrice[i][n+1]/matrice[i][i];
    }
    for (i=1; i<=n; i++){
        float s=0;
        for (j=1; j<=n; j++){
            s+=matrice[i][j]*tentativo[j];
        }
        s-=matrice[i][i]*tentativo[i];
        tentativo2[i]=(matrice[i][n+1]-s)/matrice[i][i];
        dif=abs(tentativo2[i]-tentativo[i]);
        if (dif>Rmax) Rmax=dif;
    }
    if(dif>eps){
        cout << "\n Il sistema non ammette una soluzione accettabile\n\n";
        //pausa
        pausa();
        return ;
    }
    else{
        cout << "\n La soluzione e': \n\n";
        for (i=1; i<=n; i++){
            tentativo[i]=tentativo2[i];
            cout << tentativo[i] << endl;
        }
    }
    cout<<endl<<endl;

    //pausa
    pausa();
}

void GSeidel(){
    //Controllo diagonale dominante
    for(i=1;i<=n;i++)
	{
		float s=0;
		j=1;
		for(j=1;j<=n;j++)
		{
			s+=abs(matrice[i][j]);
		}
		s-=abs(matrice[i][i]);
		if(abs(matrice[i][i])<s)
		{
			cout << " La matrice non soddisfa la condizione di dominanza della diagonale\n\n";
			//pausa
            pausa();
		}
	}

	//Creazione vettore tentativo
	double tentativo[n], tentativo1[n];
	for(i=1;i<=n;i++) tentativo[i]=(matrice[i][n+1])/(matrice[i][i]);

	//Calcolo del vettore soluzione
	double eps, diff, rmax=0;
    cout << " Inserire il fattore di tolleranza per la soluzione\n\n";
    cin >> eps;
    cout << endl << endl;
    for(i=1;i<=n;i++){
    	double s=0;
    	for(j=1;j<=n;j++){
    		s+=matrice[i][j]*tentativo[j];
		}
		s-=matrice[i][i]*tentativo[i];
		tentativo1[i]=(matrice[i][n+1]-s)/matrice[i][i];
		diff=abs(tentativo1[i]-tentativo[i]);
		if(diff>rmax) rmax=diff;
	}
	//Stampa del vettore soluzione
    if(diff>eps){
			cout << " Il sistema non può essere risolto con un'approssimazione valida\n\n";
			//pausa
            pausa();
		}
    else{
    	cout << " Il vettore soluzione e':\n\n";
		for(i=1;i<=n;i++){
    		tentativo[i]=tentativo1[i];
    		cout << tentativo[i] << endl;
		}
		cout << endl;
	}

	//pausa
	pausa();

    return ;

}

void edit(){//richiama un editor di testo per inserire una matrice
    if(SO)system("edit NuovaMatrice.txt");
    else system("nano NuovaMatrice.txt");
    return;
}
//  M A I N ---------------------------------------
int main(){


    char c;
    //------------------------------------
    int h=1;
    cout<<"Quale Sistema Operativo stai utilizzando?\n0 per Unix-Linux-OSX\n1 per Microsoft Windows\n\n";
    cin>>SO;
    //------------------------------------
    while(1) {

    if(SO)system("cls");
    else system("clear");

    cout << "\tMetodi risolutivi per sistemi lineari\n\n";


    cout<<"a:acquisici da file\ns:stampa\ng:Gauss\nk:Gauss-Seidel\nj:Jacobi\nl:matrice da editor di testo\ne:esci\n\n";
    cin>>c;
    switch(c){
        case 'a': acquisisciMatrice();break;
        //case 'i': inserisci();break;
        case 's': stampa();break;
        case 'g': Gauss();break;
        case 'j': Jacobi();break;
        case 'k': GSeidel();break;
        case 'l': edit();break;
        case 'e': return 0;
        default:cout<<"\tCARATTERE NON VALIDO!\n";
    }//chiusura dello switch

    }//chiusura del while
}
