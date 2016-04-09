%-------GENERA UNA MATRICE A DIAGONALE DOMINANTE DELL'ORDINE SCELTO--------
%
% IMPORTANTE: dopo aver generato la matrice controllare che sia
% effettivamente a diagonale dominante usando la funzione domdiag.m, che
% restituisce 1 in caso affermativo e 0 in caso negativo

fprintf('\nLo script genera una matrice a diagonale dominante della dimensione scelta\n');
n=input('Inserire la dimensione (>0): ');
%Controllo che la dimensione inserita sia valida (>0)
while(n<=0)
    fprintf('dimensione non valida!\n');
    n=input('Inserire la dimensione (>0): ');
end

C=rand(n);  %genera una matrice con coefficienti casuali compresi tra 0 e 1
B=C;        %variabile di appoggio
B(eye(size(C))~=0)=0; %???
z=max([sum(B,1)' sum(B,2)],[],2); %???
B(eye(size(B))~=0)=z; %???
D=eye(n);   %genera la matrice identità
B=B+D*0.06; %aumento un po' il valore degli elementi sulla diagonale in 
            %modo da rendere la relazione di dominanza strettamente vera
B=B*100;    %gli elementi generati sono minori di 1, così li rendo maggiori
A=floor(B)  %arrotondo i valori della matrice considerandone solo la parte intera
clear B; %(
clear z; %
clear C; %   pulizia generale
clear D; %
clear n; %                     )