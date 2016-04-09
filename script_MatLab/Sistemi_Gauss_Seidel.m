A = input ('Inserisci la matrice A dei coefficienti: ')
b = input ('Inserisci il vettore B dei termini noti: ')
[n,m] = size(b);

D = zeros(n,n);
for i=1:3
    D(i,i) = A(i,i);
end;
D


L = [0, 0, 0; -A(2,1), 0, 0; -A(3,1), -A(3,2), 0]

U = [0, -A(1,2), -A(1,3); 0, 0, -A(2,3); 0, 0, 0]

%x(k) = inv(D-L)*U*x(k-1)+(D-L)*B
x = input('Inserisci il vettore di tentativo iniziale: ');

k = input('Inserisci il numero di iterazioni da eseguire: ');
Vet = zeros(3,k);
Vet(1,1) = x(1);
Vet(2,1) = x(2);
Vet(3,1) = x(3);
%fprintf ('Vettore di tentativo n. 0');
%Vet(:,1)

for i=1:3 %scandisce indice i (=righe)
    Vet(1,i+1)= (b(1)-A(1,2)*Vet(2,i)-A(1,3)*Vet(3,i))/A(1,1);
    Vet(2,i+1)= (b(2)-A(2,1)*Vet(1,i+1)-A(2,3)*Vet(3,i))/A(2,2);
    Vet(3,i+1)= (b(3)-A(3,1)*Vet(1,i+1)-A(3,2)*Vet(2,i+1))/A(3,3);
    Vet(:,i+1)
    differenza = [Vet(1,i+1)-Vet(1,i); Vet(2,i+1)-Vet(2,i); Vet(3,i+1)-Vet(3,i)]
    ValoreMaxDiffernza = max(abs(differenza))
    CondConvergenza = (max(abs(differenza)))/(max(Vet(:,i+1)))
end

Tg=inv(D-L)*U;
autovalori = eig (Tg)
AutovaloriValoreAssoluto = zeros(3,1);

for i = 1:3
    AutovaloriValoreAssoluto(i) = abs(autovalori(i));
end;

AutovaloriValoreAssoluto
massimoAutovalori = max(abs(AutovaloriValoreAssoluto))

if massimoAutovalori>1
    fprintf('Il metodo di Gauss Siedel non converge in quanto il massimo è maggiore di 1');
else
    fprintf('Il metodo di Gauss Siedel converge in quanto il massimo è minore di 1')
end
