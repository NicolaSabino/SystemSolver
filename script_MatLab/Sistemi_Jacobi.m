A = input ('Inserisci la matrice A dei coefficienti: ')
B = input ('Inserisci il vettore B dei termini noti: ')

D = zeros(3,3);
for i=1:3
    D(i,i) = A(i,i);
end;
D

L = [0, 0, 0; -A(2,1), 0, 0; -A(3,1), -A(3,2), 0]

U = [0, -A(1,2), -A(1,3); 0, 0, -A(2,3); 0, 0, 0]

k = input ('Inserisci il numero di iterazioni da eseguire ');

Vet = zeros(3,k);
Vet(1,1)=B(1)/A(1,1);
Vet(2,1)=B(2)/A(2,2);
Vet(3,1)=B(3)/A(3,3);
fprintf('Vettore di tentativo n. 1')
Vet(:,1)

for i = 1:k
    Vet(1,i+1) = (-A(1,2)*Vet(2,i)-A(1,3)*Vet(3,i)+B(1))/A(1,1);
    Vet(2,i+1) = (-A(2,1)*Vet(1,i)-A(2,3)*Vet(3,i)+B(2))/A(2,2);
    Vet(3,i+1) = (-A(3,1)*Vet(1,i)-A(3,2)*Vet(2,i)+B(3))/A(3,3);
    fprintf('Vettore di tentativo n. %1.0f',i+1)
    Vet(:,i+1)
    differenza = [Vet(1,i+1)-Vet(1,i); Vet(2,i+1)-Vet(2,i); Vet(3,i+1)-Vet(3,i)]
    ValoreMaxDiffernza=max(abs(differenza))
    CondConvergenza = (max(abs(differenza)))/max(abs(Vet(:,i+1)))
end

%differenza = [Vet(1,k)-Vet(1,k-1); Vet(2,k)-Vet(2,k-1); Vet(3,k)-Vet(3,k-1)]
%ErrRelativo = norm(differenza)/norm(Vet(:,k))
%max (abs(differenza))

Tj = (inv(D))*(L+U);
autovalori = eig(Tj)
lamda = zeros(3, 1);

for i=1:3
    lamda(i)=abs(autovalori(i));
end;

lamda
massimo = max(lamda)

if massimo>1
    fprintf('Il metodo di Jacobi non converge in quanto il massimo è maggiore di 1');
else
    fprintf('Il metodo di Jacobi converge in quanto il massimo è minore di 1')
end