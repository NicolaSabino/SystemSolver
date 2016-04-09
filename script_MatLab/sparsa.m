 %Genera la matrice sparsa n*n a diagonale dominante
% usata come esempio in questo capitolo

function a=sparsa(n)

e=ones(n,1);
% per avere una matrice a diagonale dominante, diag>=5
diag=6;
b=[e, -e, diag*e, -e, 2*e];
d=[-n/2, -1, 0, 1, n/2];
a=spdiags(b,d,n,n);