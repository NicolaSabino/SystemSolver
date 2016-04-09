function result = gauss(M)
% GAUSS	Risolve un sistema lineare con il metodo di Guass.
%		Prende come parametro la matrice completa (coefficienti 
%		e termini noti), stampa tutti i passaggi e restituisce
%		il vettore delle soluzioni.
%		Esempio di utilizzo:
% x = GAUSS(M)
% See also GAUSS_PARTIAL_PIVOTING, GAUSS_SCALED_PIVOTING, LU_FACTOR.
	dim = size(M);
	% Controllo se la dimensione della matrice è del tipo nxn+1
	if(dim(2) - dim(1) == 1)
        fprintf('Benvenuto all esame di Analisi Numerica!\n');
        round(M*10^5)/10^5
		n = dim(1);
        for i = 1:n
            % Gauss
            for j = i:n-1
                m = M(j+1, i)/M(i,i);
                if (m ~= 0)
                    fprintf('m[%d, %d] = %f\n', j+1, i, m);
                    fprintf('[E(%d) - %f*E(%d)] -> [E(%d)]\n', j+1, m, i, j+1);
                    M(j+1,:) = M(j+1, :) - m*M(i,:);
                    round(M*10^5)/10^5
                end
            end
        end

		% Risoluzione all'indietro
		% Inizializzazione del vettore soluzione
		x = zeros(n, 1);
		% fliplr inverte l'ordine degli elementi di un vettore
		for i = fliplr(1:n)
			somma = 0;
            text = '';
			for j = i+1:n
                text = strcat(text, ' + (', num2str(M(i,j)), ')*(', num2str(x(j,1)), ') ');
                somma = somma + M(i,j)*x(j,1);
            end
            text = strcat('x[', num2str(i), '] = (', num2str(M(i, n+1)), ' - (', text, '))/(', num2str(M(i,i)), ')');
			x(i,1) = (M(i, n+1) - somma)/M(i,i);
            fprintf('%s\n', text);
		end

		% Stampo la soluzione
		fprintf('Soluzione \n');
		for i = 1:n
			fprintf('x[%d] = %f \n', i, x(i,1));
		end
		result = x;
		return;

	else
		fprintf('Dimensione non ammessa! Hai inserito anche il termine dei coefficienti?\n');
		return;
	end
end
