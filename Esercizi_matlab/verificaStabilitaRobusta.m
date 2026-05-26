%%------DEFINIZIONE--------

s = tf('s');

% Nominal system
G_0 = 10/(s+1);

% Multiplicative uncertainty
deltaG_m = 1/(s+5);  % keep in mind that it's k/(s+5)



%%------SVOLGIMENTO--------

% Nominal closed-loop function
T_0 = feedback(G_0, 1);  % the 'feedback' command makes T=G/(1+G)

% Controlliamo se la closed-loop function e la deltaG_m sono asintoticamente stabili
poles_T_0 = eig(T_0);
poles_deltaG_m = eig(deltaG_m);
if max(real([poles_T_0 ; poles_deltaG_m])) > 0
    disp('Small gain Theorem is not applicable')
    return
else
    disp('Small gain Theorem is applicable')
end


%calcoliamo i valori di k che fanno si che il sistema rimanga robustamente
%stabili

passo = 0.01; % Questa è la risoluzione. Più è piccolo, più il k sarà preciso.
k_cond = 0;   % Una bandierina (flag) per fermare il ciclo.
k = passo;    % Partiamo da un k piccolissimo.

while k_cond == 0
    DeltaK = k * deltaG_m;             % Moltiplica l'incertezza per il k attuale.
    cond2 = norm(T_0 * DeltaK, Inf); % Calcola il picco massimo del modulo del prodotto.
    
    if cond2 < 1
        k = k + passo;              % Se siamo sotto 1, il sistema è stabile. Aumentiamo k!
    else
        k_cond = 1;                 % Se cond2 >= 1, abbiamo trovato il limite. Stop!
    end
end

fprintf('Il valore di k può variare tra %f e %f\n', -k, k);