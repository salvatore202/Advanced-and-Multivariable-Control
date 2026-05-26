clear all
close all
clc

%% SISTEMA con incertezza moltiplicativa (SISO)

s = tf('s');
% Sistema nominale
%G = 10/(s+1);
G = 10/(s+1)^2;
%G = 10/(s+1)^3;

% incertezza 
Delta = 1/(s+5);



% CASO MIMO

% Sistema nominale
 G = [5/(s+1)^2   2/(s+2);
      1/((s+3)*(s+1))   4/(s+1)];

% Incertezza moltiplicativa
 Delta = [1/(s+5) 0;
          0         2/(s+8)];



%% verifica applicabilità del teorema del piccolo guadagno

% Funzione a ciclo chiuso nominale (complementary sensitivity)
T0 = feedback(G,eye(length(G)));    % G/(1+G)

% verifica della stabilità di T0 e di Delta
polesT0=eig(T0);
polesDelta=eig(Delta);
if max(real([polesT0 ; polesDelta]))>0
    disp('Teorema del piccolo guadagno non applicabile perchè i due sistemi ad anello aperto non sono asint. stabili')
    return
end


%% Applicazione del teorema del piccolo guadagno

% Trova il valore di k limite con il teorema del piccolo guadagno
passo=0.01;
k_cond=0;  % ci serve per far continuare il ciclo
k=passo;
normT0 = norm(T0, Inf);  % norma infinito di T0
while k_cond==0
    DeltaK = k*Delta;
    normDelta = norm(DeltaK, Inf); 
    cond = normT0 * normDelta; % di solito più conservativa
    cond2 = norm(T0*DeltaK, Inf); 
    if cond2<1
        k=k+passo;
    else
        k_cond=1;
    end
end
fprintf('Il massimo valore di k che verifica il teorema del piccolo guadagno è  k= %+5.2f  ||T0||*||Delta|| = %.6f ||T0*Delta|| = %.6f\n', ...
        k, cond,cond2);

%% Test risposta per k maggiore e minore del limite

k_values=[k/2 k-passo 2*k 10*k ];

legendEntries = {};
for k = k_values
    DeltaK = k*Delta;
    cond2 = norm(DeltaK*T0, Inf);
    fprintf('k = %+5.2f  ||T0*Delta|| = %.6f\n', k, cond2);
    % Sistema reale e chiusura del loop
    G_real = G*(eye(length(G)) + DeltaK);       % plant reale con incertezza moltiplicativa
    T_real = feedback(G_real, eye(length(G_real))); % funzione a ciclo chiuso reale
    % Step response
    figure
    subplot(1,2,1)
    step(T_real);
    legendEntries = sprintf('k=%+.2f (norma=%.2f)', k, cond2);
    legend(legendEntries, 'Location', 'Best','Fontsize',14);
    title('Step responses','FontSize',14);
    xlabel('Time [s]'); ylabel('y(t)');
       % Nyquist diagrams of L(s)
    subplot(1,2,2)
    nyquist(G_real);
    legendEntries = sprintf('k=%+.2f (norma=%.2f)', k, cond2);
    legend(legendEntries, 'Location', 'Best','Fontsize',14);
    title('Nyquist diagrams','FontSize',14);
end


