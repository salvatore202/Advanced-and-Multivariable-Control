% Pole-placement stato (MIMO), osservatore di Luenberger,
% controllo tramite stima (principio di separazione).
% Sistema d'esempio: ordine 3, 2 ingressi, 2 uscite.

clear; close all; clc;
%%
% -------------------------
% 1) Definizione del sistema
% -------------------------
A = [0 1 0;
     0 0 1;
    -2 -3 -4];        % 3x3
B = [0 0;
     0 1;
     1 0];           % 3x2
C = [1 0 0;
     0 0 1];         % 2x3
D = zeros(2,2);

%%
% Verifica controllabilità e osservabilità
Ctr = [B, A*B, A^2*B];       % controllability matrix (3x6), ctrb(A,B)
Ob  = [C; C*A; C*A^2];       % observability matrix (6x3), obsv(A,C)
rankCtr = rank(Ctr); 
rankOb  = rank(Ob);
fprintf('Rank controllability: %d (deve essere 3)\n', rankCtr);
fprintf('Rank observability:   %d (deve essere 3)\n', rankOb);

% -------------------------
% 2) Progetto retroazione di stato (K) - pole placement
% -------------------------
% Autovalori desiderati per il controllore (scelti a sinistra nel piano reale)
poles_ctrl = [-2, -3, -4];     % vettore di 3 poli reali, puoi cambiare
% Nota: usa place per sistemi MIMO
K = -place(A, B, poles_ctrl);   % K è 2x3 (per u = -K*x)
fprintf('K (state feedback) =\n'); disp(K);

% Autovalori effettivi del sistema in retroazione di stato
eig_Acl = eig(A + B*K);
fprintf('Autovalori A+BK (controller closed-loop):\n'); disp(eig_Acl);

% -------------------------
% 3) Progetto osservatore di Luenberger (L)
% -------------------------
% Scegli poli dell'osservatore più a sinistra per avere stima rapida
poles_obs = [-40, -50, -60];   % più veloci rispetto al controllore
% Osservatore: progetto L come place(A', C', poles_obs)' (dual)
L = place(A', C', poles_obs)';  % L è 3x2
fprintf('L (observer gain) =\n'); disp(L);

% Autovalori dell'errore di stima
eig_Aobs = eig(A - L*C);
fprintf('Autovalori A-LC (observer):\n'); disp(eig_Aobs);

% -------------------------
% 4) Principio di separazione: verifica autovalori complessivi
% -------------------------
% Matrice di stato aumentata nello spazio [x; e] (con e = x - xhat)
A11 = A + B*K;
A12 = B*K;
A21 = zeros(size(A));
A22 = A - L*C;
A_aug = [A11 A12; A21 A22];
eig_aug = eig(A_aug);
n_aug=length(eig_aug);

fprintf('Autovalori del sistema separato (unione controller+observer):\n');
disp(eig_aug);

% Controlla che gli autovalori dell'aug siano l'unione di eig(A11) e eig(A22)
fprintf('Autovalori controller: '); disp(sort(real(eig(A11))));
fprintf('Autovalori observer:   '); disp(sort(real(eig(A22))));
fprintf('Autovalori aug:        '); disp(sort(real(eig_aug)));

% -------------------------
% 5) Simulazione closed-loop con osservatore (u = -K*xhat)
% -------------------------
% Dinamica implementata:
% x_dot      = A*x + B*u
% xhat_dot   = A*xhat + B*u + L*(y - C*xhat)
% u = -K*xhat
%
% Sistema per integrazione: z = [x; xhat] in R^{6}

% Parametri simulazione
tspan = 0:0.01:8;    % da 0 a 8 s
x0 = [0.5; -0.8; 0.3];           % stato iniziale
xhat0 = [0; 0; 0];               % stima iniziale (inizialmente nulla)
z0 = [x0; x0-xhat0];

sys=ss(A_aug,zeros(n_aug,1),eye(n_aug),zeros(n_aug,1));
[y,tt,zzz] = initial(sys,z0);
x = zzz(:,1:3)';
e = zzz(:,4:6)';    
xhat = x-e;
y = C * x;
yhat = C * xhat;


% -------------------------
% 6) Grafici
% -------------------------
figure('Name','Stati reali vs stimati');
for i=1:3
    subplot(3,1,i);
    plot(tt, x(i,:), 'LineWidth', 1.5); hold on;
    plot(tt, xhat(i,:), '--', 'LineWidth', 1.3);
    xlabel('Tempo (s)'); ylabel(sprintf('x_%d',i));
    legend('x','x_{hat}'); grid on;
end

figure('Name','Uscite reali vs stimate');
for i=1:2
    subplot(2,1,i);
    plot(tt, y(i,:), 'LineWidth', 1.5); hold on;
    plot(tt, yhat(i,:), '--', 'LineWidth', 1.3);
    xlabel('Tempo (s)'); ylabel(sprintf('y_%d',i));
    legend('y','y_{hat}'); grid on;
end

figure('Name','Errore di stima (componenti dello stato)');
for i=1:3
    subplot(3,1,i);
    plot(tt, e(i,:), 'LineWidth', 1.4);
    xlabel('Tempo (s)'); ylabel(sprintf('e_%d',i));
    grid on;
end

% -------------------------
% 7) Stampa riepilogativa
% -------------------------
fprintf('\nRiepilogo:\n');
fprintf('Sistema: n = %d, m = %d, p = %d\n', size(A,1), size(B,2), size(C,1));
fprintf('Polo controller scelti: '); disp(poles_ctrl);
fprintf('Polo observer scelti:   '); disp(poles_obs);
fprintf('Autovalori chiusi (A-BK): '); disp(eig(A-B*K));
fprintf('Autovalori osservatore (A-LC): '); disp(eig(A-L*C));
fprintf('Autovalori complessivi (aug): '); disp(eig_aug);

