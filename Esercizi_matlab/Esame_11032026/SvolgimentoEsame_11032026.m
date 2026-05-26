% Svolgimento esame del 11.03.2026

clear; close all; clc;

fprintf('__________ Salvatore Raiola __________ DE6000008\n')

%% 
%  =======================================================================
%  Definizioni 
%  =======================================================================

% valori
m1 = 1.0;
m2 = 1.2;
k1 = 20;
k2 = 25;
c1 = 1;
c2 = 1.2;
k = 15;
c = 0.8;
k3 = 40;

% variabili simboliche 
syms q1 dq1 ddq1 q2 dq2 ddq2 u1 u2 real

% Stato e ingresso
x = [q1;
     q2;
     dq1;
     dq2];

u = [u1;
     u2];

ddq1_sym = (1/m1)*(-k1*q1 - c1*dq1 + k*(q2-q1) + k3*(q2-q1)^3 + c*(dq2-dq1) + u1);
ddq2_sym = (1/m2)*(-k2*q2 - c2*dq2 - k*(q2-q1) - k3*(q2-q1)^3 - c*(dq2-dq1) + u2);

ddx_sym = [ddq1_sym;
           ddq2_sym];



%% 
% =======================================================================
% 1.1 - Punti di equilibrio
% =======================================================================

fprintf('\n\n\n=========== 1.1 - Punti di equilibrio ===========\n');

eq_ddq1 = subs(ddq1_sym, {dq1, dq2, u1}, {0,0,0}) == 0;
eq_ddq2 = subs(ddq2_sym, {dq1, dq2, u2}, {0,0,0}) == 0;

% Risolviamo il sistema di equazioni rispetto a q1 e q2
sol_q1q2 = solve([eq_ddq1, eq_ddq2], [q1, q2], 'ReturnConditions',true);

fprintf('\nPunti di equilibrio:\n');
disp(sol_q1q2);



%% 
% ========================================================================
% 1.2 - Punto di equilibrio x* con u1*=1N e u2*=-0.5N
% ========================================================================

fprintf('\n\n\n============ 1.2 - Punto di equilibrio x* con u1*=1N e u2*=-0.5N ============\n');

eq_ddq1 = subs(ddq1_sym, {dq1, dq2, u1}, {0,0,1}) == 0;
eq_ddq2 = subs(ddq2_sym, {dq1, dq2, u2}, {0,0,-0.5}) == 0;

% Risolviamo il sistema di equazioni rispetto a q1* e q2*
sol_q1q2 = solve([eq_ddq1, eq_ddq2], [q1, q2]);

fprintf('\nPunti di equilibrio per u1*=1N e u2*=-0.5N:\n');
q1_star = double(sol_q1q2.q1);
q2_star = double(sol_q1q2.q2);
fprintf('            q1*: %f', q1_star);
fprintf('\n            q2*: %f\n', q2_star); 

x_star = {q1_star;
          0;
          q2_star;
          0};



%% 
% ========================================================================
% 1.3 - Linearizzazione
% ========================================================================

fprintf('\n\n\n============ 1.3 - Linearizzazione ============\n');


f = [dq1; 
     ddq1_sym;
     dq2;
     ddq2_sym];

A_sym = jacobian(f, x);
B_sym = jacobian(f, u);

A = double(subs(A_sym, x, x_star));
B = double(subs(B_sym, x, x_star));

fprintf('\n');
disp('Matrice A (Sist. Linearizzato):'); disp(A);
disp('Matrice B (Sist. Linearizzato):'); disp(B);



%% 
% ========================================================================
% 1.4 - Verifica Controllabilità e osservabilità
% ========================================================================

fprintf('\n\n\n============ 1.4 - Verifica Controllabilità e osservabilità ============\n\n');

C = eye(4);

D = zeros(4, 2);

% Controllabilità
if rank(ctrb(A, B)) == size(A, 1)
    disp('Il sistema linearizzato è completamente controllabile');
else
    disp('Il sistema linearizzato NON è completamente controllabile');
end

% Osservabilità
if rank(obsv(A, C)) == size(A, 1)
    disp('Il sistema linearizzato è completamente osservabile');
else
    disp('Il sistema linearizzato NON è completamente osservabile');
end



%% 
% ========================================================================
% 2 - Assegnamento autovalori
% ========================================================================

fprintf('\n\n\n============ 2 - Assegnamento autovalori ============\n\n');

% autovalori desiderati
eps = 1e-4; % Piccola perturbazione numerica per usare il comando place 
p = [-5, -5-eps, -10, -10-eps];

% parametri schema Simulink
q1_0 = 0.1;
q2_0 = 0.3;
q1_ref = 0.5;
q2_ref = 1;

x_0 = [q1_0;
       q2_0;
       0;
       0];

x_ref = [q1_ref;
         q2_ref];
        

% guadagno
K = place(A, B, p);
disp('Matrice di Guadagno:');
disp(K);

% verifica dei poli effettivamente ottenuti
disp('Poli del sistema a ciclo chiuso:');
disp(eig(A - B*K));

% compensazione del guadagno
A_cl = A - B*K;        % Matrice a ciclo chiuso
C_N = [1, 0, 0, 0;
       0, 1, 0, 0];
N = -inv(C_N * (A_cl \ B));


% Apertura automatica del modello
modello_verifica = 'VerificaAssegnamentoAutoval_11032026';

disp('Apriamo il modello Simulink per verificare il funzionamento del controllore...')
if exist(modello_verifica, 'file')
    open_system(modello_verifica);
    disp(['Modello ', modello_verifica, ' aperto con successo.']);
else
    error('File Simulink non trovato. Assicurati che sia nella cartella corrente.');
end

