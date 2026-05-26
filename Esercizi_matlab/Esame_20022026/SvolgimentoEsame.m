

% Svolgimento traccia d'esame - Controllo Avanzato e Applicazioni - appello
% del 20/02/2026

clear; close all; clc;

%% Definizioni e inizializzazione
    
%definiamo i valori dei parametri
M = 0.8;      % Massa carrello
m = 0.2;      % Massa pendolo
l = 0.3;      % Lunghezza asta (centro di massa)
I = 0.006;    % Inerzia baricentrica
g = 9.81;     % Gravità
b = 0.08;     % Attrito carrello
c = 0.002;    % Attrito giunto
J = I + m*l^2; % Inerzia totale rispetto al pivot

% Definiamo le variabili simboliche
syms x dx ddx theta dtheta ddtheta u real

matr_M = [ (M + m)    , m*l*cos(theta); 
          m*l*cos(theta) ,     J      ];

vett_F = [ u - b*dx + m*l*sin(theta)*dtheta^2; 
            - c*dtheta + m*g*l*sin(theta)   ];

acc = matr_M\vett_F;

ddx_sym = acc(1);  % la prima riga di acc
ddtheta_sym = acc(2);  % la seconda riga di acc




%% 1.1 - Determinazione punti di equilibrio considerando u=0

% 1. Poniamo le velocità a zero: dx = 0, dth = 0
% 2. Poniamo l'ingresso a zero: u = 0
% 3. Cerchiamo dove le accelerazioni (ddx, ddth) si annullano

eq_ddx = subs(ddx_sym, {dx, dtheta, u}, {0, 0, 0}) == 0;  % sostituiamo 0 a dx, dtheta e u nell'espressione di ddx
eq_ddtheta = subs(ddtheta_sym, {dx, dtheta, u}, {0, 0, 0}) == 0;

% Risolviamo il sistema di equazioni rispetto a theta
sol_theta = solve([eq_ddx, eq_ddtheta], theta, 'ReturnConditions',true);

disp('Punti di equilibrio per theta:');
disp(sol_theta);



%% 1.2 - Linearizzazione

% Definiamo lo stato X = [x; dx; th; dth]
X = [x; dx; theta; dtheta];

% Definiamo la funzione di stato f = [x_dot; x_ddot; th_dot; th_ddot]
f = [dx; ddx_sym; dtheta; ddtheta_sym];

% Matrice A: derivata di f rispetto allo stato X
A_sym = jacobian(f, X);

% Matrice B: derivata di f rispetto all'ingresso u
B_sym = jacobian(f, u);

p_eq_instabile = {x, dx, theta, dtheta, u};
val_instabili = {0, 0, 0, 0, 0};

A = double(subs(A_sym, p_eq_instabile, val_instabili));
B = double(subs(B_sym, p_eq_instabile, val_instabili));

disp('Matrice A (Sist. Linearizzato):'); disp(A);
disp('Matrice B (Sist. Linearizzato):'); disp(B);



%% 1.3 - Verifica Osservabilità e Controllabilità

% Osservabilità
C = eye(4);      % Ora l'uscita ha 4 elementi (tutto lo stato)
D = zeros(4, 1); % D deve avere lo stesso numero di righe di C

% Matrice di Osservabilità: [C; C*A; C*A^2; C*A^3]
Ob = obsv(A, C);

fprintf('\nRango matrice Osservabilità: %d\n', rank(Ob));

if rank(Ob) == size(A, 1)
    disp('Il sistema linearizzato è completamente osservabile');
else
    disp('Il sistema linearizzato NON è completamente osservabile');
end

% Controllabilità

Co = ctrb(A, B);

fprintf('\nRango matrice Controllabilità: %d\n', rank(Co));

if rank(Co) == size(A, 1)
    disp('Il sistema linearizzato è completamente controllabile');
else
    disp('Il sistema linearizzato NON è completamente controllabile');
end



%% 2 - Assegnamento degli autovalori

% Condizioni iniziali (convertiamo i gradi in radianti!)
x0 = 0.1;           % m
dx0 = 0;            % m/s
th0 = 15 * pi/180;  % rad (15 gradi)
dth0 = 0;           % rad/s

% Vettore stato iniziale
X_initial = [x0; dx0; th0; dth0];

% Vettore degli autovalori desiderati
p = [-2, -2.5, -3, -3.5];


% Matrice di guadagno
K_2 = place(A, B, p);

disp('Matrice di guadagno K calcolata:');
disp(K_2);

% Verifica degli autovalori del sistema a ciclo chiuso
autovalori_chiuso = eig(A - B*K_2);
disp('Autovalori a ciclo chiuso (verifica):');
disp(autovalori_chiuso);

% --- Apertura automatica del modello ---
modello_verifica = 'VerificaAssegnamentoAutoval';

if exist(modello_verifica, 'file')
    open_system(modello_verifica);
    disp(['Modello ', modello_verifica, ' aperto con successo.']);
else
    error('File Simulink non trovato. Assicurati che sia nella cartella corrente.');
end



%% 3 - Approccio alla Lyapunov (richiede cvx)

% Stato desiderato
x_ref_val = 0.2;
X_ref = [x_ref_val; 0; 0; 0];

[n, m] = size(B);

epsLMI = 1e-6; % Margine di stabilità
cvx_clear;      % Pulisce variabili CVX precedenti
cvx_begin sdp quiet
    cvx_precision high
     variable Q(n,n) symmetric
    variable N(m,n)
     minimize(0) % feasibility
     subject to
        Q >= epsLMI*eye(n);
        A*Q + Q*A' + B*N + N'*B' <= -epsLMI*eye(n);
cvx_end

if ~strcmp(cvx_status,'Solved')
    error('CVX non ha trovato soluzione per (W,N): %s', cvx_status);
end

Qval = Q;
Nval = N;

% ricostruzione K
K_3 = Nval / Qval;
fprintf('K trovato (m x n):\n'); disp(K_3);

% autovalori
eig_cl = eig(A + B*K_3);
fprintf('Autovalori A+B*K:\n'); disp(eig_cl.');