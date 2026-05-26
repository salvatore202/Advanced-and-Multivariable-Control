%% Progettazione di un Filtro di Kalman in Tempo Continuo
% Sistema reale massa-molla-smorzatore
% Live Script MATLAB (.mlx)

%% Descrizione del sistema
% m x¨ + c x˙ + k x = u + w
% y = x + v

m = 1;      % massa [kg]
c = 0.5;    % smorzamento [Ns/m]
k = 2;      % rigidezza [N/m]

%% Modello in forma di stato
A = [0 1;
    -k/m -c/m];

B = [0;
     1/m];

C = [1 0];
D = 0;

sys = ss(A,B,C,D);

%% Modello dei rumori
% Rumore di processo e di misura
Q = [0.01 0;
     0    0.1];

R = 0.05;

%% Calcolo del guadagno di Kalman (tempo continuo)
% Risoluzione dell'equazione di Riccati algebrica continua
[P,~,~] = care(A', C', Q, R);
K = P*C'/R;

K

%% Dinamica dello stimatore di Kalman
% x_hat_dot = (A - K C) x_hat + B u + K y

A_kf = A - K*C;
B_kf = [B K];   % ingressi: [u  y]
C_kf = eye(2);  % stimo posizione e velocità
D_kf = zeros(2,2);

KF = ss(A_kf, B_kf, C_kf, D_kf);



%% Simulazione
T = 0:0.01:10;

u = sin(T);                     % ingresso noto
x_true = lsim(sys, u, T);        % stato reale (senza rumore)
y = x_true(:,1) + 0.2*randn(size(T')); % misura rumorosa

inputKF = [u' y];

x_hat = lsim(KF, inputKF, T);

%% Grafico: posizione
figure
plot(T, x_true(:,1), 'k', 'LineWidth', 1.5)
hold on
plot(T, y, 'r--')
plot(T, x_hat(:,1), 'b', 'LineWidth', 1.5)
grid on
legend('Posizione reale','Misura rumorosa','Stima Kalman')
xlabel('Tempo [s]')
ylabel('x [m]')
title('Filtro di Kalman continuo - posizione')

%% Grafico: velocità stimata
figure
plot(T, x_hat(:,2), 'b', 'LineWidth', 1.5)
grid on
xlabel('Tempo [s]')
ylabel('dx/dt [m/s]')
title('Velocità stimata (non misurata)')

%% Controllo ottimo LQR (tempo continuo)
% Obiettivo: progettare un controllo di stato u = -L x
% che minimizzi il costo:
% J = ∫ (x' Q_lqr x + u' R_lqr u) dt

Q_lqr = [10 0;
         0  1];

R_lqr = 0.5;

%% Calcolo del guadagno LQR
[L, S, eigCL] = lqr(A, B, Q_lqr, R_lqr);

L

%% Sistema in anello chiuso (stato reale)
A_cl = A - B*L;

sys_cl = ss(A_cl, B, C, D);

%% Principio di separazione: LQG
% Controllo basato sulla stima:
% u = -L x_hat

A_lqg = [A - B*L,        B*L;
         zeros(size(A)), A - K*C];

B_lqg = [B;
         B];

C_lqg = [C zeros(size(C))];

D_lqg = 0;

sys_lqg = ss(A_lqg, B_lqg, C_lqg, D_lqg);

%% Simulazione LQG
T = 0:0.01:10;

x0 = [0.5; 0];            % stato iniziale
xhat0 = [0; 0];           % stima iniziale

x_init = [x0; xhat0];

u = zeros(size(T));       % nessun riferimento esplicito

[y_lqg, T, x_lqg] = lsim(sys_lqg, u, T, x_init);

x_real = x_lqg(:,1:2);
x_hat  = x_lqg(:,3:4);

%% Grafico: controllo LQG
figure
plot(T, x_real(:,1), 'k', 'LineWidth', 1.5)
hold on
plot(T, x_hat(:,1), 'b--', 'LineWidth', 1.5)
grid on
legend('Posizione reale','Posizione stimata')
xlabel('Tempo [s]')
ylabel('x [m]')
title('Controllo LQG (LQR + Kalman)')

%% Osservazioni
% - Il controllo LQR stabilizza il sistema
% - Il filtro di Kalman fornisce la stima degli stati
% - Il principio di separazione garantisce stabilità
% - LQG = LQR + Kalman continuo
