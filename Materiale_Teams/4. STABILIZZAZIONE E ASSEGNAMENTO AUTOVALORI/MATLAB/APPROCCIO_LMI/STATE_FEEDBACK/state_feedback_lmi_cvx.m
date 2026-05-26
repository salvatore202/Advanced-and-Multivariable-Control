% State-feedback stabilizzante via LMI (variabili W,N) usando CVX
% Sistema d'esempio: ordine 3, 2 ingressi, 2 uscite
% Richiede CVX.

clear; close all; clc;

% --- Sistema 
A = [0 1 0;
     0 0 1;
    2 3 -4];
B = [0 0;
     0 1;
     1 0];
C = [1 0 0;
     0 0 1];    % non usata per questo esempio
[n, m] = size(B);

fprintf('State-feedback LMI: n=%d, m=%d\n', n, m);

% Parametri LMI
epsLMI = 1e-6;

% Risolvo LMI per (W,N)
cvx_clear;
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
% ricostruisci K
K = Nval / Qval;

fprintf('K trovato (m x n):\n'); disp(K);
eig_cl = eig(A + B*K);
fprintf('Autovalori A+B*K:\n'); disp(eig_cl.');
fprintf('Costanti di tempo A+B*K:\n'); disp(-1./real(eig_cl)');

% Simulazione: risposta a condizione iniziale (solo per verificare stabilità)
t = 0:0.01:5;
x0 = [1; -0.5; 0.2];
% dinamica chiusa: xdot = (A + B*K)*x
sys_cl = ss(A + B*K, [], eye(n), []);
[y,t,x] = initial(sys_cl, x0, t);

figure('Name','State-feedback: stati (chiuso)');
for i=1:n
    subplot(n,1,i);
    plot(t, x(:,i), 'LineWidth', 1.4); grid on;
    ylabel(sprintf('x_%d', i));
    if i==1, title('Stati in anello chiuso con state-feedback da LMI'); end
end
xlabel('Tempo (s)');

fprintf('Esecuzione completata.\n');
