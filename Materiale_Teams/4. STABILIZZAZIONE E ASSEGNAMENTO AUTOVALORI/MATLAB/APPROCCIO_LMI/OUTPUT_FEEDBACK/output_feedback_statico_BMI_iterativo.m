clear; close all; clc;

%% ---------------------------
% 1) Definizione del sistema
% ----------------------------
A = [0 1 0;
     0 0 1;
    -2 3 -4];
B = [0 0;
     0 1;
     1 0];
C = [1 0 0;
     0 0 1];

[n, m] = size(B);
[p, ~] = size(C);

fprintf('Sistema: ordine n = %d, ingressi m = %d, uscite p = %d\n', n, m, p);

%% ---------------------------
% 2) Parametri iterazione
% ----------------------------
maxIter = 50;
tolK = 1e-4;
epsLMI = 1e-6;


%% ---------------------------
% 3) Loop iterativo
% ----------------------------
for it = 1:maxIter
    fprintf('\n--- Iterazione %d ---\n', it);

    % STEP 1: Fissa un P definita positiva

    P=diag(randn(n,1));
    P=P*P';

lim = 10;
nn = 60;
x = linspace(-lim, lim, nn);
y = linspace(-lim, lim, nn);
z = linspace(-lim, lim, nn);
[X, Y, Z] = meshgrid(x, y, z);
% ---- Calcolo di f(x) = x' P x ----
F = P(1,1)*X.^2 + P(2,2)*Y.^2 + P(3,3)*Z.^2 ...
    + 2*P(1,2)*X.*Y + 2*P(1,3)*X.*Z + 2*P(2,3)*Y.*Z;
level = 1;  % scegli il valore costante per l'isolevel
figure;
isosurface(X, Y, Z, F, level);
grid
axis equal
keyboard

    %P=eye(n);

    cvx_begin sdp quiet
        variable K(m,p) 
        minimize(0) % nessuna funzione obiettivo, solo feasibility
        subject to
            (A + B*K*C)'*P + P*(A + B*K*C) <= -epsLMI * eye(n);
    cvx_end

    if ~strcmp(cvx_status, 'Solved')
        warning('Iter %d: impossibile trovare P per K corrente (status: %s).', it, cvx_status);
    end
    Pval = P;

    % Verifica stabilità del sistema chiuso
    eigs_cl = eig(A + B*K*C);
    fprintf('Autovalori(A + BKC):\n');
    disp(eigs_cl.');

    if all(real(eigs_cl) < 0)
        fprintf('✅ Sistema stabilizzato!\n');
    else
        fprintf('❌ Sistema ancora instabile.\n');
    end
end

%% ---------------------------
% 4) Risultati finali
% ----------------------------
if all(real(eig(A + B*K*C)) < 0)
    fprintf('\n=== RISULTATO FINALE ===\n');
    fprintf('K stabilizzante:\n');
    disp(K);
    fprintf('Autovalori del sistema chiuso:\n');
    disp(eig(A + B*K*C));
else
    warning('Nessuna soluzione stabilizzante trovata dopo %d iterazioni.', it);
end
