% Esercitazione per l'esame di Controllo avanzato e Applicazioni prof.
% Ambrosino Roberto


clear; close all; clc;

fprintf('__________ Salvatore Raiola __________ DE6000008\n')




%% ========================================================================
%  0. DEFINIZIONI 
%  ========================================================================


    % Parametri

        A1 = 1;
        A2 = 1;
        A3 = 1;

        a12 = 0.6;
        a23 = 0.4;
        a2 = 0.3;
        a3 = 0.25;



    % Variabili Simboliche

        syms h1 dh1 h2 dh2 h3 dh3 u1 u2 real



    % Stato e Ingresso

        x = [h1;
             h2;
             h3];
        
        u = [u1;
             u2];



    % Modello dinamico

        dh1_sym = (1/A1)*(u1 - a12*sqrt(h1-h2));
        dh2_sym = (1/A2)*(a12*sqrt(h1-h2) - a23*sqrt(h2-h3) - a2*sqrt(h2));
        dh3_sym = (1/A3)*(u2 + a23*sqrt(h2-h3) - a3*sqrt(h3));
        
        dx_sym = [dh1_sym;
                  dh2_sym;
                  dh3_sym];





%% ========================================================================
%  1. PUNTI DI EQUILIBRIO E LINEARIZZAZIONE 
%  ========================================================================


fprintf("\n\n\n========== 1. PUNTI DI EQUILIBRIO E LINEARIZZAZIONE ==========\n\n");


    % 1.1 Valori degli ingressi u1* e u2* che rendono x* un punto di equilibrio

        x_star = [3.274;
                  2;
                  1.6];

        eq_dh1 = subs(dh1_sym, {h1, h2}, {3.274, 2}) == 0;
        eq_dh3 = subs(dh3_sym, {h2, h3}, {2, 1.6}) == 0;

        sol_u1u2 = solve([eq_dh1, eq_dh3], [u1, u2], 'ReturnConditions',true);

        u1_star = double(sol_u1u2.u1);
        u2_star = double(sol_u1u2.u2);

        fprintf('1.1) Valori degli ingressi u1* e u2* che rendono x* un punto di equilibrio\n\n');
        fprintf('    u1*: %f\n', u1_star);
        fprintf('    u2*: %f\n\n', u2_star);

        u_star = [u1_star;
                  u2_star];

    % 1.2 Calcolo del modello linearizzato in (x*, u*)

        A_sym = jacobian(dx_sym, x);
        B_sym = jacobian(dx_sym, u);

        A = double(subs(A_sym, [x; u], [x_star; u_star]));
        B = double(subs(B_sym, [x; u], [x_star; u_star]));

        n = size(A, 1);
        C = eye(n);
        D = zeros(n, size(B, 2));

        fprintf('\n1.2) Modello linearizzato in (x*, u*)\n\n');
        disp('   Matrice A :'); disp(A);
        disp('   Matrice B :'); disp(B);
        disp('   Matrice C :'); disp(C);
        disp('   Matrice D :'); disp(D);


    % 1.3 Verifica stabilità, raggiungibilità e osservabilità
       
        fprintf('\n1.3) Verifica stabilità, raggiungibilità e osservabilità\n\n');
        autoval_A = eig(A);
        disp('   Autovalori di A (sistema linearizzato): '); disp(autoval_A);

        if all(real(autoval_A)<0)
            disp('  Il sistema è LOCALMENTE ASINTOTICAMENTE STABILE in (x*, u*)');
        elseif any(real(autoval_A) > 0)
            disp('  Il sistema è LOCALMENTE INSTABILE in (x*, u*)');
        else
            disp('  Caso con autovalori nulli o coincidenti. Non è possibile fare un analisi con il Metodo degli autovalori ');
        end

        % Raggiungibilità
        if rank(ctrb(A, B))==n
            disp('  Il sistema linearizzato è completamente Raggiungibile');
        else
            disp('  Il sistema linearizzato NON è completamente Raggiungibile');
        end

        % Osservabilità
        if rank(obsv(A,C)) == size(A, 1)
            disp('  Il sistema linearizzato è completamente Osservabile');
        else
            disp('  Il sistema linearizzato NON è completamente Osservabile');
        end




%% ========================================================================
%  2. CONTROLLO STABILIZZANTE MEDIANTE APPROCCIO ALLA LYAPUNOV (LMI)
%  ========================================================================

fprintf("\n\n\n========== 2. CONTROLLO STABILIZZANTE MEDIANTE APPROCCIO ALLA LYAPUNOV (LMI) ==========\n");
    
    % Definizioni preliminari

        x0 = [0.3;
              0.2;
              0.1];
    
        h1_ref = 0.5;
        h2_ref = 0.4;
    
        x_ref = [h1_ref;
                 h2_ref];
    
        m = size(B, 2);


    % Calcolo di L e Q con cvx

        cvx_begin sdp quiet
            variable Q(n, n) symmetric
            variable L(m, n)
            
            % feasibility
            minimize(0) 
            
            subject to
                % Vincolo 1: Q definita positiva (LMI: Q > 0)
                Q >= 1e-6 * eye(n);
                
                % Vincolo 2: Diseguaglianza di Lyapunov per il ciclo chiuso
                A*Q + Q*A' - B*L - L'*B' <= -1e-6 * eye(n);
        cvx_end
    

    % Calcolo del guadagno K
    
        if strcmp(cvx_status, 'Solved')
            K = L/Q;
            fprintf('\n2.1) Guadagno K calcolato:\n\n');
            disp(K);
        else
            error('CVX non ha trovato una soluzione. Verifica la controllabilità del sistema.');
        end


    % sistema a ciclo chiuso

        A_cl = A - B*K; 
    

    % Definizione della matrice C di "uscita controllata"

        C_ctrl = [1 0 0;  % Estrae h1
                  0 1 0]; % Estrae h2
        

    % Calcolo del guadagno statico a ciclo chiuso

        G_statico = C_ctrl * (-(A - B*K) \ B);

        
    % Calcolo della matrice N (compensazione)

        if det(G_statico) ~= 0
            N = inv(G_statico);
            fprintf('\n2.2) Matrice di compensazione N calcolata:\n\n');
            disp(N);
        else
            error('G_statico non invertibile. Controlla il modello.');
        end
                
    
    % Apertura automatica del modello

        modello_verifica = 'EsercitazioneEsame_LMI_sim';

        fprintf('\n2.3) Apriamo il modello Simulink per verificare il funzionamento del controllore...\n\n')
        if exist(modello_verifica, 'file')
            open_system(modello_verifica);
            disp(['    Modello ', modello_verifica, ' aperto con successo.']);
        else
            error('File Simulink non trovato. Assicurati che sia nella cartella corrente.');
        end




%% ========================================================================
%  3. CONTROLLO LQR CON AZIONE INTEGRALE 
%  ========================================================================



fprintf("\n\n\n========== 3. CONTROLLO LQR CON AZIONE INTEGRALE ==========\n");


    % Definizioni preliminari
        
        fprintf('\n3.1) Sistema Aumentato\n\n');

        h2_LQR_ref = 2.8;
        t_ass = 1;

        Ch2 = [0 1 0];

        A_aum = [A    ,  zeros(n,1); 
                 -Ch2 ,          0];
        fprintf('   A_aum:\n');
        disp(A_aum);


        fprintf('\n   B_aum:\n');
        B_aum = [B;
                 zeros(1,m)];
        disp(B_aum);


    % Scelta Q e R
        
        fprintf('\n3.2) Scelta Q e R\n\n');

        Q = diag([1, 10, 1, 100]);  % 4x4
        fprintf('   Q:\n');
        disp(Q);

        R = eye(m) * 0.1;
        fprintf('   R:\n');
        disp(R);


    % Sintesi LQR
        fprintf('\n3.3) Sintesi LQR\n\n');

        K_LQR = lqr(A_aum, B_aum, Q, R);
        fprintf('   K_LQR:\n');
        disp(K_LQR);

    % Guadagni

        fprintf('\n3.4) Guadagni\n\n');

        K_LQR_x = K_LQR(:, 1:3); % Guadagno per la retroazione dello stato (2x3)
        fprintf('   Guadagno per la retroazione dello stato:\n');
        disp(K_LQR_x);

        K_LQR_i = - K_LQR(:, 4);   % Guadagno per l'azione integrale (2x1)
        fprintf('   Guadagno per azione integrale:\n')
        disp(K_LQR_i);
        
        




        

