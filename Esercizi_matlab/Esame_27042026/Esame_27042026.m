clear; close all; clc;

fprintf('\nCONTROLLO AVANZATO E APPLICAZIONI ____________________ ESAME DEL 27/04/2026 ____________________ SALVATORE ___ RAIOLA ___ DE6000008')


%% 0) DEFINIZIONI 

    % Variabili simboliche 

        syms x1 x2 theta dx1 dx2 dtheta ddx1 ddx2 ddtheta u1 u2 real
        

    % Parametri

        m1 = 1.5;
        m2 = 1.0;
        mp = 0.25;
        l = 0.40;
        J = 0.013;
        k = 35;
        b = 4;
        g = 9.81;


    % Ingresso, Stato e Uscita

        u = [u1;
             u2];

        x = [x1;
             x2;
             theta;
             dx1;
             dx2;
             dtheta];

        y = [x1;
             x2;
             theta];


    % Modello Dinamico non lineare

        eq1 = (m1+mp)*ddx1 + mp*l*cos(theta)*ddtheta - mp*l*sin(theta)*dtheta^2 + k*(x1-x2) + b*(dx1-dx2) == u1;
        eq2 = m2*ddx2 + k*(x2-x1) +b*(dx2-dx1) == u2;
        eq3 = (J+mp*l^2)*ddtheta + mp*l*cos(theta)*ddx1 - mp*g*l*sin(theta) == 0;

        eqs = [eq1;
               eq2;
               eq3];

        sol = solve(eqs, [ddx1, ddx2, ddtheta]); 

        dx_sym = [dx1;
                  dx2;
                  dtheta;
                  sol.ddx1;
                  sol.ddx2;
                  sol.ddtheta];

        eq_ddx1_sym = sol.ddx1;
        eq_ddx2_sym = sol.ddx2;
        eq_ddtheta_sym = sol.ddtheta;

        f = dx_sym;



 %% 1. PUNTI DI EQUILIBRIO E LINEARIZZAZIONE


fprintf('\n====================== 1. PUNTI DI EQUILIBRIO E LINEARIZZAZIONE ======================\n');


    % 1.1) Punti di equilibrio del sistema 

        eq_ddx1 = subs(eq_ddx1_sym, {dx1, dx2, dtheta, u1, u2}, {0, 0, 0, 0, 0}) == 0;
        eq_ddx2 = subs(eq_ddx2_sym, {dx1, dx2, dtheta, u1, u2}, {0, 0, 0, 0, 0}) == 0;
        eq_ddtheta = subs(eq_ddtheta_sym, {dx1, dx2, dtheta, u1, u2}, {0, 0, 0, 0, 0}) == 0;
        
        % Fissiamo x1 = 0 per eliminare la traslazione libera e trovare x2 e theta
        sol_peql = solve([eq_ddx1, eq_ddx2, eq_ddtheta, x1 == 0], [x1, x2, theta], 'ReturnConditions',true);

        x1_bar = sol_peql.x1;
        x2_bar = sol_peql.x2;
        theta_bar = sol_peql.theta;
        dx1_bar = 0;
        dx2_bar = 0;
        dtheta_bar = 0;

        fprintf('1.1) Punti di equilibrio del sistema\n\n'); 

        disp(table(x1_bar, x2_bar, theta_bar, dx1_bar, dx2_bar, dtheta_bar));



    % 1.2) Sistema linearizzato nel punto di lavoro con x=0 e u=0

        x_star = [0;
                  0;
                  0;
                  0;
                  0;
                  0];
        
        u_star = [0;
                  0];

        A_sym = jacobian(f, x);
        disp('Matrice A (simbolica):');
        A_sym

        B_sym = jacobian(f, u);
        disp('Matrice B (simbolica):');
        B_sym

        A = double(subs(A_sym, [x; u], [x_star; u_star]));
        B = double(subs(B_sym, [x; u], [x_star; u_star]));
        
        n = size(A, 1);
        m = size(B, 2);

        C = eye(n);

        D = zeros(n, m);

        poles = eig(A);

        fprintf('\n\n1.2) Sistema linearizzato attorno al punto di lavoro corrispondente a x=0 e u=0\n\n');
        A
        B
        C
        D

        disp('Autovalori A: '); fprintf('\n'); disp(poles);



    % 1.3) Verifica controllabilità e Osservabilità

        fprintf('\n\n1.3) Verifica controllabilità e Osservabilità\n\n');

        if rank(ctrb(A, B)) == n
            disp(' Il Sistema lineaarizzato è completamente controllabile');
        else
            disp('  Il Sistema lineaarizzato NON è completamente controllabile');
        end

        fprintf('\n');

        if rank(obsv(A, C)) == n
            disp(' Il Sistema lineaarizzato è completamente osservabile');
        else
            disp('  Il Sistema lineaarizzato NON è completamente osservabile');
        end




    

%% 2. APPROCCIO ALLA LYAPUNOV


 fprintf('\n====================== 2. APPROCCIO ALLA LYAPUNOV ======================\n');


    % 2.1) Risoluzione LMI

        cvx_clear;
        cvx_begin sdp quiet
            cvx_precision high
    
            variable Q(n,n) symmetric
            variable L(m,n)
    
            minimize(0)
    
            subject to
                Q >= 1e-6*eye(n);
                A*Q + Q*A' + B*L + L'*B' <= -1e-6*eye(n);
        cvx_end
    
        if ~(strcmp(cvx_status, 'Solved'))
            disp('CVX non ha trovato soluzione per: '); disp(cvx_status);
        end
    
        Qval = Q;
        Lval = L;
    
        K_LMI = L/Q;
    
        disp('K (stabilizzazione con LMI):');
        fprintf('\n');
        disp(K_LMI);


    % 2.2) Verifica MATLAB

        Acl_LMI = A+B*K_LMI;

        poles_cl_LMI = eig(Acl_LMI);

        disp('Autovalori a ciclo chiuso: '); 
        fprintf('\n');
        disp(poles_cl_LMI);

        sys_cl_LMI = ss(Acl_LMI, B, C, D);
        
        % supponiamo uno stato iniziale diverso da 0
        x0_LMI = [1;
                  2;
                  3;
                  4;
                  5;
                  6];
       
        figure;
        initial(sys_cl_LMI, x0_LMI);
        grid on;
        title('Verifica Convergenza stato - Stabilizzazione LMI');





%% 3. ASSEGNAMENTO DEGLI AUTOVALORI


 fprintf('\n====================== 3. ASSEGNAMENTO DEGLI AUTOVALORI ======================\n');


    % Autovalori desiderati per avere tempo di assestamento < 10s

        eig_des = [-1.5, -1.6, -1.7, -1.8, -1.9, -2.0];


    % Condizioni iniziali

        x0_EIG = [0.05;
                  -0.05;
                  10;
                  0;
                  0;
                  0];


    % Valori di riferimento costante

        xref_EIG = [0.15;
                    0.15];


    % 3.1) Guadagno del controllore 

        K_EIG = place(A, B, eig_des);
        K_EIG


    % 3.2) Parametri per la simulazione 

        C_inseguimento_EIG = [1, 0, 0, 0, 0, 0;
                              0, 1, 0, 0, 0, 0];

        % Compensazione
        N = inv(C_inseguimento_EIG*(-inv(A-B*K_EIG))*B);

        disp('Compensazione:'); 
        N        

        fprintf('\n');

        fprintf('\n\nVerificare il risultato aprendo il file Simulink: VerificaAssegnamentoAutovalori_27042026.slx\n\n');




%% 4. CONTROLLO LQR CON AZIONE INTEGRALE


 fprintf('\n====================== 4. CONTROLLO LQR CON AZIONE INTEGRALE ======================\n');


    % Sistema Aumentato

        C_LQR = C_inseguimento_EIG;
    
        A_a = [A      , zeros(6,2);
               - C_LQR, zeros(2,2)];
    
        B_a = [B;
               zeros(2,2)];
    
        D_a = zeros(8,2);
    
        fprintf('4.1) Sistema Aumentato: \n\n');
        A_a
        B_a
        C_LQR
        D_a


    % Sintesi LQR con Tempo di Assestamento specificato

        Ts = 5;

        alpha = 4/Ts;

        A_a_hat = A_a + alpha*eye(size(A_a));

        q_i_stato = 1;
        q_1_integrali = 100;

        Q_a = diag([q_i_stato, q_i_stato, q_i_stato, q_i_stato, q_i_stato, q_i_stato, q_1_integrali, q_1_integrali]);
        Q_a

        R_a = 1*eye(2);
        R_a

        [K_LQR, P, e_bar] = lqr(A_a_hat, B_a, Q_a, R_a);

        K_LQR


    % Verifica 

        A_cl_LQR = A_a - B_a*K_LQR;

        poles_cl_LQR = eig(A_cl_LQR);
        disp('Poli a ciclo chiuso:'); disp(poles_cl_LQR);


    % Guadagni
        
        K_LQR_x = K_LQR(:, 1:n);       % Prendi tutte le righe, colonne da 1 a 3
        K_LQR_i = K_LQR(:, n+1:end);   % Prendi tutte le righe, colonne da 4 a 5

        
     
                 
                   
                   