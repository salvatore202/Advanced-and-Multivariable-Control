%[text] # PROF. AMBROSINO R.\_\_\_ESAME\_\_\_20/04/2026\_\_\_SALVATORE\_RAIOLA\_DE6000008
%[text] 
clear; close all; clc;
%[text] ## 
%%
%[text] ## 0)  DEFINIZIONI 
%[text] #### 
%[text] #### Variabili simboliche 
syms h1 h2 h3 dh1 dh2 dh3 u1 u2 w real
%[text] #### 
%[text] #### Parametri
A1 = 1;
A2 = 1;
A3 = 1;

a12 = 0.6;
a23 = 0.4;
a2 = 0.3;
a3 = 0.25;
%[text] #### 
%[text] #### Modello Dinamico 
% Ingresso
u = [u1;
     u2];

% Stato
x = [h1;
     h2;
     h3];

% Uscita
y = x;

% Equazioni del sistema
eq_1 = -A1*dh1 + u1 - a12*sqrt(h1-h2) == 0;
eq_2 = -A2*dh2 + a12*sqrt(h1-h2) - a23*sqrt(h2-h3) - a2*sqrt(h2) == 0;
eq_3 = -A3*dh3 + u2 + a23*sqrt(h2-h3) - a3*sqrt(h3) == 0;

eqs = [eq_1;
       eq_2;
       eq_3];

sol = solve(eqs, [dh1, dh2, dh3], "ReturnConditions",true);

% dx/dt
f = [sol.dh1;
     sol.dh2;
     sol.dh3];
% y
g = x;
%[text] ## 
%%
%[text] ## 1)  PUNTI DI EQUILIBRIO E LINEARIZZAZIONE
%[text] 
eq_dh1_sym = sol.dh1;
eq_dh2_sym = sol.dh2;
eq_dh3_sym = sol.dh3;

x_star_1 = 3.274;
x_star_2 = 2;
x_star_3 = 1.6;

% Punto di Equilibrio
x_star = [x_star_1;
          x_star_2;
          x_star_3];
%[text] ### 
%[text] ### 1.1)  Valori di u1\* e u2\* che rendono x\* un punto di equilibrio 

eq_dh1 = subs(eq_dh1_sym, {h1, h2, h3}, {x_star_1, x_star_2, x_star_3}) == 0;
eq_dh3 = subs(eq_dh3_sym, {h1, h2, h3}, {x_star_1, x_star_2, x_star_3}) == 0;

u_star_sol = solve([eq_dh1, eq_dh3], [u1, u2], 'ReturnConditions',true);

u1_star = double(u_star_sol.u1);
u2_star = double(u_star_sol.u2);

disp(table(x_star_1, x_star_2, x_star_3, u1_star, u2_star)); %[output:59b43b97]
%[text] 
%[text] ### 
%[text] ### 1.2)  Modello Linearizzato nell'intorno di (x\*, u\*)
u_star = [u1_star;
          u2_star];

A_sym = jacobian(f, x);
B_sym = jacobian(f, u);
C_sym = jacobian(g, x);
D_sym = jacobian(g, u);
%[text] ### 
%[text] ### 1.3)  Matrici del Sistema
A = double(subs(A_sym, [x; u], [x_star; u_star]));
B = double(subs(B_sym, [x; u], [x_star; u_star]));
C = double(subs(C_sym, [x; u], [x_star; u_star]));
D = double(subs(D_sym, [x; u], [x_star; u_star]));
n = size(A,1);
m = size(B,2);

sys = ss(A, B, C, D) %[output:49bbd0e2]
%[text] ### 
%[text] 
%[text] ### 1.4) Verifica Stabilità x\*, raggiungibilità e osservabilità sistema
%[text] #### 
%[text] #### Stabilità
% Approccio agli autovalori

autovalori_A = eig(A);

if (all(autovalori_A < 0)) %[output:group:7a3d89c6]
    disp('Il punto di equiibrio è asintoticamente stabile'); autovalori_A %[output:57a57d30] %[output:0a0396e9]
elseif (any(autovalori_A == 0))
    disp('Il punto di equilibrio è stabile'); autovalori_A
else
    disp('Il punto di equilibrio è instabile'); autovalori_A
end %[output:group:7a3d89c6]
%[text] #### 
%[text] #### Raggiungibilità
if (rank(ctrb(A,B)) == n) %[output:group:84df5856]
    disp('Il sistema è completamente raggiungibile'); %[output:2bda1780]
else
    disp('Il sistema NON è completamente raggiungibile');
end %[output:group:84df5856]
%[text] #### 
%[text] #### Osservabilità
if (rank(obsv(A,C)) == n) %[output:group:1bbeaa31]
    disp('Il sistema è completamente osservabile'); %[output:732c43f7]
else
    disp('Il sistema NON è completamente osservabile');
end %[output:group:1bbeaa31]
%[text] 
%[text] 
%%
%[text] ## 2)  CONTROLLO STABILIZZANTE MEDIANTE APPROCCIO ALLA LYAPUNOV
%[text] #### 
%[text] #### Specifiche
% Stato iniziale
h0 = [0.3;
      0.2;
      0.1];

C_sel = C(1:2,:);

% Riferimento
h1_ref_lmi = 0.5;
h2_ref_lmi = 0.4;

h_ref = [h1_ref_lmi;
         h2_ref_lmi];
%[text] #### 
%[text] #### Risoluzione con cvx
cvx_clear; %[output:1038c82c]
cvx_begin sdp quiet
    cvx_precision high
    variable Q(n,n) symmetric
    variable L(m,n)
    minimize(0);
    subject to
        Q >= (1e-6)*eye(n);
        A*Q + Q*A' + B*L + L'*B' <= (1e-6)*eye(n);
cvx_end;

K_LMI = -L/Q
%[text] #### 
%[text] #### Compensazione del Guadagno
gain = dcgain(ss((A-B*K_LMI), B, C_sel, D(1:2, 1:2)));
N = pinv(gain)
%[text] 
%[text] #### Verifica Simulazione
%[text] Schema Simulink: [**Verifica\_LMI**](matlab:open_system('Verifica_LMI.slx'))
%[text] 
%[text] 
%%
%[text] ## 3)  CONTROLLO LQR CON AZIONE INTEGRALE
%[text] 
%[text] #### Specifiche
% Riferimento
h2_ref_lqr = 2.8;

% Tempo di Assestamento
Ts = 1;
%[text] 
%[text] #### Sistema Aumentato
C_LQR = [0, 1, 0];

p_a = size(C_LQR, 1);

A_a = [A, zeros(n, p_a); %[output:group:06f8114d] %[output:1a250746]
       -C_LQR, 0] %[output:group:06f8114d] %[output:1a250746]
n_a = size(A_a, 1);

B_a = [B; %[output:group:2f31f91e] %[output:7d045af8]
       zeros(p_a,m)] %[output:group:2f31f91e] %[output:7d045af8]
%[text] 
%[text] #### Sintesi LQR
Q_x = eye(n);
Q_z = 10*eye(p_a);

Q = blkdiag(Q_x, Q_z);

R = eye(m);

% Trasliamo lo spettro di A_a per rispettare la specifica sul Ts
alpha = 4/Ts;
A_hat = A_a + alpha*eye(n_a);

[K_LQR, ~, ~] = lqr(A_hat, B_a, Q, R) %[output:2be88114]
K_LQR_x = K_LQR(:, 1:n) %[output:3ad4ef1c]
K_LQR_i = K_LQR(:, n+1:end) %[output:8db5c3c6]
%[text] 
%[text] #### Verifica Simulazione
%[text] Schema Simulink: [**Verifica\_LQR**](matlab:open_system('Verifica_LQR.slx'))
%[text] 
%[text] 
%%
%[text] ## **4)  CONTROLLO H-INF**
%[text] 
%[text] #### Sistema con Disturbo
% Equazione con disturbo
eq_2 = -A2*dh2 + a12*sqrt(h1-h2) - a23*sqrt(h2-h3) - a2*sqrt(h2) + w == 0;

eqs_w = [eq_1;
         eq_2;
         eq_3];

sol_w = solve(eqs_w, [dh1, dh2, dh3], "ReturnConditions",true);

% dx/dt
f_w = [sol_w.dh1;
       sol_w.dh2;
       sol_w.dh3];

ro = 0.1;

% y
g_w = [x;
      ro*u];

A_w_sym = jacobian(f_w, x);
B_u_sym = jacobian(f_w, u);
B_w_sym = jacobian(f_w, w)
C_w_sym = jacobian(g_w, x);
D_w_sym = jacobian(g_w, u);

A_w = double(subs(A_w_sym, [x; u], [x_star; u_star]));
B_u = double(subs(B_u_sym, [x; u], [x_star; u_star]));
B_w = double(subs(B_w_sym, [x; u], [x_star; u_star]));
C_w = double(subs(C_w_sym, [x; u], [x_star; u_star]));
D_w = double(subs(D_w_sym, [x; u], [x_star; u_star]));

Btot = [B_u B_w];

sys_w = ss(A_w, B_u, C_w, D_w)
B_w
%[text] #### 
%[text] #### Sintesi H-inf
[K_HINF, ~, gamma] = hinfsyn(sys_w, 5, 2);
K_HINF

if gamma>1
    disp('MIXSYN  FAILED:');
else
    fprintf("\n---- MIXSYN RESULT ----\n");
end

fprintf("Achieved gamma = %.4f\n", gamma);

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:59b43b97]
%   data: {"dataType":"text","outputData":{"text":"    <strong>x_star_1<\/strong>    <strong>x_star_2<\/strong>    <strong>x_star_3<\/strong>    <strong>u1_star<\/strong>    <strong>u2_star<\/strong> \n    <strong>________<\/strong>    <strong>________<\/strong>    <strong>________<\/strong>    <strong>_______<\/strong>    <strong>________<\/strong>\n\n     3.274         2          1.6       0.67723    0.063246\n\n","truncated":false}}
%---
%[output:49bbd0e2]
%   data: {"dataType":"text","outputData":{"text":"\nsys =\n \n  A = \n            x1       x2       x3\n   x1  -0.2658   0.2658        0\n   x2   0.2658  -0.6881   0.3162\n   x3        0   0.3162   -0.415\n \n  B = \n       u1  u2\n   x1   1   0\n   x2   0   0\n   x3   0   1\n \n  C = \n       x1  x2  x3\n   y1   1   0   0\n   y2   0   1   0\n   y3   0   0   1\n \n  D = \n       u1  u2\n   y1   0   0\n   y2   0   0\n   y3   0   0\n \nContinuous-time state-space model.\n<a href=\"matlab:disp(char([10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 65 58 32 91 51 215 51 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 66 58 32 91 51 215 50 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 67 58 32 91 51 215 51 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 68 58 32 91 51 215 50 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 69 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 79 102 102 115 101 116 115 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 83 99 97 108 101 100 58 32 48 10 32 32 32 32 32 32 32 32 83 116 97 116 101 78 97 109 101 58 32 123 51 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 32 32 83 116 97 116 101 80 97 116 104 58 32 123 51 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 32 32 83 116 97 116 101 85 110 105 116 58 32 123 51 215 49 32 99 101 108 108 125 10 32 32 32 32 73 110 116 101 114 110 97 108 68 101 108 97 121 58 32 91 48 215 49 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 73 110 112 117 116 68 101 108 97 121 58 32 91 50 215 49 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 79 117 116 112 117 116 68 101 108 97 121 58 32 91 51 215 49 32 100 111 117 98 108 101 93 10 32 32 32 32 32 32 32 32 73 110 112 117 116 78 97 109 101 58 32 123 50 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 32 32 73 110 112 117 116 85 110 105 116 58 32 123 50 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 32 73 110 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 79 117 116 112 117 116 78 97 109 101 58 32 123 51 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 32 79 117 116 112 117 116 85 110 105 116 58 32 123 51 215 49 32 99 101 108 108 125 10 32 32 32 32 32 32 79 117 116 112 117 116 71 114 111 117 112 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10 32 32 32 32 32 32 32 32 32 32 32 32 78 111 116 101 115 58 32 91 48 215 49 32 115 116 114 105 110 103 93 10 32 32 32 32 32 32 32 32 32 85 115 101 114 68 97 116 97 58 32 91 93 10 32 32 32 32 32 32 32 32 32 32 32 32 32 78 97 109 101 58 32 39 39 10 32 32 32 32 32 32 32 32 32 32 32 32 32 32 32 84 115 58 32 48 10 32 32 32 32 32 32 32 32 32 84 105 109 101 85 110 105 116 58 32 39 115 101 99 111 110 100 115 39 10 32 32 32 32 32 83 97 109 112 108 105 110 103 71 114 105 100 58 32 91 49 215 49 32 115 116 114 117 99 116 93 10]))\">Model Properties<\/a>\n","truncated":false}}
%---
%[output:57a57d30]
%   data: {"dataType":"text","outputData":{"text":"Il punto di equiibrio è asintoticamente stabile\n","truncated":false}}
%---
%[output:0a0396e9]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"autovalori_A","rows":3,"type":"double","value":[["-0.9690"],["-0.3390"],["-0.0609"]]}}
%---
%[output:2bda1780]
%   data: {"dataType":"text","outputData":{"text":"Il sistema è completamente raggiungibile\n","truncated":false}}
%---
%[output:732c43f7]
%   data: {"dataType":"text","outputData":{"text":"Il sistema è completamente osservabile\n","truncated":false}}
%---
%[output:1038c82c]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'cvx_clear'."}}
%---
%[output:1a250746]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"A_a","rows":4,"type":"double","value":[["-0.2658","0.2658","0","0"],["0.2658","-0.6881","0.3162","0"],["0","0.3162","-0.4150","0"],["0","-1.0000","0","0"]]}}
%---
%[output:7d045af8]
%   data: {"dataType":"matrix","outputData":{"columns":2,"name":"B_a","rows":4,"type":"double","value":[["1","0"],["0","0"],["0","1"],["0","0"]]}}
%---
%[output:2be88114]
%   data: {"dataType":"matrix","outputData":{"columns":4,"name":"K_LQR","rows":2,"type":"double","value":[["13.5851","251.7205","7.1912","-706.4912"],["7.1912","301.5251","15.9478","-841.6954"]]}}
%---
%[output:3ad4ef1c]
%   data: {"dataType":"matrix","outputData":{"columns":3,"name":"K_LQR_x","rows":2,"type":"double","value":[["13.5851","251.7205","7.1912"],["7.1912","301.5251","15.9478"]]}}
%---
%[output:8db5c3c6]
%   data: {"dataType":"matrix","outputData":{"columns":1,"name":"K_LQR_i","rows":2,"type":"double","value":[["-706.4912"],["-841.6954"]]}}
%---
