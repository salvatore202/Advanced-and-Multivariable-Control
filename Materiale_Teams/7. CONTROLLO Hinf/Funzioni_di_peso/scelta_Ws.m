close all
clear all
clc

s = tf('s');

%% Ws(s)

% Specifiche
e_max = 0.01;   % errore massimo 1%
Ts = 2;         % tempo di assestamento
PO = 10;        % sovraelongazione massima %

% Derivazioni
A = e_max;
wc = 4/Ts;

zeta = 0.6;                   % da tabella approx per PO = 10%
Ms = 1/(2*zeta);              % picco massimo
M =Ms;

Ws = (s/M + wc) / (s + A*wc);
figure
bode(Ws,{1e-5,1e5})
title('Ws')

%% Wt(s)

omega_t = 30;    % frequanza minima del rumore
Mt = 100;        % desidero ridurre T di 40 dB in alta freq
At = 1e-4;       % modulo in bassa frequenza

Wt = ( s + omega_t*At ) / ( s/Mt + omega_t );

figure
bode(Wt,{1e-5,1e5})
title('Wt')