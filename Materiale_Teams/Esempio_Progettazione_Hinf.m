%% =======================================
% H∞ MIMO CONTROL EXAMPLE (mixsyn)
% Mixed sensitivity design: Ws (tracking), Wt (noise), Wu (control)
% =======================================

clear; close all; clc;
s = tf('s');

%% -----------------------------
% Tracking specs
e_ss_max = 0.02;     % max steady-state error (1%)
Ts_des   = 0.3;      % desired settling time
PO_max   = 50;       % max overshoot (%)

% Noise rejection specs
noise_atten_dB = 100; % desired attenuation in HF (dB)
omega_noise    = 30; % frequency from which to attenuate

% Control effort penalization
Wu_gain = 0.01;

%% -----------------------------
% MIMO Plant (2x2 example)
G = [  1/(s+1)     0.35/(s+2) ;
       0.25/(s+3)  1.1 /(s+4) ];
[ny, nu] = size(G);

%% =======================================
% 1) DESIGN OF Ws(s) — TRACKING
% =======================================

% Crossover frequency estimate from settling time
wc = 4 / Ts_des;

% Overshoot -> damping ratio -> Ms approx
PO = PO_max / 100;
if PO == 0
    zeta = 1;
else
    lnPO = log(PO);
    zeta = -lnPO / sqrt(pi^2 + lnPO^2);
    if ~(zeta>0 && zeta<1)
        zeta = 0.6;
    end
end

M_s = 1/(2*zeta);   % peak sensitivity approx
As  = e_ss_max;     % Ws(0) = 1/As → enforce steady-state error

% Definition of Ws:

Ws = ( s/M_s + wc ) / ( s + wc*As );

fprintf("\n---- Ws DESIGN ----\n");
fprintf("Ws(0)    = %.3f   (target 1/e_ss_max = %.3f)\n", evalfr(Ws,0), 1/e_ss_max);
fprintf("Ws(inf) ≈ %.3f    (1/Ms = %.3f, Ms = %.3f)\n", 1/M_s, 1/M_s, M_s);
fprintf("wc = %.3f rad/s\n", wc);

%% =======================================
% 2) DESIGN OF Wt(s) — NOISE REJECTION
% =======================================

% Convert attenuation in dB → linear gain
alpha = 10^(noise_atten_dB/20);

omega_t = omega_noise;
Mt = alpha;
At = 1e-4;

Wt = ( s + omega_t*At ) / ( s/Mt + omega_t );

fprintf("\n---- Wt DESIGN ----\n");
fprintf("Noise attenuation target = %.1f dB  -> Mt = %.1f\n", noise_atten_dB, Mt);
fprintf("Wt(inf) = %.3f\n", 1/Mt);

%% =======================================
% 3) CONTROL-EFFORT WEIGHT Wu(s)
% =======================================

Wu = Wu_gain * eye(nu);
fprintf("\n---- Wu DESIGN ----\n");
fprintf("Wu gain = %.3f\n", Wu_gain);

%% =======================================
% PLOT WEIGHTS
figure('Name','Weights Ws and Wt');
bode(Ws, {1e-5,1e5}); hold on;
bode(Wt, {1e-5,1e5});
grid on; legend('Ws','Wt'); title('Weights');

%% =======================================
% 4) H∞ SYNTHESIS (mixsyn)
% =======================================



[K, CL, gamma] = mixsyn(G, Ws*eye(ny), Wu, Wt*eye(ny));
if gamma>1
    disp('MIXSYN  FAILED:');
else
    fprintf("\n---- MIXSYN RESULT ----\n");
end

fprintf("Achieved gamma = %.4f\n", gamma);


%% =======================================
% 4) H∞ SYNTHESIS alternativo (augw e hinfsyn)
% =======================================


% P è strutturato come:
% [ z ]   = P11 [ w ] + P12 [ u ]
% [ y ]         P21 [ w ] + P22 [ u ]

P = augw(G, Ws*eye(ny), Wu, Wt*eye(ny));
[K,CL,gamma] = hinfsyn(P, 2, 2);
if gamma>1
    disp('MIXSYN  FAILED:');
else
    fprintf("\n---- MIXSYN RESULT ----\n");
end

fprintf("Achieved gamma = %.4f\n", gamma);


%% =======================================
% 5) EXTRACT LOOP FUNCTIONS (S, T)
% =======================================
K=tf(K);
L = G * K;
I_n = eye(ny);

S = feedback(I_n, L);   % S = (I + GK)^(-1)
T = I_n - S;            % T = GK (I+GK)^(-1)

%% =======================================
% 6) H∞ NORMS CHECK (should be < 1)
% =======================================

[hn_WsS, ~] = hinfnorm(Ws * S);
[hn_WtT, ~] = hinfnorm(Wt * T);
[hn_WuK, ~] = hinfnorm(Wu*K*S);

fprintf("\n---- H∞ Norms ----\n");
fprintf("||Ws*S||∞ = %.4f\n", hn_WsS);
fprintf("||Wt*T||∞ = %.4f\n", hn_WtT);
fprintf("||Wu*K*S||∞ = %.4f\n", hn_WuK);

%% =======================================
% 7) FREQUENCY DOMAIN PLOTS
% =======================================

w = logspace(-2,5,800);

figure('Name','Sensitivity S and Complementary Sensitivity T');
sigma(S, w); hold on; sigma(T, w);
grid on; legend('S','T'); title('S & T');

figure('Name','Weighted channels');
sigma(Ws*S, w); hold on;
sigma(Wt*T, w);
grid on; legend('Ws*S','Wt*T'); title('Weighted Channels');
