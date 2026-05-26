% ===================================
% Script: confronto_norme.m
% ===================================

clear; clc; close all;

% --- Griglia ---
x1_range = -2:0.01:2;
x2_range = -2:0.01:2;
[x1,x2] = meshgrid(x1_range, x2_range);



% --- Calcolo delle norme ---
norm_inf = max(abs(x1), abs(x2));         % Norma infinito
norm_2   = sqrt(x1.^2 + x2.^2);           % Norma 2


% --- Curve di livello ---
levels = 0:0.25:2;   % Valori delle curve

figure;

subplot(1,2,1);
contour(x1, x2, norm_inf, levels, 'ShowText','on','LineWidth',1.5);
axis equal; grid on;
title('Norma infinito: ||x||_\infty');
xlabel('x'); ylabel('y');


subplot(1,2,2);
contour(x1, x2, norm_2, levels, 'ShowText','on','LineWidth',1.5);
axis equal; grid on;
title('Norma 2: ||x||_2');
xlabel('x'); ylabel('y');
