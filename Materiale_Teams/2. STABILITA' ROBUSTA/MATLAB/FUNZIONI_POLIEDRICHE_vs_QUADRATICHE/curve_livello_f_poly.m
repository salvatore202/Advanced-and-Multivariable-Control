% ================================
% Script: curve_livello_poliedrica.m
% ================================

clear; clc; close all;

% --- Definizione dei piani lineari ---
% Ogni riga di A è [a1, a2], B è il termine noto
% Esempio: f(x,y) = max{ x, y, -x-y }
P = [ 1  0;    % piano 1:  x
      0  1;    % piano 2:  y
     -1 -1 ];  % piano 3: -x - y


% --- Griglia su cui valutare la funzione ---
xrange = -2:0.05:2;
yrange = -2:0.05:2;
[x,y] = meshgrid(xrange, yrange);

% --- Calcolo della funzione poliedrica ---
f_vals = -inf(size(x));
for i = 1:size(P,1)
    f_i = abs(P(i,1)*x + P(i,2)*y);
    f_vals = max((f_vals), f_i);
end

% --- Curve di livello ---
figure;
contour(x, y, f_vals, 'ShowText','on','LineWidth',1.5);
axis equal;
grid on;
xlabel('x'); ylabel('y');
title('Curve di livello di una funzione poliedrica');
