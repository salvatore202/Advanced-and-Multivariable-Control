clc; 
clear; 
close all;

k0 = 1;
k1 = 1;
b = 1;
M = 1;

% Sistema non lineare
odefun = @(t, x)[x(2); -k0/M*x(1) - k1/M*(x(1))^3 - b/M*x(2)*abs(x(2))];

tspan = [0 50];
x0 = [1; 0]; % Posizione e velocità iniziale
[t, x] = ode45(odefun, tspan, x0);

% Plot risultati
figure;
plot(t, x(:, 1), 'b-', 'LineWidth', 2);
hold on;
plot(t, x(:, 2), 'r-', 'LineWidth', 2);
xlabel('Tempo (s)');
ylabel('Spostamento (m) / Velocita (m/s)');
title('Risposta del sistema non lineare');
legend('Spostamento','Velocità');

% Plot pp
figure;
plot(x(:, 1), x(:, 2), 'k-', 'LineWidth', 2);
xlabel('Spostamento (m)');
ylabel('Velocita (m/s)');
title('Ritratto di fase del sistema non lineare');

% Linearizzazione attorno all'equilibrio
x1_bar = 0;   % posizione di equilibrio
x2_bar = 0;   

% Matrice Jacobiana (lineare)
J = [ 0                     1;
     -k0/M             -b/M*0];
B = [0; 0];
C = eye(2);
D = 0;

disp('Matrice J')
disp(J)

% Calcolo degli autovalori e autovettori
[eigenVectors, eigenValues] = eig(J);
disp('Autovalori:')
disp(diag(eigenValues))
disp('Autovettori:')
disp(eigenVectors)

% % Stabilità del sistema
% stability = all(real(diag(eigenValues)) < 0);
% if stability
%     disp('Il sistema è stabile.');
% else
%     disp('Il sistema è instabile.');
% end

%Approccio alla Lyapunov
[X1, X2] = meshgrid(-1:.01:1, -1:.01:1);
V = 0.5*M*X2.^2 + 0.5*k0*X1.^2 + 1/4*k1*X1.^4;  %X1. ogni elemento di X1
figure
mesh(X1, X2, V) %crea grafico 3D

V_dot = -b*abs(X2.^3);
% Visualizzazione della Lyapunov funzione e della sua derivata
figure;
mesh(X1, X2, V_dot) %crea grafico 3D