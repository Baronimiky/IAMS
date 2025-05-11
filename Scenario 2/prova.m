clc; clear; close all;

%% PARAMETRI ORBITALI
mu = 398600; % km^3/s^2, parametro gravitazionale della Terra
r0 = 7000; % km, raggio orbitale iniziale
v0 = sqrt(mu / r0); % km/s, velocità orbitale iniziale

% Tempo di simulazione (1 orbita)
T_orbita = 2 * pi * sqrt(r0^3 / mu);
tspan = linspace(0, T_orbita, 500);

% Condizioni iniziali [x, y, z, vx, vy, vz]
y0 = [r0; 0; 0; 0; v0 * cos(pi/6); v0 * sin(pi/6)]; % Orbita inclinata di 30°

%% EQUAZIONI DEL MOTO ORBITALE (3D)
orbit_eq = @(t, y) [y(4); y(5); y(6); 
                     -mu * y(1) / norm([y(1), y(2), y(3)])^3; 
                     -mu * y(2) / norm([y(1), y(2), y(3)])^3;
                     -mu * y(3) / norm([y(1), y(2), y(3)])^3];

% Risoluzione numerica dell'orbita
[t, Y] = ode45(orbit_eq, tspan, y0);

%% CREAZIONE DELLA FIGURA E DELLA TERRA
figure;
hold on;
[Xs, Ys, Zs] = sphere(50); % Genera una sfera per la Terra
surf(Xs * 6378, Ys * 6378, Zs * 6378, 'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
axis equal;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
grid on;
title('Orbita Satellitare con Manovra');
view(3);

%% LINEA ANIMATA E SATELLITE COME PALLINO
h = animatedline('Color', 'r', 'LineWidth', 2);
satellite = plot3(Y(1,1), Y(1,2), Y(1,3), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8); % Pallino rosso

%% ANIMAZIONE DELL'ORBITA INIZIALE
for i = 1:length(Y)
    addpoints(h, Y(i,1), Y(i,2), Y(i,3));
    set(satellite, 'XData', Y(i,1), 'YData', Y(i,2), 'ZData', Y(i,3)); % Aggiorna la posizione del satellite
    drawnow;
end

%% MANOVRA ORBITALE (Aumento di Velocità del 10%)
disp('ESECUZIONE DELLA MANOVRA...');
pause(1); % Simula il tempo necessario alla manovra

% Punto della manovra (a metà orbita)
idx_manoeuvre = round(length(t) / 2);
Y(idx_manoeuvre, 4:6) = Y(idx_manoeuvre, 4:6) * 1.1; % Incremento del 10% della velocità

% Nuova simulazione post-manovra
[t2, Y2] = ode45(orbit_eq, tspan, Y(idx_manoeuvre, :));

%% AGGIORNAMENTO DEL PLOT CON LA NUOVA ORBITA
h2 = animatedline('Color', 'g', 'LineWidth', 2);
for i = 1:length(Y2)
    addpoints(h2, Y2(i,1), Y2(i,2), Y2(i,3));
    set(satellite, 'XData', Y2(i,1), 'YData', Y2(i,2), 'ZData', Y2(i,3)); % Aggiorna il pallino
    drawnow;
end

legend('Terra', 'Orbita Iniziale', 'Orbita dopo Manovra', 'Satellite');
disp('MANOVRA COMPLETATA!');
