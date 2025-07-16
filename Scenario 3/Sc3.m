clear all 
close all 
clc

%% IPERBOLE DI FUGA - SCENARIO 3 (caso coplanare) nel caso di ottimizzazione delta V

% Costanti
G = 6.672591e-11;              % [Nm^2/kg^2] costante gravitazionale
epsilon = deg2rad(23.45);  % Obliquità eclittica

% Terra
R_T = 6378;                    % [km] raggio medio Terra
mu_T = 398600;                 % [km^3/s^2] parametro gravitazionale Terra
M_T = 5.9741e24;               % [kg] massa Terra

% Sole
mu_S = 132712440018;           % [km^3/s^2]
M_S = 1.9891e30;               % [kg]

% Stato iniziale: orbita di parcheggio (da SC1)
rrf = [-3457.33; 6440.71; 2824.78];            % [km]
vvf = [-8.1090; -4.3270; -0.0599];             % [km/s]

% Velocità eccesso iperbolico (da SC2) nel caso di ottimizzazione deltav
v1T_n = [10.6091; 37.5509; 0.4811];                  % [km/s] trasferimento eliocentrico
vv1i=[-13.4511; 26.7331; -0.0019];
v_inf1_vec1 = (v1T_n-vv1i);                           % v_infinity geocentrica

% Sphere of Influence della Terra
r_soi_T = norm(rrf) * (M_T / M_S)^(2/5);             % [km]

%Applico funzione HypFuga
% esempio che mi funziona con il plot incorporato nella funz HypFuga
% rrf = [7000; 0; 0];                % posizione iniziale [km]
% vvf = [0; 7.546; 0];
% % Velocità eccesso iperbolico (v_inf)
% v_inf1_vec = [3; 0; 0];           % v_inf = 3 km/s, direzione tangenziale

[DeltaV1, v_inf1, v_H1, v_circ, a_H, e_H, Delta_H1] = HypFuga(rrf, vvf, v_inf1_vec1, mu_T, true);



%% === GRAFICO ORBITA DI PARCHEGGIO E IPERBOLE DI FUGA ===
[a, e, i, om, OM, th] = car2par(rrf, vvf, mu_T);
[rP,vP] = par2car(a,e,i,OM,om,0,mu_T);
figure;
hold on; axis equal; grid on;
axis([-8e4 8e4 -8e4 8e4 -2e4 8e4]); % [xmin xmax ymin ymax zmin zmax]
xlabel('x [km]', 'FontSize', 12);
ylabel('y [km]', 'FontSize', 12);
zlabel('z [km]', 'FontSize', 12);
title('Orbita di Parcheggio e Iperbole di Fuga Coplanare', 'FontSize', 14);

% === Orbita di Parcheggio ===
theta_vals = linspace(0, 2*pi, 400);
r_p = a * (1 - e^2) ./ (1 + e * cos(theta_vals));
x_p = r_p .* cos(theta_vals);
y_p = r_p .* sin(theta_vals);
z_p = zeros(size(x_p));
R_p = rotation_matrix(OM, i, om);
orb_parking = R_p * [x_p; y_p; z_p];
plot3(orb_parking(1,:), orb_parking(2,:), orb_parking(3,:), ...
      '-', 'Color', [0 0.4 1], 'LineWidth', 2, 'DisplayName', 'Orbita di parcheggio');

% === Iperbole di Fuga ===
theta_HH = acos(-1/e_H);
theta_vals_H = linspace(-theta_HH + 0.01, theta_HH - 0.01, 500);
r_H_vals = a_H * (1 - e_H^2) ./ (1 + e_H * cos(theta_vals_H));
x_H = r_H_vals .* cos(theta_vals_H);
y_H = r_H_vals .* sin(theta_vals_H);
z_H = zeros(size(x_H));
R_H = rotation_matrix(OM, i, om);
orb_hyper = R_H * [x_H; y_H; z_H];
plot3(orb_hyper(1,:), orb_hyper(2,:), orb_hyper(3,:), ...
      '-', 'Color', [1 0 0], 'LineWidth', 2.5, 'DisplayName', 'Iperbole di fuga');

% === Terra ===
[X,Y,Z] = sphere(50);
surf(6378*X, 6378*Y, 6378*Z, ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
     'FaceColor', [0.5 0.7 1], 'DisplayName', 'Terra');

% === Punto di Intersezione r_H ===
plot3(rP(1), rP(2), rP(3), 'ko', 'MarkerSize', 8, ...
      'MarkerFaceColor', 'k', 'DisplayName', 'Punto r_H');

% === Vista 3D e illuminazione ===
view(35, 25); % Angolo di visualizzazione 3D
camlight headlight; lighting gouraud;

% === Vettore r_H (dal centro alla posizione di intersezione con l'iperbole)
quiver3(0,0,0, rP(1), rP(2), rP(3), ...
        'Color', [1 0.6 0], 'LineWidth', 1.5, 'MaxHeadSize', 0.5, ...
        'DisplayName', 'intersection radius (r_H)');
% === Direzione v_infty (facoltativo) ===
    T_eclip_to_eci = [1 0 0;
        0 cos(epsilon) sin(epsilon);
        0 -sin(epsilon) cos(epsilon)];
    v_inf1_vec = T_eclip_to_eci * v_inf1_vec1;
 quiver3(rP(1), rP(2), rP(3), v_inf1_vec(1), v_inf1_vec(2), v_inf1_vec(3), ...
         300, 'Color', 'm', 'LineWidth', 1.5, 'MaxHeadSize', 1.5, 'DisplayName', 'v_\infty');

%% IPERBOLE DI ARRIVO - SCENARIO 3 (caso coplanare) nel caso di ottimizzazione deltav

% Costanti
G = 6.672591e-11;              % [Nm^2/kg^2] costante gravitazionale

% Sole
mu_S = 132712440018;           % [km^3/s^2]
M_S = 1.9891e30;               % [kg]

% Asteroide
m_a = 1.75181e11;              % [kg]
D27 = 0.55;                    % [km] diametro
R_ast = D27/2;                   % [km] raggio
mu_a = G * m_a / (1000^3);     % [km^3/s^2] convertito da m^3/s^2

% dati dall'orbita eliocentrica asteroide
a_f = 1.066293 *1.496*1e8;
e_f = 0.315212;
i_f = 8.51*pi/180;
OM_f = 205.80*pi/180;
om_f = 172.04*pi/180;
rcs=a_f;

% se vogliamo essere precisi, si usa la posizione esatta sapendo anche il
% th da sc2
% [rra,vva] = par2car(a_f,e_f,i_f,OM_f,om_f,th,mu_S);
% rcs = norm(rra);

Ra_SOI = rcs *(m_a/M_S)^(2/5);  % [km] Sphere of Influence asteroide

% Range pericentro: da impatto diretto a ingresso marginale nella SOI
r_pH = R_ast:0.1:Ra_SOI;   % [km] passo 0.1 km

% Velocità eccesso iperbolico (da SC2)
v2T_n = [10.3523; -20.5084; 1.0517];                  % [km/s] trasferimento eliocentrico
vv2f=[10.2049; -18.0308; 3.0936];
v_inf2_vec = (v2T_n-vv2f);                           % v_infinity geocentrica
v_inf2=norm(v_inf2_vec);

% Preallocazione
DeltaV2 = zeros(size(r_pH));
Delta_H = zeros(size(r_pH));

%% Calcoli
a_H = -mu_a / v_inf2^2;
eH=[];

for k = 1:length(r_pH)
    r_oc = r_pH(k);
    if r_oc < R_ast
    error('Il pericentro interseca la superficie asteroide!');
    end
    
    % Velocità iperbolica al pericentro
    v_pH = sqrt(v_inf2^2 + 2 * mu_a / r_oc);
    
    % Velocità in orbita circolare (ipotetica cattura)
    v_oc = sqrt(mu_a / r_oc);
    
    % Delta-V per inserimento
    DeltaV2(k) = abs(v_pH - v_oc);
    
    % Parametri iperbolici
    e_H = 1 - r_oc / a_H;
    Delta_H(k) = -a_H * sqrt(e_H^2 - 1);  % parametro di impatto
    eH=[eH,e_H];
end

% Plot Delta-V vs r_pH
figure;
plot(r_pH-R_ast, DeltaV2,'Marker','o','Color','r');
xlabel('r_{pH} [km]');
ylabel('\Deltav_2 [km/s]');
title('\Deltav_2 vs rp-Ra');
xline(Ra_SOI-R_ast, '--k', 'SOI');

grid on;

% Plot Delta (imp. param) vs r_pH
figure;
plot(r_pH-R_ast, Delta_H, 'Marker','o','Color','b');
xlabel('r_{pH} [km]');
ylabel('\Delta_H [km]');
title('Parametro di impatto \Delta_H vs rp-Ra');
xline(Ra_SOI-R_ast, '--k', 'SOI');

grid on;


%% Plot di 5 orbite iperboliche

r_pH = linspace(R_ast,Ra_SOI,5);
theta_inf_offset = 0.01;
% Direzione unità lungo v_inf
theta_vinf = atan2(v_inf2_vec(2), v_inf2_vec(1));
R_rot = [cos(theta_vinf), -sin(theta_vinf); sin(theta_vinf), cos(theta_vinf)];

% Plot
figure; hold on; axis equal; grid on;
axis([-0.4e1 1.4e1 -0.6e1 0.6e1 -1e1 1e1]);
title('Asteroid Inertial Reference Frame');
xlabel('X [km]'); ylabel('Y [km]');

% Asteroide (riempito)
theta_c = linspace(0, 2*pi, 200);
fill(R_ast*cos(theta_c), R_ast*sin(theta_c), [0.4 0.4 0.4], 'EdgeColor', 'none');

% Iperboli orientate (interrotte al pericentro)
for i = 1:length(r_pH)
    rp = r_pH(i);
    e_H = 1 - rp / a_H;
    if e_H <= 1, continue; end  % non iperbolica
    theta_inf = acos(-1 / e_H);
 
    % Solo tratto da -θ_inf a θ = 0 (cioè fino al pericentro)
    theta_hyp = linspace(-theta_inf + theta_inf_offset, 0, 400);
    r = a_H * (e_H^2 - 1) ./ (1 + e_H * cos(theta_hyp));
    x = r .* cos(theta_hyp);
    y = r .* sin(theta_hyp);

    % Ruota lungo direzione v_inf2_vec
    coords_rot = R_rot * [x; y];
    plot(coords_rot(1,:), coords_rot(2,:), 'b', 'LineWidth', 1.2);

    % Punto pericentro (manovra)
    plot(0, r_pH, 'ko', 'Marker', 'o');
end

% Orbite circolari di cattura
colors = {'g--', 'b--', 'm--', 'c--', 'r--'};
for i = 1:length(r_pH)
    rc = r_pH(i);
    plot(rc*cos(theta_c), rc*sin(theta_c), colors{i}, 'LineWidth', 1.0);
    yline(Ra_SOI,'--k', 'SOI','LineWidth',2)
end

legend('Asteroide', 'Iperboli di arrivo', 'Manovre', 'Location','best');



%% Fuel fraction vs r_pH (Tsiolkovsky)
r_pH = R_ast:0.1:Ra_SOI;
Isp = 300;                  % [s] impulso specifico
g0 = 9.81e-3;               % [km/s^2] accelerazione gravità standard
fuel_frac = 1 - exp(-DeltaV2 ./ (Isp * g0));   % massa persa in proporzione

%Plot Delta-V + Fuel fraction
% assume un impulso specifico Isp=300 s (valore tipico per motori chimici)

figure;
yyaxis left
plot(r_pH, DeltaV2, 'r', 'LineWidth', 2);
ylabel('\Deltav_2 [km/s]');
xlabel('r_{pH} [km]');
title('\Deltav_2 e Fuel Fraction vs r_{pH}');
grid on;

yyaxis right
plot(r_pH, fuel_frac, 'b--', 'LineWidth', 2);
ylabel('Frazione di massa consumata [-]');
legend('\Deltav_2', 'Frazione carburante');

%% conclusioni 

deltaVTOT=DeltaV1+DeltaV2;

min(deltaVTOT); % 24.3287 km/s
max(deltaVTOT); % 24.3288 km/s
