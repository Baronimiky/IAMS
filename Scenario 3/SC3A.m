
clear all 
close all
clc

%% IPERBOLE DI FUGA - SCENARIO 3 (caso Fuori Piano) nel caso di ottimizzazione delta V
clc
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

[a,e,i,OM, om, th] = car2par(rrf,vvf,mu_T);

[rP,vP] = par2car(a,e,i,OM,om,0,mu_T); %Raggio e velocità al pericentro dell'orbita di parcheggio

% Velocità eccesso iperbolico (da SC2) nel caso di ottimizzazione deltav
v1T_n = [10.6091; 37.5509; 0.4811];                  % [km/s] trasferimento eliocentrico
vv1i=[-13.4511; 26.7331; -0.0019];
v_inf_eclip = (v1T_n-vv1i);                           % v_infinity geocentrica

Dist=149.61e6;
% % Sphere of Influence della Terra
r_soi_T = Dist * (M_T / M_S)^(2/5);             % [km]

%% 5.a Cambio rH 
ehvec=[];
for j=1:360

    [rP,vP] = par2car(a,e,i,OM,om,deg2rad(j),mu_T);
    % 1. Trasformazione da ECLIP (ECLITTICO) a ECI (PLANETOCENTRICO)

    r_H=rP; %raggio di intersezione tra orbita di parcheggio e hyp di fuga = raggio pericentro OP
    T_eclip_to_eci = [1 0 0;
        0 cos(epsilon) sin(epsilon);
        0 -sin(epsilon) cos(epsilon)];
    v_inf = T_eclip_to_eci' * v_inf_eclip;

    % 2. Versori
    r_ver_inf = v_inf / norm(v_inf); %Direzione dell'asintoto uscente dell'iperbole
    r_ver_H = r_H / norm(r_H);

    % 3. Calcolo dell'angolo α tra i due versori
    alpha = acos(dot(r_ver_H,r_ver_inf));

    % 4. Semiasse maggiore a_H
    v_inf_mod = norm(v_inf);
    a_H = -mu_T / v_inf_mod^2;

    % 5. Risoluzione sistema non lineare per e_H
    % Equazione implicita: cos(θ∞) = -1 / e_H ; θ∞ = α + θ_H
    % Funzione implicita in e_H
    %f=@(e_H) acos(a_H*(1-e_H^2)/(norm(r_H)*e_H)-1/e_H)-acos(-1/e_H)+alpha;

    f = @(e_H) acos(-1/e_H) - alpha - acos((a_H*(1 - e_H.^2)/norm(r_H) - 1) / e_H);

    %fplot(f)
    %yline(0)
    % Risolvo numericamente con fzero

    e_0 = 2;  % stima iniziale
    e_H = fsolve(f,e_0);
    ehvec=[ehvec,e_H];
end
figure
plot(1:360,ehvec)
xline(257,Label='Inizio valori non complessi')
xline(320,Label='Fine valori non complessi')
xlabel('Theta orbita parcheggio', 'FontSize', 12);
ylabel('eccentricità iperbole', 'FontSize', 12);
%% verifica vincoli 6. Calcolo θ_H e θ_∞
%Verifica vincoli
% Se i vincoli non sono rispettati: non esiste un’iperbole in grado di soddisfare 
% sia la 𝒗∞ desiderata che il raggio di intersezione 𝒓𝑯  scelto

eHvalidi=ehvec(257:320);
[e_H,rot]=min(eHvalidi); %eH minimo = 2.9740 per una rotazione di 262° di theta

[roP,voP] = par2car(a,e,i,OM,om,deg2rad(262),mu_T);
thInf=[];
rRPH=[];
for j=257:320
    [r_H,voP] = par2car(a,e,i,OM,om,deg2rad(j),mu_T);
    e_H=ehvec(j);
    theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
    theta_inf = alpha + theta_H;
    RPH=a_H*(1-e_H);
    thInf=[thInf,theta_inf];
    rRPH=[rRPH,RPH];
end
    % Verifica vincoli
    if theta_inf <= pi/2 || theta_inf >= pi
        error('Vincolo violato: theta_inf non è compreso tra π/2 e π');
    end

    if RPH <= R_T
        error('Vincolo violato: r_PH è troppo basso (sotto il perigeo minimo consentito)');
    end
%applicando i vari angoli l'eH adatta va dai 12.1627 a 13.4152, corrisponde
%a anomalia teta di 308°

%% scelta ufficiale e_H 12.1627
e_H=(eHvalidi(52));
%e_H=7.0753;
[roP,voP] = par2car(a,e,i,OM,om,deg2rad(308),mu_T);
r_H=roP;
theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
theta_inf = alpha + theta_H;
RPH=a_H*(1-e_H);


%% 7. Momento angolare, Nodi, Inclinazione, RAAN

% 1. Momento angolare specifico
h_H = cross(r_H, v_inf);
h_hat = h_H / norm(h_H);

% 2. Inclinazione
i_H = acos(h_hat(3));  % componente z

% 3. Nodo ascendente
K = [0 0 1]';  % versore asse z
N_H = cross(K, h_H);
N_hat = N_H / norm(N_H);

% 4. RAAN Ω_H
i_hat = [1 0 0]';
if N_H(2) >= 0
    Omega_H = acos(dot(N_hat, i_hat));
else
    Omega_H = 2*pi - acos(dot(N_hat, i_hat));
end

%5. vettore eccentricità
e_vec = (cross(v_inf, h_H) / mu_T) - (r_H / norm(r_H));
e_hat = e_vec / norm(e_vec);

%% SLIDE 10: Calcolo ω_H tramite angolo β

% Angolo β tra linea dei nodi e r_H
cos_beta = dot(N_hat, r_H / norm(r_H));
sin_beta = dot(cross(N_hat, r_H / norm(r_H)), h_hat);
beta = atan2(sin_beta, cos_beta);

% Argomento del pericentro
omega_H = mod(beta - theta_H, 2*pi);  % normalizzato in [0, 2π]

%% Output parametri orbitali dell’iperbole
fprintf('\n=== Parametri Orbitali Iperbole ===\n');
fprintf('e_H = %.4f\n', e_H);
fprintf('a_H = %.2f km\n', a_H);
fprintf('Inclinazione i_H = %.2f deg\n', rad2deg(i_H));
fprintf('RAAN Ω_H = %.2f deg\n', rad2deg(Omega_H));
fprintf('Argomento pericentro ω_H = %.2f deg\n', rad2deg(omega_H));

%% vettore velocità corrispondente a RH
[r_p, V_H] = par2car(a_H, e_H, i_H, Omega_H, omega_H, theta_H, mu_T);
DeltaV1=abs(norm(V_H)-norm(voP))

%%
% === GRAFICO ORBITA DI PARCHEGGIO E IPERBOLE DI FUGA ===

figure;
hold on; axis equal; grid on;
axis([-5e4 5e4 -5e4 5e4 -2e4 8e4]); % [xmin xmax ymin ymax zmin zmax]
xlabel('x [km]', 'FontSize', 12);
ylabel('y [km]', 'FontSize', 12);
zlabel('z [km]', 'FontSize', 12);
title('Orbita di Parcheggio e Iperbole di Fuga', 'FontSize', 14);

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
R_H = rotation_matrix(Omega_H, i_H, omega_H);
% orb_hyper = R_H * [x_H; y_H; z_H];
% plot3(orb_hyper(1,:), orb_hyper(2,:), orb_hyper(3,:), ...
%       '-', 'Color', [1 0 0], 'LineWidth', 2.5, 'DisplayName', 'Iperbole di fuga');

% Angoli prima e dopo il punto di intersezione
theta_pre  = linspace(-theta_HH + 0.01, theta_H, 300);     % parte non percorsa
theta_post = linspace(theta_H, theta_HH - 0.01, 300);      % parte percorsa

% === Parte non percorsa (tratteggiata) ===
r_pre  = a_H * (1 - e_H^2) ./ (1 + e_H * cos(theta_pre));
x_pre  = r_pre .* cos(theta_pre);
y_pre  = r_pre .* sin(theta_pre);
z_pre  = zeros(size(x_pre));
orb_pre = R_H * [x_pre; y_pre; z_pre];

plot3(orb_pre(1,:), orb_pre(2,:), orb_pre(3,:), ...
      'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% === Parte percorsa (continua) ===
r_post  = a_H * (1 - e_H^2) ./ (1 + e_H * cos(theta_post));
x_post  = r_post .* cos(theta_post);
y_post  = r_post .* sin(theta_post);
z_post  = zeros(size(x_post));
orb_post = R_H * [x_post; y_post; z_post];

plot3(orb_post(1,:), orb_post(2,:), orb_post(3,:), ...
      'r-', 'LineWidth', 2.5, 'DisplayName', 'Iperbole percorsa');


% === Terra ===
[X,Y,Z] = sphere(50);
surf(6378*X, 6378*Y, 6378*Z, ...
     'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
     'FaceColor', [0.5 0.7 1], 'DisplayName', 'Terra');

% === Punto di Intersezione r_H ===
plot3(r_H(1), r_H(2), r_H(3), 'ko', 'MarkerSize', 8, ...
      'MarkerFaceColor', 'k', 'DisplayName', 'Punto r_H');

% === Direzione v_infty (facoltativo) ===
 quiver3(r_H(1), r_H(2), r_H(3), v_inf(1), v_inf(2), v_inf(3), ...
         300, 'Color', 'm', 'LineWidth', 1.5, 'MaxHeadSize', 1.5, 'DisplayName', 'v_\infty');

% === Vista 3D e illuminazione ===
view(30, 25); % Angolo di visualizzazione 3D
camlight headlight; lighting gouraud;

% === Vettore r_H (dal centro alla posizione di intersezione con l'iperbole)
quiver3(0,0,0, r_H(1), r_H(2), r_H(3), ...
        'Color', [1 0.6 0], 'LineWidth', 1.5, 'MaxHeadSize', 0.5, ...
        'DisplayName', 'Raggio di intesezione (r_H)');

% Calcolo versore eccentricità se non già noto

r_pH_vec = (RPH+100) * e_vec / norm(e_vec);  % vettore posizione del pericentro iperbolico

% === Vettore pericentro iperbolico
quiver3(0,0,0, r_pH_vec(1), r_pH_vec(2), r_pH_vec(3), ...
        'Color', [0.2 0.8 0.2], 'LineWidth', 1.5, 'MaxHeadSize', 0.5, ...
        'DisplayName', 'Raggio Pericentro Iperbole (r_{pH})');

% === Legenda e formattazione ===
legend('Location', 'northeastoutside', 'FontSize', 10);
% === Sistema di riferimento ECI ===
quiver3(0,0,0, 1e4,0,0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off'); % asse X
quiver3(0,0,0, 0,1e4,0, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off'); % asse Y
quiver3(0,0,0, 0,0,1e4, 'k', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off'); % asse Z
text(1.05e4,0,0,'X','FontSize',11);
text(0,1.05e4,0,'Y','FontSize',11);
text(0,0,1.05e4,'Z','FontSize',11);
%%
clc
th2_ott = 0;
for theta_vals = 0:0.1:2*pi
    r_p = a_H * (1 - e_H^2) ./ (1 + e_H * cos(theta_vals));
    x_p = r_p .* cos(theta_vals);
    y_p = r_p .* sin(theta_vals);
    z_p = 0;
    R_p = rotation_matrix(Omega_H, i_H, omega_H);
    orb_hyper = R_p * [x_p; y_p; z_p];
    norm(orb_hyper)
    if norm(orb_hyper)>= r_soi_T && th2_ott == 0 
            th2_ott = theta_vals
    end
end

%theta_vals = linspace(0, 2*pi, 400);
% r_p = a * (1 - e^2) ./ (1 + e * cos(theta_vals));
% x_p = r_p .* cos(theta_vals);
% y_p = r_p .* sin(theta_vals);
% z_p = zeros(size(x_p));
% R_p = rotation_matrix(OM, i, om);



%% conclusioni
%deltav2 trovato nel caso coplanare
DeltaV2=[3.21376088307636	3.21379049603695	3.21381017553883	3.21382446625499	3.21383544989177	3.21384423229354	3.21385146242788	3.21385754954886	3.21386276621076	3.21386730187685	3.21387129288641	3.21387484018972	3.21387802034512	3.21388089260847	3.21388350365658	3.21388589082463	3.21388808438016	3.21389010915518	3.21389198574004	3.21389373137175	3.21389536060493	3.21389688582572	3.21389831765023	3.21389966523702	3.21390093653471	3.21390213847998	3.21390327715747	3.21390435792973	3.21390538554383	3.21390636421927	3.21390729772107	3.21390818942082	3.21390904234800	3.21390985923339	3.21391064254586	3.21391139452382];
deltaVTOT=DeltaV1+DeltaV2;
format long
min(deltaVTOT) % 22.881674933325520 km/s
max(deltaVTOT) % 22.881825444772979 km/s


