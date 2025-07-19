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
v1T_n = [-13.6985; 29.1856; -1.4761];                  % [km/s] trasferimento eliocentrico
vv1i=[-13.4511; 26.7331; -0.0019];
v_inf_eclip = (v1T_n-vv1i);                           % v_infinity geocentrica

DistT_S=149.61e6;              % [km]
% % Sphere of Influence della Terra
r_soi_T = DistT_S * (M_T / M_S)^(2/5);             % [km]

% 1. Trasformazione da ECLIP (ECLITTICO) a ECI (PLANETOCENTRICO)
r_H=rP; %raggio di intersezione tra orbita di parcheggio e hyp di fuga = raggio pericentro OP
T_eclip_to_eci = [1 0 0;
    0 cos(epsilon) sin(epsilon);
    0 -sin(epsilon) cos(epsilon)];
v_inf = T_eclip_to_eci * v_inf_eclip;

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
%ehvec=[ehvec,e_H];

%% verifica vincoli 6. Calcolo θ_H e θ_∞
%Verifica vincoli
% Se i vincoli non sono rispettati: non esiste un’iperbole in grado di soddisfare 
% sia la 𝒗∞ desiderata che il raggio di intersezione 𝒓𝑯  scelto

theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
theta_inf = alpha + theta_H;
RPH=a_H*(1-e_H)

   
    % if theta_inf <= pi/2 || theta_inf >= pi
    %     error('Vincolo violato: theta_inf non è compreso tra π/2 e π');
    % end
    % 
    % if RPH <= R_T
    %     error('Vincolo violato: r_PH è troppo basso (sotto il perigeo minimo consentito)');
    % end



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

    %f = @(e_H) acos(-1/e_H) - alpha - acos((a_H*(1 - e_H.^2)/norm(r_H) - 1) / e_H);

    %fplot(f)
    %yline(0)
    % Risolvo numericamente con fzero
    % 
    % e_0 = 2;  % stima iniziale
    % e_H = fsolve(f,e_0);
    % ehvec=[ehvec,e_H];
    safe_acos = @(x) acos(min(max(x, -1), 1));
f = @(e_H) safe_acos(-1/e_H) - alpha - safe_acos((a_H*(1 - e_H.^2)/norm(r_H) - 1) / e_H);

% Ricerca nell’intervallo [1.01, 20]
try
    e_H = fzero(f, [1.01, 20]);
    if isreal(e_H)
        ehvec = [ehvec, e_H];
    else
        ehvec = [ehvec, NaN];  % oppure scarta
    end
catch
    ehvec = [ehvec, NaN];  % soluzione non trovata
end

end
%% 6.b
clc
thInf=[];
rRPH=[];
ehvalidi=[];

% Verifica vincoli
for j=1:360
    [r_H,voP] = par2car(a,e,i,OM,om,deg2rad(j),mu_T);
    norm(r_H)
    e_H=ehvec(j);
    theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
    theta_inf = alpha + theta_H;
    % theta_inf=acos(-1/e_H);
    % theta_H=theta_inf-alpha;
    RPH=a_H*(1-e_H);
    
    
    if RPH>R_T && pi/2 < theta_inf < pi
        ehvalidi=[ehvalidi,e_H];
        thInf=[thInf,theta_inf];
        rRPH=[rRPH,RPH];
    else
        ehvalidi=[ehvalidi,NaN];
        thInf=[thInf,NaN];
        rRPH=[rRPH,NaN];
    end
end
figure
plot(1:360,ehvec)

%ottengo eH validi per  58°<theta<67 e 241<teta<271

%% andamento delta V1

DeltaV1= zeros(size(360));
for j=1:360
    [r_H,voP] = par2car(a,e,i,OM,om,deg2rad(j),mu_T);
    e_H=ehvalidi(j);
    theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
    theta_inf = alpha + theta_H;
    % theta_inf=acos(-1/e_H);
    % theta_H=theta_inf-alpha;
    RPH=a_H*(1-e_H);
% 7. Momento angolare, Nodi, Inclinazione, RAAN

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

% SLIDE 10: Calcolo ω_H tramite angolo β

% Angolo β tra linea dei nodi e r_H
cos_beta = dot(N_hat, r_H / norm(r_H));
sin_beta = dot(cross(N_hat, r_H / norm(r_H)), h_hat);
beta = atan2(sin_beta, cos_beta);

% Argomento del pericentro
omega_H = mod(beta - theta_H, 2*pi);  % normalizzato in [0, 2π]

[r_p, V_H] = par2car(a_H, e_H, i_H, Omega_H, omega_H, theta_H, mu_T);
DeltaV1(j)=abs(norm(V_H)-norm(voP));
end

figure()
plot(1:360,DeltaV1);
% Plot Delta-V vs theta
figure;
plot(1:360,DeltaV1,'Marker','o','Color','r');
xlabel('\theta [km]');
ylabel('\DeltaV_1 [km/s]');
title('\DeltaV_2 vs \thetaOP');



[deltav1,j]=min(DeltaV1); %per theta = 58°

%% scelta theta

[r_H,voP] = par2car(a,e,i,OM,om,deg2rad(j),mu_T);
    e_H=ehvalidi(j);
    theta_H = acos((a_H*(1 - e_H^2)/norm(r_H) - 1) / e_H);
    theta_inf = alpha + theta_H;
    % theta_inf=acos(-1/e_H);
    % theta_H=theta_inf-alpha;
    RPH=a_H*(1-e_H);
% 7. Momento angolare, Nodi, Inclinazione, RAAN

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

% SLIDE 10: Calcolo ω_H tramite angolo β

% Angolo β tra linea dei nodi e r_H
cos_beta = dot(N_hat, r_H / norm(r_H));
sin_beta = dot(cross(N_hat, r_H / norm(r_H)), h_hat);
beta = atan2(sin_beta, cos_beta);

% Argomento del pericentro
omega_H = mod(beta - theta_H, 2*pi);  % normalizzato in [0, 2π]

[r_p, V_H] = par2car(a_H, e_H, i_H, Omega_H, omega_H, theta_H, mu_T);
DeltaV1=abs(norm(V_H)-norm(voP))



%%
% === GRAFICO ORBITA DI PARCHEGGIO E IPERBOLE DI FUGA ===

figure;
hold on; axis equal; grid on;
axis([-5e4 6e4 -5e4 7e4 -2e4 2e4]); % [xmin xmax ymin ymax zmin zmax]
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
     'FaceColor', [0.5 0.7 1], 'DisplayName', 'Terra','HandleVisibility', 'off');

% === Punto di Intersezione r_H ===
plot3(r_H(1), r_H(2), r_H(3), 'ko', 'MarkerSize', 4, ...
      'MarkerFaceColor', 'k', 'DisplayName', 'Punto r_H','HandleVisibility', 'off');

r_end = orb_post(:, 208);  
vsc=4;
% === Direzione v_infty (facoltativo) ===
 quiver3(r_end(1), r_end(2), r_end(3), vsc*v_inf(1), vsc*v_inf(2), vsc*v_inf(3),500, 'Color', 'k', 'LineWidth', 3, 'MaxHeadSize', 3, 'DisplayName', 'v_\infty');

% === Vista 3D e illuminazione ===
view(30, 25); % Angolo di visualizzazione 3D
camlight headlight; lighting gouraud;

% === Vettore r_H (dal centro alla posizione di intersezione con l'iperbole)
quiver3(0,0,0, r_H(1), r_H(2), r_H(3), ...
        'Color', [1 0.6 0], 'LineWidth', 1.5, 'MaxHeadSize', 0.5, ...
        'DisplayName', 'Raggio di intesezione (r_H)');

% Calcolo versore eccentricità se non già noto
[r_pH_vec, V_H] = par2car(a_H, e_H, i_H, Omega_H, omega_H, 0, mu_T);

%r_pH_vec = (RPH) * e_vec / norm(e_vec);  % vettore posizione del pericentro iperbolico

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
view(55,19)
%% calcolo delta T1
[theta_exit, delta_T1] = TOFHyp(a_H, e_H, mu_T, 0,r_soi_T) %2.7821e5 s tempo per uscire dal soi terrestre
deltat = TOF(a,e,th,deg2rad(j),mu_T); %9.9771e2s
DetlT1=deltat+delta_T1; %2.79201e5


%% conclusioni

%deltav2 trovato nel caso coplanare
DeltaV2=3.2139;
deltaVTOT=DeltaV1+DeltaV2;%4.660
norm(r_H)
Delta_H1 = -a_H * sqrt(e_H^2 - 1);% 2.587115616448778e+04
