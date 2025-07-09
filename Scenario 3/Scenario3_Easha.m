clear all
close all
clc
%%
% DATI
G = 6.67259*1e-11; % N*m^2/Kg^2

% Terra
R_T = 6378; % km    
mu_T = 398600; % km3/s2
M_T =  5.974*1e+24; % kg

% Sole 
mu_S = 132712440018; % km^3/s^2
M_S = 1.989*1e+30; % kg

% Asteroide
m_a = 1.7518*1e11; % Kg
D27 = 0.55;       % Km
R27= D27/2;         % Km
v_inf2 = 3.0295;   % km/s
mu_a = G*m_a/(1000^3);

% Fuga / SC1
rrf = [-12662.895000;-5218.469500;2462.179800]; %[x,y,z] km
vvf = [1.364000;-5.338000;-4.296000];           %[x,y,z] km/s

% Arrivo / SC2
% a27=1.066293;     % AU
% e27=0.315212;
% i27=8.51;         % deg
% OMf=205.80;       % deg
% omf=172.04;       % deg

v1T_n = [10.6091;37.5509;0.4811];
v2T_n = [-36.2068;-28.8691;-1.2861];
vinf1 = norm(v1T_n);
vinf2 = norm(v2T_n);
r_soiT = norm(rrf)*(M_T/M_S)^(2/5);


%% 1 - Iperbole di fuga
%  Raggio di pericentro uguale al raggio di pericentro dell orbita di
%  parcheggio (?)

% 1. Parametri orbitali orbita di parcheggio
[a_f, e_f, i_f, om_f, OM_f, th_f] = car2par(rrf, vvf, mu_T);

% 2. Calcolo raggio pericentro orbita di parcheggio / iperbole
rpf = a_f*(1 + e_f);
r_of = rpf;
v_of = sqrt(mu_T/r_of);

% 3. Calcolo velocità al pericentro dell'iperbole
v_fuga1 = sqrt(2*mu_a/r_of);
v_PH1 = sqrt(vinf1^2+v_fuga1^2);

% 4. Calcolo parametri iperbole
a_H1 = -mu_T/vinf1^2;
e_H1 = 1 - r_of/a_H1;
Delta_n1 = -a_H1*sqrt(e_H1^2 - 1);
DeltaV = abs(v_of - v_PH1)



%% 2 - Iperbole di arrivo
% SETUP
v_oc_n =[];
v_PH2_n =[];
Delta_n_n =[];
Delta_v_n =[];

%  1. Calcolo del semiasse maggiore dell’iperbole
a_H = -mu_a/vinf2^2;
if a_H > 0
    error('Semiasse maggiore positivo');
end

% 2. Definizione di un range di valori per il raggio di pericentro dell'iperbole
r_PH_n = R27 + [0 500 1000 1500 2000]; % km

% 3. Ciclo for: per ogni r_PH
for r_PH = r_PH_n

    % 3a. Caratterizzazione dell3orbita circolare di arrivo (posizione e velocità)
    r_oc = r_PH;
    v_oc = sqrt(mu_a/r_oc);
    v_oc_n = [v_oc_n, v_oc];

    % 3b. Determinazione dell%eccentricità dell iperbole:
    e_H = 1 - r_PH/a_H;
    if e_H < 1
        error('eccentricità minore di 1');
    end

    % 3c. Calcolo della velocità al pericentro dell'iperbole:
    v_fuga2 = sqrt(2*mu_a/r_PH);
    v_PH2 = sqrt(v_inf2^2+v_fuga2^2);
    v_PH2_n = [v_PH2_n, v_PH2];

    % 3d. Determinazione del parametro di impatto per ogni r_PH
    Delta_n = -a_H*sqrt(e_H^2 - 1);
    Delta_n_n = [Delta_n_n, Delta_n];

    % 3e. Calcolo del Dv di arrivo per ogni r_PH
    Delta_v = abs(v_oc - v_PH2);
    Delta_v_n = [Delta_v_n, Delta_v];
end

plot(Delta_n_n, (r_PH_n - ones(1,5).*R27), 'Marker','o',)


