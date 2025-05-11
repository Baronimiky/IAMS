clear all
close all
clc
data= load("DatiSC1-2025.txt");
Or_i = data(27,2:7);
Or_f = data(27,8:13);

% orbita iniziale
a_i = Or_i(1,1); % semiasse maggiore[km] 
e_i = Or_i(1,2); % eccentricità [-]
i_i = Or_i(1,3); % inclinazione del piano orbitale[rad] inc_i
OM_i = Or_i(1,4); % ascensione retta del nodo ascendente[rad] RAAN
om_i = Or_i(1,5); % anomalia del pericentro[rad] w
th_i = Or_i(1,6); % posizione del satellite nell'orbita [rad] theta_i
% orbita finale
rx_f = Or_f(1,1); % [km]
ry_f = Or_f(1,2); % [km]
rz_f = Or_f(1,3); % [km]
vx_f = Or_f(1,4); % [km/s]
vy_f = Or_f(1,5); % [km/s]
vz_f = Or_f(1,6); % [km/s]

rr_f = [rx_f; ry_f; rz_f]; % raggio punto di arrivo orbita
vv_f = [vx_f; vy_f;vz_f]; % velocità punto di arrivo orbita

mu_t = 398600; % costante gravitazionale terra [km^3/s^2]
mu_s = 132712440018; % costante gravitazionale sole [km^3/s^2]

[rr_i,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(rr_f,vv_f,mu_t);

dth = 0.01; % discriminante theta
th0 = 0; % theta iniziale
thf = 2*pi; %theta finale

Rt = 6378; % raggio terra
Terra3d(Rt); % plot terra 3d
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th0,thf,dth,th_i,mu_t); % plot orbita iniziale
plotOrbit(a_f,e_f,i_f,OM_f,om_f,th0,thf,dth,th_f,mu_t); % plot orbita finale


%[DeltaV, thi, thf] = changePericenterArg (a_i, e_i, om_i, om_f, mu_t)

%[DeltaV, omf, theta] = changeOrbitalPlane(a_i,e_i,i_i,OM_i,om_i,i_f,OM_f,mu_t)

%deltat = TOF(a_i,e_i,th0,thf,mu_t)

[DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pa', mu_t)