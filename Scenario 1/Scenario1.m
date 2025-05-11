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

%[rr_i,vv_i] = par2car(a_i,e_i,i_i,OM_i,om_i,th_i,mu_t);

[a_f,e_f,i_f,OM_f,om_f,th_f] = car2par(rr_f,vv_f,mu_t);

dth = 0.01; % discriminante theta
th0 = 0; % theta iniziale
thf = 2*pi; %theta finale

Rt = 6378; % raggio terra
Terra3d(Rt); % plot terra 3d
plotOrbit(a_i,e_i,i_i,OM_i,om_i,th0,thf,dth,th_i,mu_t);
% plot orbita iniziale
plotOrbit(a_f,e_f,i_f,OM_f,om_f,th0,thf,dth,th_f,mu_t); % plot orbita finale

% manovra standard

% Conversione delle coordinate cartesiane in parametri orbitali
[a_f, e_f, i_f, OM_f, om_f, th_f] = car2par(rr_f,vv_f,mu_t);
%
% Cambio del piano orbitale
[DeltaV_CP, om_2, theta_CP] = changeOrbitalPlane(a_i, e_i, i_i, OM_i, om_i, i_f, OM_f,2,mu_t);
DeltaT1 = TOF(a_i, e_i, th_i, theta_CP,mu_t); % tempo manovra
% plot punto di partenza manovra
[rr_1] = par2car(a_i,e_i,i_i,OM_i,om_i,theta_CP,mu_t); % punto partenza manovra
plot3(rr_1(1),rr_1(2),rr_1(3),'k+','LineWidth',2);

plotOrbit(a_i,e_i,i_f,OM_f,om_2,0,2*pi,dth,theta_CP,mu_t) % prima manovra

% Cambio dell'argomento del pericentro
[DeltaV_CA, thi_CPA, thf_CPA] = changePericenterArg(a_i, e_i, om_i, om_f, mu_t);

Deltat2_1 = TOF(a_i,e_i,theta_CP,thi_CPA(1),mu_t); % tempo manovra
Deltat2_2 = TOF(a_i,e_i,theta_CP,thi_CPA(2),mu_t); % tempo manovra

if Deltat2_1 > Deltat2_2
    DeltaT2 = Deltat2_2;
    DeltaT3 = TOF(a_i,e_i,thf_CPA(2),0,mu_t); % Tempo di volo per ulteriore trasferimento
    thf_CPAs = thf_CPA(2);
else
    DeltaT2 = Deltat2_1;
    DeltaT3 = TOF(a_i,e_i,thf_CPA(1),0,mu_t); % Tempo di volo per ulteriore trasferimento
     thf_CPAs = thf_CPA(1);
end
% Selezionare theta opportuni
plotOrbit(a_i,e_i,i_f,OM_f,om_f,0,2*pi,dth,thf_CPAs,mu_t) % seconda manovra
%%
% Trasferimento bitangente
[DeltaV1_BT, DeltaV2_BT, DeltaT4] = bitangentTransfer(a_i, e_i, a_f, e_f, 'pa', mu_t);

% Ultimo tempo di volo
DeltaT5 = TOF(a_f, e_f, pi, th_f,mu_t);

plotOrbit(a_f,e_f,i_f,OM_f,om_f,0,2*pi,dth,pi,mu_t) % terza manovra
legend('Orbita iniziale','orbita finale','Cambio piano','Cambio angolo','trasferimento bitangente');

Delta_VTot = abs(DeltaV_CP) + abs(DeltaV_CA) +abs(DeltaV1_BT) +abs(DeltaV2_BT) % deltav totale [m/s]
Delta_TTot = DeltaT1 + DeltaT2 + DeltaT3 + DeltaT4 + DeltaT5 % delta Tempo totale [s]