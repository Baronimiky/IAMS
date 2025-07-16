function [omi]=InverseCP(Ii,OMi,omf,If,OMf)
% Inverse Change of Plane maneuver
%
% [omi] = InverseCP(Ii, OMi, omf, If, OMf)
%
% -------------------------------------------------------------------------
% Scopo:
% Calcola l’anomalia del pericentro iniziale `omi` che permette una manovra
% di cambio di piano orbitale tra due orbite con parametri noti, in modo da
% ottenere una specifica anomalia del pericentro finale `omf`, incl. e RAAN.
%
% -------------------------------------------------------------------------
% Argomenti di input:
% Ii       [1x1]    inclinazione orbitale iniziale          [rad]
% OMi      [1x1]    RAAN iniziale (ascensione retta del nodo ascendente) [rad]
% omf      [1x1]    anomalia del pericentro finale desiderata [rad]
% If       [1x1]    inclinazione orbitale finale            [rad]
% OMf      [1x1]    RAAN finale                             [rad]
%
% -------------------------------------------------------------------------
% Argomenti di output:
% omi      [1x1]    anomalia del pericentro iniziale calcolata [rad]
%
% -------------------------------------------------------------------------
% Note:
% - La funzione sfrutta la geometria tra i due piani orbitali per determinare
%   la posizione relativa dei pericentri.
% - Si assume che il cambiamento di piano avvenga in un singolo punto, e si
%   considera solo la geometria dei due piani orbitali.
% - Utilizza la formula dell’angolo tra piani orbitali (alpha).
%

% Porta tutti gli angoli in [0, 2*pi]
Ii  = PositiveAngles(Ii);
OMi = PositiveAngles(OMi);
omf = PositiveAngles(omf);
If  = PositiveAngles(If);
OMf = PositiveAngles(OMf);

% Differenze angolari normalizzate portate tra [-pi,pi]
dOM = mod(OMf - OMi + pi, 2*pi) - pi;

% Angolo tra i due piani orbitali
alpha = acos(cos(Ii)*cos(If) + sin(Ii)*sin(If)*cos(dOM));

sin_alpha = sin(alpha);

x_ui = (cos(alpha)*cos(Ii) - cos(If)) / (sin_alpha*sin(Ii));
x_uf = (cos(Ii) - cos(alpha)*cos(If)) / (sin_alpha*sin(If));

ui = atan2(sqrt(1 - x_ui^2), x_ui);
uf = atan2(sqrt(1 - x_uf^2), x_uf);

% Correggo i segni
if dOM < 0
    ui = -ui;
    uf = -uf;
end
if (If - Ii) < 0
    ui = -ui;
    uf = -uf;
end

ui = PositiveAngles(ui);
uf = PositiveAngles(uf);

% Calcolo dell'anomalia vera al punto di manovra
thf = mod(uf - omf, 2*pi);

% Calcolo dell'anomalia del pericentro iniziale
omi = mod(ui - thf, 2*pi);

omi = PositiveAngles(omi);

end