function deltat = TOF(a,e,th1,th2,mu)

%TOF Time Of Flight

%Input arg:
% a   [1x1]   semi-major axis         [km]
% e   [1x1]   eccentricity            [-]
% th1 [1x1]   Initial true anomaly    [rad]
% th2 [1x1]   Final true anomaly      [rad]
% mu  [1x1]   gravitational paramete  [km^3/s^2]
%----------------------------------------------
%Output arg:
% deltat   [1x1]   time of flight       [s]
%
%-----------------------------------------------

%anomalie eccentriche

T=2*pi*sqrt((a.^3)/mu);
%
% Nota:
%   Se theta2 risulta minore di theta1, si assume che theta2 appartenga al giro
%   successivo e si aggiunge 2*pi. La conversione in anomalia eccentrica
%   utilizza atan2 per evitare problemi numerici vicino a 0 e 2*pi.

    % Se theta2 è minore di theta1, si assume che theta2 appartenga al giro successivo.
        % orbita ellittica: conversione stabile
    E1=2*atan(sqrt((1-e)/(1+e))*tan(th1/2));
    E2=2*atan(sqrt((1-e)/(1+e))*tan(th2/2));

    if(E1<0)
        E1=E1+2*pi;
    end
    if(E2<0)
        E2=E2+2*pi;
    end

    M1 = E1 - e*sin(E1);
    M2 = E2 - e*sin(E2);

    deltat = real(sqrt(a^3 / mu) * (M2 - M1));

    if th2 < th1
        deltat = deltat + T;
    end
    if e>1 || e<0
        deltat = 1e100;
    end

end