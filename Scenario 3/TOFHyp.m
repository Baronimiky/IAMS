function [theta_exit, delta_t] = TOFHyp(aH, eH, mu, theta0,r_SOI)
% deltaT_exit_SOI  Tempo per superare la SOI terrestre in orbita iperbolica
%
% Sintassi
%   [theta_exit, delta_t, r_SOI] = deltaT_exit_SOI(aH, eH, muE, theta0)
%
% Input
%   aH     : semiasse maggiore dell'iperbole (km)  [negativo]
%   eH     : eccentricità (>1)
%   muE    : costante gravitazionale del corpo centrale (km^3/s^2)
%            (di solito 3.986e5 per la Terra)
%   theta0 : anomalia vera di partenza (rad)  (default = 0, pericentro)
%
% Output
%   theta_exit : anomalia vera in cui r = r_SOI (rad)
%   delta_t    : tempo trascorso da theta0 a theta_exit (s)
%   r_SOI      : raggio della SOI terrestre calcolato (km)
%
% Note
%   La SOI è calcolata con la formula classica di Laplace:
%       r_SOI = a_Earth * (m_Earth/m_Sun)^(2/5)
%   dove a_Earth = 1 AU = 1.495978707e8 km e m_Earth/m_Sun ≈ 3.003489e-6.

    if nargin < 4, theta0 = 0; end
    if aH >= 0, error('Per un''iperbole devi fornire aH < 0'); end
    if eH <= 1, error('Eccentricità iperbolica: eH deve essere > 1'); end

    % ------------------ 1. r(θ) = r_SOI  -------------------------------
    A   = -aH;                            % |a|
    cos_theta_exit = (A*(eH^2 - 1)/r_SOI - 1) / eH;

    if abs(cos_theta_exit) > 1
        error('La SOI non viene mai superata con questi parametri.');
    end

    theta_candidate = acos(cos_theta_exit);   % rad ∈ [0, π]

    % ramo “in uscita”: dobbiamo assicurarci theta_exit > theta0
    theta_exit = theta_candidate;
    if theta_exit <= theta0
        theta_exit = 2*pi - theta_candidate;
    end

    % ------------------ 3. Δt tramite anomalia iperbolica --------------
    theta2H = @(th) asinh( sin(th)*sqrt(eH^2 - 1) ./ (1 + eH*cos(th)) );
    H0 = theta2H(theta0);
    Ht = theta2H(theta_exit);

    Mh  = @(H) eH*sinh(H) - H;            % M_h(H)
    delta_t = sqrt(A^3 / mu) * (Mh(Ht) - Mh(H0));

    % valore assoluto per sicurezza
    delta_t = abs(delta_t);
end