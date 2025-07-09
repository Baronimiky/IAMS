function plotOrbit(a,e,i,OM,om,th0,thf,dth,mu,color)
    if nargin < 10  % Se il colore non è fornito
        color = 'b';  % Default: blu
    end

    theta = th0:dth:thf;
    r_plot = zeros(3, length(theta));
    for k = 1:length(theta)
        r_plot(:,k) = par2car(a,e,i,OM,om,theta(k),mu);
    end

    plot3(r_plot(1,:), r_plot(2,:), r_plot(3,:), color, 'LineWidth', 1.5);
end
