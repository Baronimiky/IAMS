function asteroid(Ra, pos)
% asteroid(Ra, pos)
% Plotta un asteroide con raggio Ra centrato nel punto pos = [x0, y0, z0]
% Ra: raggio dell'asteroide (default 6.955e6 m)
% pos: coordinate del centro (default [0, 0, 0])

    if nargin < 1
        Ra = 6.955e6;
    end
    if nargin < 2
        pos = [0, 0, 0];  % centro sull'origine se non specificato
    end

    x0 = pos(1);
    y0 = pos(2);
    z0 = pos(3);

    % Sfera di riferimento
    [X, Y, Z] = sphere(50);

    % Crea la superficie dell'asteroide spostata in (x0, y0, z0)
    HAST = surf(Ra*X + x0, Ra*Y + y0, Ra*Z + z0);

    % Carica la texture
    topoAST = imread('asteroid.png');

    % Applica la texture
    set(HAST, 'FaceColor', 'texture', ...
              'CData', im2double(topoAST), ...
              'EdgeColor', 'none');

    set(gca, 'Color', 'w');
    axis equal
    rotate3d;
end
