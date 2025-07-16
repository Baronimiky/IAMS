function Terra3d(Rt, pos)
% Terra3d(Rt, pos)
% Plotta la Terra in 3D con raggio Rt centrata nel punto pos = [x0, y0, z0]
% Rt: raggio della Terra (default 6378 km)
% pos: coordinate del centro della Terra (default [0, 0, 0])

    if nargin < 1
        Rt = 6378;
    end
    if nargin < 2
        pos = [0, 0, 0];  % centro sull'origine se non specificato
    end

    x0 = pos(1);
    y0 = pos(2);
    z0 = pos(3);

    % Carica i dati topografici
    load('topo.mat','topo','topomap1');

    % Mappa 2D di riferimento (non traslata nello spazio)
    contour(0:359, -89:90, topo, [0 0], 'b')
    axis equal
    hold on
    image([0 360], [-90 90], topo, 'CDataMapping', 'scaled');
    colormap(topomap1)

    % Sfera 3D (Terra)
    [x, y, z] = sphere(100);
    props.AmbientStrength = 0.1;
    props.DiffuseStrength = 1;
    props.SpecularColorReflectance = .5;
    props.SpecularExponent = 20;
    props.SpecularStrength = 1;
    props.FaceColor = 'texture';
    props.EdgeColor = 'none';
    props.FaceLighting = 'phong';
    props.Cdata = topo;

    % Sposta la Terra in (x0, y0, z0)
    surface(x*Rt + x0, y*Rt + y0, z*Rt + z0, props);

end