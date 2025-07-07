function sole(Rm)

[X, Y, Z] = sphere(50);          % Reference Sphere

if nargin==0
      Rm = 6.955 *1e+6;
end  

        % Place SUN
        HSUN = surf(Rm*X,Rm*Y,Rm*Z);
        % http://stereo.gsfc.nasa.gov/img/stereoimages/preview/euvisdoCarringtonMap.jpg
        % Load Sun Image
        topoSUN = imread('euvisdoCarringtonMap.jpg');
        % Set it on SUN
        set(HSUN,'facecolor','texture',...
                'cdata',im2double(topoSUN),...
                'edgecolor','none');

        set(gca,'color','w');
        axis equal
        rotate3d;
end

