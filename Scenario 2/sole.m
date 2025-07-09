function sole(Rm)

[X, Y, Z] = sphere(50);          % sfera

if nargin==0
      Rm = 6.955 *1e+5;
end  

        % Plotta il sole
        HSUN = surf(Rm*X,Rm*Y,Rm*Z);
        % Load Sun Image
        topoSUN = imread('sun.jpg');
        % immagine sul sole
        set(HSUN,'facecolor','texture',...
                'cdata',im2double(topoSUN),...
                'edgecolor','none');

        set(gca,'color','w');
        axis equal
        rotate3d;
end

