function [XYZ, craterTbl, handles] = craters(NameValueArgs)
%CRATERS creates a cratered sphere.
%   CRATERS() produces and plots a cratered sphere on the current
%   axes centered at (0,0,0) with radius=1 using default
%   parameter values. Crater locations are randomly sampled
%   uniformly across the spherical surface. Crater radii are
%   are defined by a gramma distribution and their depths are
%   proportional to radii. Crater profiles are defined by a
%   generalized logistic function.
%
%   XYZ = CRATERS(__) returns an N*N*3 array of 3D coordinates for
%   N^2 points of the cratered sphere.
%
%   [XYZ, CRATERTBL] = CRATERS(__) returns a table containing
%   the (X,Y,Z) coordinates of each crater center, crater radius,
%   and crater depth. These coordinates lie on the surface of the
%   sphere prior to crater impacts.
%
%   [XYZ, CRATERTBL, HANDLES] = CRATERS(__) returns a structure
%   containing graphics handles to objects created or accessed
%   within this function which may include,
%       .fig: figure handle
%       .ax: axes handle
%       .surf: surface handle to cratered sphere
%       .light1: light object 1
%       .light2: light object 2
%
%   NAME-VALUE ARGUMENTS
%   sphPos - sphere position
%       [0,0,0] (default) | 1x3 vector
%   sphRad - sphere radius
%       1 (default) | scalar >0
%   sphScale - sphere size scaling options
%       'pre' (default) | 'post
%        - pre: sphere is scaled to sphRad prior to impacts
%        - post: sphere is scaled to sphRad after impacts
%        * Crater sizes are set according to pre-impact scaling.
%   nCraters - maximum number of craters
%       1000 (default) | scalar >=0
%       Fewer craters may appear due to randomized sizes.
%   distribution - distribution of craters
%       'random' (default) | 'fiblat' (fibonacci lattice)
%   craterRadLim - crater radii limits
%       [0.02, 0.5] (default) | 1x2 vector  of real numbers >=0
%   depthScale - scales crater depths (D) by D=depthScale*radius
%       0.2 (default) | scalar >=0
%   n - number of points that defines the surface is n^2
%       500 (default) | scalar >1
%       n should be sufficiently large enough to see craters
%   destination - specifies plot destination
%       axis handle (default: gca) | figure handle | 'outputs' to
%       prevent plotting.
%   demo - generates built-in demos
%       'asteroid' | 'fiblat' | 'moon' | 'none' (default)
%       - asteroid: heavily impacted sphere
%       - fiblat: unformly distributed impacts of equal radii
%       - moon: add cosmetics to the default craters() object
%       - none: uses inputs and default values 
%
% Source: <a href = "https://www.mathworks.com/matlabcentral/fileexchange/101649-craters">Matlab File exchange</a>
% Contact author: <a href = "https://www.mathworks.com/matlabcentral/profile/authors/3753776">Adam Danz</a>
% Copyright (c) 2021  All rights reserved
%
%   See also CRATERSEXAMPLES, SURF.

% Compatible with Matlab R2020b and later
%   Limited by mustbetextscalar (R20b)

% Version history
% 0.0.0  01-Nov-2021 279 character version submitted to Mini Hack contest
%       https://www.mathworks.com/matlabcentral/communitycontests/contests/4/entries/9623
% 1.0.0  07-Nov-2021 Upload to Matlab central file exchange (Matlab R2021b)
% 1.0.1  07-Nov-2021 Updated URLs to FEX page

%% Input validation
arguments
    NameValueArgs.sphPos (1,3) double       {mustBeFinite, mustBeNonNan, mustBeReal} = [0,0,0]
    NameValueArgs.sphRad (1,1) double       {mustBeFinite, mustBeNonNan, mustBeReal, mustBePositive} = 1
    NameValueArgs.sphScale string           {mustBeTextScalar} = 'pre';
    NameValueArgs.nCraters (1,1) double     {mustBeFinite, mustBeNonNan, mustBeReal, mustBeNonnegative} = 1000
    NameValueArgs.distribution string       {mustBeTextScalar} = 'random';
    NameValueArgs.craterRadLim (1,2) double {mustBeFinite, mustBeNonNan, mustBeReal, mustBeNonnegative} = [.02, .5];
    NameValueArgs.depthScale (1,1) double   {mustBeFinite, mustBeNonNan, mustBeReal, mustBeNonnegative} = 0.2
    NameValueArgs.n (1,1) double            {mustBeFinite, mustBeNonNan, mustBeReal, mustBeGreaterThan(NameValueArgs.n,1)} = 500
    NameValueArgs.destination % validated below
    NameValueArgs.demo string               {mustBeTextScalar} = 'none';
end

% Override NameValueArgs if demo is specified
demolist = {'none','demolist','asteroid','fiblat','moon'}; 
NameValueArgs.demo = validatestring(NameValueArgs.demo, demolist);
switch NameValueArgs.demo
    case 'none'
        % Do nothing
    case 'demolist'
        % undocumented; returns list of demos
        XYZ = setdiff(demolist, {'none','demolist'});
        nargoutchk(0,1)
        return
    case 'asteroid'
        NameValueArgs.sphPos = [0,0,0];
        NameValueArgs.sphRad = 1;
        NameValueArgs.nCraters = 100;
        NameValueArgs.distribution = 'random';
        NameValueArgs.craterRadLim = [1, 1.5];
        NameValueArgs.depthScale = 0.35;
        NameValueArgs.n = 500;
    case 'fiblat'
        NameValueArgs.sphPos = [0,0,0];
        NameValueArgs.sphRad = 1;
        NameValueArgs.nCraters = 340;
        NameValueArgs.distribution = 'fiblat';
        NameValueArgs.craterRadLim = [.15, .15];
        NameValueArgs.depthScale = .2;
        NameValueArgs.n = 500;
    case 'moon'
        NameValueArgs.sphPos = [0,0,0];
        NameValueArgs.sphRad = 1;
        NameValueArgs.nCraters = 1000;
        NameValueArgs.distribution = 'random';
        NameValueArgs.craterRadLim = [.02, .5];
        NameValueArgs.depthScale = 0.2;
        NameValueArgs.n = 500;
end

% Validate destination input
% NameValueArgs.destination input can be missing, an axis handle, a figure handle,
% or a string/char "outputs" (accepts partial match).  Default NameValueArgs.destination
% when missing is gca. Fig handles is replaced with new axis handle.
if ~isfield(NameValueArgs, 'destination') || isempty(NameValueArgs.destination)
    NameValueArgs.destination = gca();
elseif ischar(NameValueArgs.destination) || isa(NameValueArgs.destination,'string')
    NameValueArgs.destination = validatestring(NameValueArgs.destination,{'outputs'});
elseif isgraphics(NameValueArgs.destination,'figure')
    NameValueArgs.destination = axes(NameValueArgs.destination);
elseif isgraphics(NameValueArgs.destination,'axes')
    % do nothing
else
    error('"Destination" input must be a figure or axes handle or the string "outputs".')
end

% Other arg checks
NameValueArgs.distribution = validatestring(NameValueArgs.distribution, {'random','fiblat'});
NameValueArgs.sphScale = validatestring(NameValueArgs.sphScale, {'pre','post'});

%% Create cratered spherical surface
% Generate unit sphere coordinates
[Xs,Ys,Zs]=sphere(NameValueArgs.n-1);

% Define crater locations
switch NameValueArgs.distribution
    % Produce craterLoc; an M*3 matrix of [x,y,z] crater locations
    case 'random'
        % Random sample from a uniform distribution over sphere surface.
        % Inspired by Roger Stafford's randsphere() from the File Exchange.
        rand3D = randn(NameValueArgs.nCraters, 3);
        craterLoc = rand3D.*repmat(1./sqrt(sum(rand3D.^2,2)),1,3);
    case 'fiblat'
        % fiblat lattice; uniform distribution of craters [footnote 4]
        goldenrat = (1 + 5^0.5)/2;
        idx = (1:NameValueArgs.nCraters)';  % col vec
        theta = 2*pi*idx/goldenrat;
        phi = acos(1 - 2*(idx+0.5)/(NameValueArgs.nCraters+1));
        craterLoc = [cos(theta) .* sin(phi), ...
            sin(theta) .* sin(phi), ...
            cos(phi)];
end
% View crater loc distribution:
%   figure(); plot3(craterLoc(:,1), craterLoc(:,2), craterLoc(:,3), 'k.'); axis equal; grid on

% Define crater radii (gamma distribution, scaled).
craterRadii = rescale(randg(1,1,NameValueArgs.nCraters),...
    NameValueArgs.craterRadLim(1)/NameValueArgs.sphRad, ...
    NameValueArgs.craterRadLim(2)/NameValueArgs.sphRad);
% View distribution of crater radii
%   figure(); histogram(craterRadii,20); xlabel('Crater radii'); ylabel('Count'); grid on;

% Compute linear distance from points along the sphere to each crater center.
% This is much faster than using spherical distance but still takes a few seconds [footnote 1].
distFcn=@(M,i)(M(:)-craterLoc(:,i)').^2;
distToCraters=sqrt(distFcn(Xs,1)+distFcn(Ys,2)+distFcn(Zs,3));  % N*M for N sphere points and M craters

% Crater depths are a function of crater radii so that larger craters make deeper impact.
% Crater profiles are defined by a logistic function [footnote 2,3].
craterDepths = NameValueArgs.depthScale * craterRadii;
D = distToCraters./craterRadii;
normProfile = 1./((1+.4.*exp(-10*D+4)).^2);
sphRadii = normProfile.*craterDepths+1-craterDepths; % [footnote 3]

% Modify sphere coordinates to add craters
if isempty(sphRadii)
    % When sphRadii is empty, there are no craters.
    Xe=Xs; Ye=Ys; Ze=Zs;
else
    % The radius of each spherical point is determined by the deepest
    % crater that affected that point, or min(sphRadii,[],2).
    [T,P]=cart2sph(Xs,Ys,Zs);
    R=reshape(min(sphRadii,[],2),NameValueArgs.n,NameValueArgs.n);
    if strcmpi(NameValueArgs.sphScale, 'post')
        % Scale R so max(R)==1
        R = R*(1/max(R(:)));
    end
    [Xe,Ye,Ze]=sph2cart(T,P,R);
end

% Scale to size and location
Xe = NameValueArgs.sphRad * Xe + NameValueArgs.sphPos(1);
Ye = NameValueArgs.sphRad * Ye + NameValueArgs.sphPos(2);
Ze = NameValueArgs.sphRad * Ze + NameValueArgs.sphPos(2);


%% Plot results
if ~ischar(NameValueArgs.destination)
    handles.fig = ancestor(NameValueArgs.destination,'figure');
    handles.ax = NameValueArgs.destination;
    handles.surf = surf(handles.ax,Xe,Ye,Ze,'EdgeColor','none','FaceColor',[.8 .8 .8]); % also try facecolor='interp'

    if ~strcmpi(NameValueArgs.demo,'none')
        % Set demo properties
        % set(handles.ax, 'Color', 'k'); % these are helpful during development.
        % handles.ax.XAxis.Color = 'w';
        % handles.ax.YAxis.Color = 'w';
        % handles.ax.ZAxis.Color = 'w';
        % ax.GridColor = 'w';
        % ax.GridAlpha = .85;
        % box(ax,'on')
        set(handles.fig, 'Color', 'k')
        axis(handles.ax, 'equal')
        axis(handles.ax, 'vis3d')
        axis(handles.ax, 'off')
        handles.light1 = light(handles.ax,'Position',[1 0 1],'Color',[1 1 1]);
        handles.light2 = light(handles.ax,'Position',[-1 0 -1],'Color',[.2 .2 .2]);

        switch NameValueArgs.demo
            % Define materialType and do demo-specific settings
            case 'asteroid'
                materialType = 'dull';
            case 'fiblat'
                materialType = 'shiny';
                view([150,15])
            case 'moon'
                handles.surf.FaceColor = 'i';
                colormap(handles.ax, gray(265))
                caxis(handles.ax, [-.4,1.1])
                materialType = 'dull';
                camup([0, .6 .8])
                camva(7)
        end
        material(handles.surf, materialType)
    end
else
    handles = struct();
end

%% Prepare outputs
if nargout>0
    XYZ = cat(3,Xe,Ye,Ze);
end
if nargout>1
    craterTbl = array2table([craterLoc,craterRadii(:),craterDepths(:)], ...
        'VariableNames',{'X','Y','Z','Radius','Depth'});
    craterTbl.Properties.Description = sprintf('Crater properties produced in %s.m.',mfilename);
end

%% Footnotes
% [1] Spherical distance
%     Compute spherical distance from sphere points to crater centers.
%     This takes way too much time for large n and NameValueArgs.nCraters.
%     Haversine formula: https://en.wikipedia.org/wiki/Great-circle_distance
%         [T,P,R]=cart2sph(X,Y,Z);          % duplicate of line below
%         [Tc,Pc,Rc]=cart2sph(craterLoc(:,1),craterLoc(:,2),craterLoc(:,3));
%         dP = P(:)-Pc(:)';
%         dT = T(:)-Tc(:)'; % isn't there a boundary issue at [0,360] ?
%         distToCraters = 2*asin(sqrt(sin(dP/2).^2 + cos(P(:)).*cos(Pc(:)').*sin(dT/2).^2));
%
% [2] Crater profiles are defined by a logistics function
%   https://en.wikipedia.org/wiki/Generalised_logistic_function
%                   A
%      P =  ------------------
%           (B + C*e^(-Dx+F))^G
%   A: signed height of curve
%   B: typically equals 1
%   C: Shifts curve rightward nonlinearly as C increases
%       When Q=1, half point of curve is at x=.5,y=.5
%   D: Growth rate; steepness
%   F: horizontal shift in direction of sign
%   G: affects growth rate and shape
%
%   Demo of each parameter:
%     x = -1:.01:2;
%     f = @(p)p(1)./((p(2)+p(3).*exp(p(4)*x+p(5))).^p(6));
%     %   1./((1+.4.*exp(-10*x+4)).^2);
%
%     bounds = [-1 2;
%         1 1;
%         0 2;
%         -10 -5;
%         -6 6;
%         .1 5]; % [min, max for each input param
%
%     defaults = [1 1 .4 -10 4 2];
%
%     figure()
%     tlo = tiledlayout(2,3);
%     for i = 1:6
%         nexttile
%         title(sprintf('param %d', i))
%         hold on
%         vals = unique(linspace(bounds(i,1), bounds(i,2), 5));
%         for j = 1:numel(vals)
%             inputs = defaults;
%             inputs(i) = vals(j);
%             y = f(inputs);
%             plot(x,y,'DisplayName',sprintf('%.2f',vals(j)))
%         end
%         legend
%         grid on
%     end
%     title(tlo, func2str(f),'interpreter','none')
%
% [3] Crater profiles demo
%   In this demo, set the crater radius and depthScale to see the crater's profile
%       distanceToCrater = 0:0.01:3;
%       craterRadius = 1.5;
%       depthScale = 0.2;
%       craterDepth = depthScale * craterRadius;
%       D = distanceToCrater/craterRadius;
%       normProfile = 1./((1+.4.*exp(-10*D+4)).^2); % make sure these values match main fcn
%       scaledProfile = normProfile * craterDepth + 1 - craterDepth;
%       figure();
%       s1 = subplot(2,1,1);
%       plot(distanceToCrater, normProfile); grid on; axis padded
%       xline(0); yline(0); xline(craterRadius, 'k--', 'Crater Radius')
%       xlabel('Distance to crater center'); ylabel('Normalized depth')
%       title('Logistic function'); axis equal
%       s2 = subplot(2,1,2);
%       plot(distanceToCrater, scaledProfile); grid on; axis padded; linkaxes([s1,s2],'y')
%       xline(0); yline(0); xline(craterRadius, 'k--', 'Crater Radius'); axis equal
%       xlabel('Distance to crater center'); ylabel('Distance to sphere center')
%       title('Crater profile function'); subtitle(sprintf('depthScale = %g', depthScale))
%
% [4] See "How to evenly distribute points on a sphere more effectively than the canonical Fibonacci Lattice"
%   http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/

%% Future ideas
% 1. Allow user to specify crater radii and depths
% 2. Add Goldberg polyhedron to the distribution options
% 3. Demo: "Luna" using real data from moon
% 4. Add additional crater profile options.
%    - profile with upper cusp
%    - consider using gamrnd instead of randg


