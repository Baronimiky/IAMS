function a=PositiveAngles(alpha)

% Normalizza l'angolo in radianti in [0,2*pi]

if (alpha>=0)
    a=mod(alpha,2*pi);
else
    a=mod(alpha,-2*pi)+2*pi;
end

if(a==0 && alpha~=0)
        a=2*pi;
end
