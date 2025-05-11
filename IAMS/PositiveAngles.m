function a=PositiveAngles(alpha)

%fa in modo che l'angolo in radianti sia sempre espresso come positivo e
%compreso tra 0 e 2pi

if (alpha>=0)
    a=mod(alpha,2*pi);
else
    a=mod(alpha,-2*pi)+2*pi;
end

if(a==0 && alpha~=0)
        a=2*pi;
end
