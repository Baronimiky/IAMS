function [dv,omf,thf]=changeOrbitalPlane(ai,ei,Ii,OMi,omi,If,OMf,thi,P,mu)

% il tempo dt è il tempo necessario a raggiungere il punto di manovra a
% partire dal punto iniziale thi

% SE OMf-OMi<=pi:
%    se P=1: prende intersezione più vicina secondo senso di movimento a
%    partire dal nodo ascendente dell'orbita iniziale (se OMf-OMi=pi prende
%    intersezione più lontana dal nodo ascendente)
%    se P=2: prende altro punto di intersezione
% SE OMf-OMi>pi:
%    fa l'opposto

Ii=PositiveAngles(Ii);
omi=PositiveAngles(omi);
OMi=PositiveAngles(OMi);
If=PositiveAngles(If);
OMf=PositiveAngles(OMf);

if nargin==9
    mu=398600.433;
end

alpha=acos(cos(Ii)*cos(If)+sin(Ii)*sin(If)*cos(OMf-OMi));


if((OMf-OMi)>=0 && (If-Ii)>=0)
    ui=acos((cos(alpha)*cos(Ii)-cos(If))/(sin(alpha)*sin(Ii)));
    uf=acos((cos(Ii)-cos(alpha)*cos(If))/(sin(alpha)*sin(If)));
    if (OMf-OMi)>pi
        ui=-ui;
        uf=-uf;
    end
end
if((OMf-OMi)>=0 && (If-Ii)<0) 
    ui=-acos((cos(If)-cos(alpha)*cos(Ii))/(sin(alpha)*sin(Ii)));
    uf=-acos((cos(alpha)*cos(If)-cos(Ii))/(sin(alpha)*sin(If)));
    if (OMf-OMi)>pi
        uf=-uf;
        ui=-ui;
    end
end
if((OMf-OMi)<0 && (If-Ii)>=0)
    ui=-acos((cos(alpha)*cos(Ii)-cos(If))/(sin(alpha)*sin(Ii)));
    uf=-acos((cos(Ii)-cos(alpha)*cos(If))/(sin(alpha)*sin(If)));
    if (OMf-OMi)<-pi
       ui=-ui;
       uf=-uf;
    end  
end
if((OMf-OMi)<0 && (If-Ii)<0) 
    ui=acos((cos(If)-cos(alpha)*cos(Ii))/(sin(alpha)*sin(Ii)));
    uf=acos((cos(alpha)*cos(If)-cos(Ii))/(sin(alpha)*sin(If)));
    if (OMf-OMi)<-pi
        uf=-uf;
        ui=-ui;
    end
end    

if(ui>=omi)
    thf=ui-omi;
else
    thf=2*pi+(ui-omi);
end

if(uf>=thf)
    omf=uf-thf;
else
    omf=2*pi+(uf-thf);
end

if (P==2)
    thf=thf+pi;
end

dv=2*sqrt(mu/(ai*(1-ei^2)))*(1+ei*cos(thf))*sin(alpha/2);

