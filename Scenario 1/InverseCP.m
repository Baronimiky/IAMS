function [omi]=InverseCP(Ii,OMi,omf,If,OMf)

%% FUNZIONE NON MIA, DA MODIFICARE


% calcola omi a partire da omf noti i parametri dei piani orbitali delle 2
% orbite

Ii=PositiveAngles(Ii);
omf=PositiveAngles(omf);
OMi=PositiveAngles(OMi);
If=PositiveAngles(If);
OMf=PositiveAngles(OMf);


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

if (uf>=omf)
    thf=uf-omf;
else
    thf=2*pi+(uf-omf);
end

if (ui>=thf)
    omi=ui-thf;
else
    omi=2*pi+(ui-thf);
end

omi=PositiveAngles(omi);