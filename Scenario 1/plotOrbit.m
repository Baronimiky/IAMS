function plotOrbit(a,e,i,OM,om,th0,thf,dth,th,mu)

 x = th0:dth:thf;
 r_plot = [];
 for j = th0:dth:thf
     [rr] = par2car(a,e,i,OM,om,j,mu); % matrice di tutti i punti della nostra orbita
     r_plot = [r_plot, rr];
 end
 hold on
% [rr_p] = par2car(a,e,i,OM,om,th,mu); % punto della nostra orbita
 plot3(r_plot(1,:),r_plot(2,:),r_plot(3,:));
 %plot3(rr_p(1),rr_p(2),rr_p(3),'k+','LineWidth',2);
end

