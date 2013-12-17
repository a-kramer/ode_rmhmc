MAPKplots

zq=quantile(z,[0.05,0.5,0.95]);

figure(1);
h=errorbar([U{:}],[Data{1,:}],[Data{2,:}]); hold on;
set(h,"marker",".");
set(h,"markersize",6);
set(h,"linestyle","none");
set(h,"linewidth",2);
plot(u,zq','-;;','linewidth',2,"color",[0.4,0,0]);
hold off;
