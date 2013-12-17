function [dS]=sensitivity2_MAPK(x,rho,u,Jf,Sf)
m=2;
n=2;
%dJf_x
dJf_x(:,:,1)=[...
0,0;...
0,0];
dJf_x(:,:,2)=[...
0,0;...
0,0];
%dJf_rho
dJf_rho(:,:,1)=[...
-(2.0)*(1.0+u(1))^(-1),-(1.0+u(1))^(-1);...
(1.0+u(1))^(-1),0];
dJf_rho(:,:,2)=[...
-1,1;...
0,-1];
%dKf_rho
dKf_rho(:,:,1)=[...
0,0;...
0,0];
dKf_rho(:,:,2)=[...
0,0;...
0,0];
 A=zeros(n,n,m);
 for k=1:m
  for l=1:n
   A(:,l,k)=dJf_x(:,:,l)*Sf(:,k)+dJf_rho(:,l,k);
  end%for
 end%for
 B=zeros(n,m,m);
 for k=1:m
  for i=1:m
   B(:,i,k)=dJf_rho(:,:,i)*Sf(:,k);
  end%for
 end%for
 for k=1:m
  dS(:,:,k)=-Jf\(A(:,:,k)*Sf+B(:,:,k)+dKf_rho(:,:,k));
 end%for
end%function
