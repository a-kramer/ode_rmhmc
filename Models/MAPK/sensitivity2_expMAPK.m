function [dS]=sensitivity2_expMAPK(x,theta,u,Jf,Sf)
m=2;
n=2;
%dJf_x
dJf_x(:,:,1)=[...
0,0;...
0,0];
dJf_x(:,:,2)=[...
0,0;...
0,0];
%dJf_theta
dJf_theta(:,:,1)=[...
-(2.0)*(1.0+u(1))^(-1)*exp(theta(1)),-(1.0+u(1))^(-1)*exp(theta(1));...
(1.0+u(1))^(-1)*exp(theta(1)),0];
dJf_theta(:,:,2)=[...
-exp(theta(2)),exp(theta(2));...
0,-exp(theta(2))];
%dKf_theta
dKf_theta(:,:,1)=[...
-x(2)*(1.0+u(1))^(-1)*exp(theta(1))-(2.0)*x(1)*(1.0+u(1))^(-1)*exp(theta(1))+u(1)*(1.0+u(1))^(-1)*exp(theta(1)),0;...
x(1)*(1.0+u(1))^(-1)*exp(theta(1)),0];
dKf_theta(:,:,2)=[...
0,x(2)*exp(theta(2))-x(1)*exp(theta(2));...
0,-x(2)*exp(theta(2))];
 A=zeros(n,n,m);
 for k=1:m
  for l=1:n
   A(:,l,k)=dJf_x(:,:,l)*Sf(:,k)+dJf_theta(:,l,k);
  end%for
 end%for
 B=zeros(n,m,m);
 for k=1:m
  for i=1:m
   B(:,i,k)=dJf_theta(:,:,i)*Sf(:,k);
  end%for
 end%for
 for k=1:m
  dS(:,:,k)=-Jf\(A(:,:,k)*Sf+B(:,:,k)+dKf_theta(:,:,k));
 end%for
end%function
