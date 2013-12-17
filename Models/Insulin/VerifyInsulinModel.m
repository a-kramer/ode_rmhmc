% calculate numerical Jacobian, Sensitivity and second order Sensitivity
% using finite differences
% and compare the result to the values in F.{Jf,Sf,dSf}
%

make_expInsulin_model


u=1;
% finite differences:
s=[+1,-1];
h=1e-8;
theta=[1.7288 0.1884 -1.7038 1.7470 -0.8761 -0.8593];
% Jacobian
Jf=zeros(F.n);
o=theta;

printf("this script calculates the missmatch between\n");
printf("F.{Jf,Sf,dSf} and the finite differences approximation\n");
printf("Parameters:\n"); theta
printf("state values: steady state.\n");

f=@(x) F.f(x,theta,u);
X=lsode(f,F.x0,[0,400]);
xs=X(2,:)'

for i=1:2
 for j=1:F.n
  x=xs;
  x(j)+=s(i)*h;
  Jf(:,j)+=s(i)*F.f(x,o,u)/(2*h);
 end%for
end%for
x=ones(F.n,1);
FJf=F.Jf(xs,theta,u)
printf("Jacobian difference norm: %g\trcond(Jf)=%g\n",norm(Jf-FJf),rcond(Jf));


% Sensitivity
Sf=zeros(F.n,F.m);
x=ones(F.n,1);


for i=1:2
 for j=1:F.m
  o=theta;
  o(j)+=s(i)*h;
  % find steady state
  f=@(x) F.f(x,o,u);
  X=lsode(f,F.x0,[0,400]);
  x=X(2,:)';
  Sf(:,j)+=s(i)*x/(2*h);
 end%for
end%for
o=theta;
printf("Sensitivity difference norm: %g\n",norm(Sf-F.Sf(x,o,u)));


% second order sensitivity
dSf=zeros(F.n,F.m,F.m);
x=ones(F.n,1);


for i=1:2
 for j=1:F.m
  o=theta;
  o(j)+=s(i)*h;
  % find steady state
  f=@(x) F.f(x,o,u);
  X=lsode(f,F.x0,[0,400]);
  x=X(2,:)';
  dSf(:,:,j)+=s(i)*F.Sf(x,o,u)/(2*h);
 end%for
end%for

FJf=F.Jf(xs,theta,u);
FSf=F.Sf(xs,theta,u);
FdSf=F.dSf(xs,theta,u,FJf,FSf);
N=zeros(1,F.m);

for j=1:F.m
 N(j)=norm(dSf(:,:,j)-FdSf(:,:,j));
end%for
printf("second order Sensitivity difference norm: ");
printf(" %g ",N);
printf("\n");