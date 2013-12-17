% mifa Insulin Dose Response Model

pkg load symbolic
symbols

n=6;
C=zeros(1,n);
C(1,5)= 13740.432 % taken from the Mifa.txt model file
rc=rows(C);
m=14;
l=1;

model_name="Mifa"
filename=sprintf("make_%s_model.m",model_name);

x=cell(n,1);
rho=cell(m,1);
u=cell(l,1);

x0="[10;0;0;0;0;0]";
# 1 IR(0) = 10                                                                                                                                                                                                                    
# 2 IRins(0) = 0                                                                                                                                                                                                                  
# 3 IRp(0) = 0                                                                                                                                                                                                                    
# 4 IRiP(0) = 0                                                                                                                                                                                                                   
# . IRi(0) = 0                                                                                                                                                                                                                    
# . IRS(0) = 10                                                                                                                                                                                                                   
# 5 IRSiP(0) = 0                                                                                                                                                                                                                  
# . X(0) = 10                                                                                                                                                                                                                     
# 6 Xp(0) = 0

% set up the state variables, sampling parameters and input parameters
for i=1:n
 x{i}=sym(sprintf('x(%i)',i));
endfor
for i=1:m
 theta{i}=sym(sprintf('theta(%i)',i));
 rho{i}=Exp(theta{i});
endfor
for i=1:l
 u{i}=sym(sprintf('u(%i)',i));
endfor

# fluxes:
# V1a = k1a · ins · [IR] + k1abasic · [IR]
# V1b = k1b · [IRins]
# V1c = k1c · [IRins]
# V1d = k1d · [IRp]
#                       k1f · [Xp]
# V1e = [IRip] · (k1e + ----------)
#                        1 + [Xp]
# 
# V1g = k1g · [IRp]
# V1r = k1r · [IRi]
# V2 = k21 · ([IRp] + k22 · [IRip]) · [IRS]
# Vm2 = km2 · [IRSip]
# V3 = k3 · [X] · [IRSip]
# Vm3 = km3 · [Xp]

# conservation relations: x{5}=10-(x{1}+x{2}+x{3}+x{4})
IRi=10-(x{1}+x{2}+x{3}+x{4});
IRS=10-x{5};
X=10-x{6};

flux={
rho{1}*u{1}*x{1} + rho{2}*x{1};
rho{3}*x{2};
rho{13}*x{2};
rho{14}*x{3};
x{4}*(rho{4}+rho{5}*x{6}/(1+x{6}));
rho{6}*x{3};
rho{7}*IRi;
rho{8}*(x{3}+rho{9}*x{4})*IRS;
rho{10}*x{5};
rho{11}*X*x{5};
rho{12}*x{6}
};

nr=length(flux);

%set up the right-hand-side of the ode (rhs)
#spc.↓dot↓
# -1- IR = [V1a] + [V1b] + [V1r] + [V1g]
# -2- IRins = [V1a] - [V1b] - [V1c]
# -3- IRp = [V1c] - [V1d] - [V1g]
# -4- IRip = [V1d] - [V1e]
# --- IRi = [V1e] - [V1r]           % 10-IR-IRins-IRp-IRip
# --- IRS = [V2] + [Vm2]            % 10-IRSip
# -5- IRSip = [V2] - [Vm2]
# --- X = [V3] + [Vm3]              % 10-X
# -6- Xp = [V3] - [Vm3]
#[1-5] V1a V1b V1c V1d V1e
#[6-11] V1g V1r V2 Vm2 V3 Vm3]
f={
% conservation relation 1:
 flux{2} + flux{7} + flux{6} - flux{1}; %1
 flux{1} - flux{2} - flux{3};           %2
 flux{3} - flux{4} - flux{6};           %3
 flux{4} - flux{5};                     %4
 % flux{5} - flux{7}; % 10 - sum(x{1:4})
 % conservation relation 2:
 % flux{9} - flux{8}; % 10 - sum(x{5:6})
 flux{8} - flux{9};  
 % conservation relation 3:          
 %flux{11}- flux{10}; % 10 - x{6}
 flux{10} - flux{11}                     %6
};

% calculate the Jacobian
for i=1:n
 for j=1:n
  Jf{i,j}=differentiate(f{i},x{j});
 endfor
endfor

% calculate K=∂f/∂ρ
for i=1:n
 for j=1:m
  Kf{i,j}=differentiate(f{i},theta{j});
 endfor
endfor

% calculate dJf/dx
dJf_x=cell(n,n,n);
for i=1:n
 for j=1:n
  for k=1:n
   dJf_x{i,j,k}=differentiate(Jf{i,j},x{k});
  endfor
 endfor
endfor

% calculate dJf/dρ
dJf_theta=cell(n,n,m);
for i=1:n
 for j=1:n
  for k=1:m
   dJf_theta{i,j,k}=differentiate(Jf{i,j},theta{k});
  endfor
 endfor
endfor

% calculate d²f/dρdρ
dKf_theta=cell(n,m,m);
for i=1:n
 for j=1:m
  for k=1:m
   dKf_theta{i,j,k}=differentiate(Kf{i,j},theta{k});
  endfor
 endfor
endfor

% print f,Jf,Kf and Sf to a file
if (exist("filename") && ischar(filename))
 fid=fopen(filename,"w");
else
 fid=stdout;
endif

% print f
fprintf(fid,"%% ode rhs f\nF.f=@(x,theta,u) [...\n");
for i=1:n
 fprintf(fid,"%s",to_char(f{i}));
 if (i<n)
  fprintf(fid,";...\n");
 endif
endfor
fprintf(fid,"];\n");

% print output matrix C
% output function: y = C(ρ) * x

fprintf(fid,"%% parameterised model output matrix\nF.C=");
fprintf(fid,"[...\n");
for i=1:rc
 fprintf(fid," %f ",C(i,:));
 if (i<rc)
  fprintf(fid,";...\n");
 endif
endfor
fprintf(fid,"];\n");


% print Jf
fprintf(fid,'%% Jacobian of f (df/dx)\nF.Jf=@(x,theta,u) [...\n');
for i=1:n
 for j=1:n
  fprintf(fid,"%s",to_char(Jf{i,j})); 
  if (j<n) 
   fprintf(fid,","); 
  endif
 endfor
 if (i<n) 
  fprintf(fid,";...\n"); 
 endif
endfor
fprintf(fid,'];\n');

% print S
fprintf(fid,'%% steady state sensitivity (dx(t→∞,ρ)/dρ)\nF.Sf=@(x,theta,u) -[...\n');
for i=1:n
 for j=1:n
  fprintf(fid,"%s",to_char(Jf{i,j})); 
  if (j<n) 
   fprintf(fid,","); 
  endif
 endfor
 if (i<n) 
  fprintf(fid,";...\n"); 
 endif
endfor
fprintf(fid,"]\\[...\n");
for i=1:n
 for j=1:m
  fprintf(fid,"%s",to_char(Kf{i,j})); 
  if (j<m) 
   fprintf(fid,","); 
  endif
 endfor
 if (i<n) 
  fprintf(fid,";...\n"); 
 endif
endfor
fprintf(fid,'];\n');

fprintf(fid,"%% model dimensions: number of state variables $n$; number of parameters $m$; number of outputs $l$;\nF.n=%i;\nF.m=%i;\nF.l=%i;\n",n,m,rc);
fprintf(fid,"%% initial conditions; either constant or function of parameters ρ and inputs u\nF.x0=%s;\n",x0);
fprintf(fid,"%% done.");


function_name=sprintf("sensitivity2_%s",model_name);
fprintf(fid,"%% link to second order sensitivity function\nF.dSf=@(x,theta,u,Jf,Sf) %s(x,theta,u,Jf,Sf);\n",function_name);
if (fid>2)
 fclose(fid);
endif

% print dS, more complicated function

filename=sprintf("%s.m",function_name);
if (exist("filename") && ischar(filename))
 fid=fopen(filename,"w");
else
 fid=stdout;
endif

fprintf(fid,"function [dS]=%s(x,theta,u,Jf,Sf)\n",function_name);

fprintf(fid,"m=%i;\nn=%i;\n",m,n);
fprintf(fid,"%%dJf_x\n");
 for l=1:n
  fprintf(fid,"dJf_x(:,:,%i)=[...\n",l);
  for i=1:n
   for j=1:n
    fprintf(fid,"%s",to_char(dJf_x{i,j,l})); 
    if (j<n) 
     fprintf(fid,","); 
    endif
   endfor
   if (i<n) 
    fprintf(fid,";...\n"); 
   endif
  endfor
 fprintf(fid,"];\n");
 endfor

fprintf(fid,"%%dJf_theta\n");
for k=1:m
 fprintf(fid,"dJf_theta(:,:,%i)=[...\n",k);
 for i=1:n
  for j=1:n
   fprintf(fid,"%s",to_char(dJf_theta{i,j,k})); 
   if (j<n) 
    fprintf(fid,","); 
   endif
  endfor
  if (i<n) 
   fprintf(fid,";...\n"); 
  endif
 endfor
 fprintf(fid,"];\n");
endfor

fprintf(fid,"%%dKf_theta\n");
for k=1:m
 fprintf(fid,"dKf_theta(:,:,%i)=[...\n",k);
 for i=1:n
  for j=1:m
   fprintf(fid,"%s",to_char(dKf_theta{i,j,k})); 
   if (j<m) 
    fprintf(fid,","); 
   endif
  endfor
  if (i<n) 
   fprintf(fid,";...\n"); 
  endif
 endfor
 fprintf(fid,"];\n");
endfor

fprintf(fid," A=zeros(n,n,m);\n");
fprintf(fid," for k=1:m\n");
fprintf(fid,"  for l=1:n\n");
fprintf(fid,"   A(:,l,k)=dJf_x(:,:,l)*Sf(:,k)+dJf_theta(:,l,k);\n");
fprintf(fid,"  end%%for\n");
fprintf(fid," end%%for\n");

fprintf(fid," B=zeros(n,m,m);\n");
fprintf(fid," for k=1:m\n");
fprintf(fid,"  for i=1:m\n");
fprintf(fid,"   B(:,i,k)=dJf_theta(:,:,i)*Sf(:,k);\n");
fprintf(fid,"  end%%for\n");
fprintf(fid," end%%for\n");


fprintf(fid," for k=1:m\n");
fprintf(fid,"  dS(:,:,k)=-Jf\\(A(:,:,k)*Sf+B(:,:,k)+dKf_theta(:,:,k));\n");
fprintf(fid," end%%for\n");
fprintf(fid,"end%%function\n");

if (fid>2)
 fclose(fid);
endif

