pkg load symbolic
symbols

n=3;
C=zeros(1,n);
C(1,3)=1;
rc=rows(C);
m=7;
l=1;

%Output Parameter Indices
o=[7];



model_name="NRInsulin"
filename=sprintf("make_%s_model.m",model_name);

x=cell(n,1);
rho=cell(m,1);
u=cell(l,1);

x0="[10;0;0]";

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

% set up the fluxes

flux={rho{1}*x{1},rho{2}*x{1}*u{1},rho{3}*(10-x{1}-x{2}),rho{4}*x{2},rho{5}*x{2}*(10-x{3}),rho{6}*x{3}};
  
nr=length(flux);

%set up the right-hand-side of the ode (rhs)
f={flux{3}-flux{1}-flux{2},flux{1}+flux{2}-flux{4},flux{5}-flux{6}};


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
% print output function:
fprintf(fid:'%% output function\n');
fprintf(fid,'F.h = @(x,theta) diag([');
for i=1:l
 fprintf(fid,' %s ',theta{o(i)});
endfor
fprintf(fid,'])*F.C*x');
% print output sensitivity
fprintf(fid,'%% output state sensitivity (dh(t→∞,ρ)/dρ)\nF.Sh=@(CSf,Cx) CSf+cat(2,zeros(%i,%i),diag(Cx));\n',l,m-l);


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

