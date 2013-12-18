pkg load symbolic
symbols

n=2;
C=zeros(1,n);
C(1,2)=1;
rc=rows(C); % output parameters
m=2;        % intrinsic parameters
s=1;

model_name="SB_MAPK_for_HMC"
filename=sprintf("%s.txt",model_name);

x=cell(n,1);
rho=cell(m,1);
u=cell(s,1);

% set up the state variables, sampling parameters and input parameters
for i=1:n
 x{i}=sym(sprintf('x_%02i',i));
endfor
for i=1:m
 theta{i}=sym(sprintf('theta_%02i',i));
 rho{i}=Exp(theta{i});
endfor
for i=1:s
 u{i}=sym(sprintf('u_%i',i));
endfor

% set up the fluxes

flux={rho{1}*u{1}/(1+u{1});rho{1}/(1+u{1})};
  
nr=length(flux);

%set up the right-hand-side of the ode (rhs)
f={flux{1} - (rho{2} + 2*flux{2})*x{1} + (rho{2} - flux{2})*x{2};
   flux{2}*x{1} - rho{2}*x{2}};


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


S=cell(n,m);
dotS=cell(n,m);
for i=1:n
 for j=1:m
  S{i,j}=sym(sprintf("S_%i_%i",i,j));
 endfor
endfor
for i=1:n
 for j=1:m
  dotS{i,j} = Kf{i,j};
  for k=1:n
   dotS{i,j} = dotS{i,j} + Jf{i,k}*S{k,j};
  endfor
 endfor
endfor

% print f,Jf,Kf and Sf to a file
if (exist("filename") && ischar(filename))
 fid=fopen(filename,"w");
else
 fid=stdout;
endif

fprintf(fid,"********** MODEL NAME\n");
fprintf(fid,"%s\n",model_name);
fprintf(fid,"********** MODEL NOTES\n");
fprintf(fid,"exponential parameters\n");
fprintf(fid,"********** MODEL STATES\n");
for i=1:n
 fprintf(fid,"d/dt(%s) = %s\n",to_char(x{i}),to_char(f{i}));
endfor
fprintf(fid,"\n");
for j=1:m
 for i=1:n
  fprintf(fid,"d/dt(%s) = %s\n",to_char(S{i,j}),to_char(dotS{i,j}));
 endfor
endfor
fprintf(fid,"\n");
% initial conditions for state variables, sensitivities and second order sensitivities
for i=1:n
 fprintf(fid,"%s(0) = 1.0\n",to_char(x{i}));
endfor 
for j=1:m
 for i=1:n
  fprintf(fid,"%s(0) = 0.0\n",to_char(S{i,j}));
 endfor
endfor 
fprintf(fid,"********** MODEL PARAMETERS\n");
for i=1:m
 fprintf(fid,"%s = 0.0\n",to_char(theta{i}));
endfor
for i=1:s
 fprintf(fid,"%s = 0.0\n",to_char(u{i}));
endfor

fprintf(fid,"********** MODEL VARIABLES\n\n");
fprintf(fid,"********** MODEL REACTIONS\n\n");
fprintf(fid,"********** MODEL FUNCTIONS\n\n");
fprintf(fid,"********** MODEL EVENTS\n\n");
fprintf(fid,"********** MODEL MATLAB FUNCTIONS\n\n");

if (fid>2)
 fclose(fid);
endif

