pkg load symbolic
symbols
n=6;
C=zeros(1,n);
C(1,5)= 13740.432; % taken from the Mifa.txt model file
rc=rows(C);
m=14;
s=1;

model_name="SBMifa_for_HMC"
filename=sprintf("%s.txt",model_name);

x=cell(n,1);
rho=cell(m,1);
Rho=cell(m,1);
u=cell(s,1);

x0={'10','0','0','0','0','0'};

% set up the state variables, sampling parameters and input parameters
for i=1:n
 x{i}=sym(sprintf('x_%02i',i));
endfor
for i=1:m
 theta{i}=sym(sprintf('theta_%02i',i));
 rho{i}=Exp(theta{i});
 Rho{i}=sym(sprintf('rho_%02i',i)); % literal
endfor
for i=1:s
 u{i}=sym(sprintf('u_%i',i));
endfor

dSFlux=cell(n,m,m);
for i=1:n
 for j=1:m
  for k=1:m
   dSFlux{i,j,k}=sym(sprintf("dSFlux_%02i_%02i_%02i",i,j,k)); % literal
  endfor
 endfor
endfor

% set up the fluxes
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
Flux=cell(1,nr);
for i=1:nr % define flux literals for later substitution here
 Flux{i}=sym(sprintf("flux_%02i",i)); %literal
endfor

%set up the right-hand-side of the ode (rhs)
f_literal={
% conservation relation 1:
 Flux{2} + Flux{7} + Flux{6} - Flux{1}; %1
 Flux{1} - Flux{2} - Flux{3};           %2
 Flux{3} - Flux{4} - Flux{6};           %3
 Flux{4} - Flux{5};                     %4
 % flux{5} - flux{7}; % 10 - sum(x{1:4})
 % conservation relation 2:
 % flux{9} - flux{8}; % 10 - sum(x{5:6})
 Flux{8} - Flux{9};  
 % conservation relation 3:          
 %flux{11}- flux{10}; % 10 - x{6}
 Flux{10} - Flux{11}                     %6
};
f=cell(1,n);
for i=1:n
 f{i}=subs(f_literal{i},Flux,flux);
 f{i}=subs(f{i},Rho,rho);
endfor


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

ch_rho=cell(1,m);
ch_Rho=cell(1,m);
for i=1:m
 ch_rho{i}=to_char(rho{i});
 ch_Rho{i}=to_char(Rho{i});
endfor

Flux_rho=cell(1,nr);
for i=1:nr
 Flux_rho{i}=to_char(flux{i});
 for a=1:m
  Flux_rho{i}=strrep(Flux_rho{i},ch_rho{a},ch_Rho{a});
 endfor
endfor

SFlux_rho=cell(n,m);
for i=1:n
 for j=1:m
  SFlux_rho{i,j}=to_char(dotS{i,j});
  for a=1:m
   SFlux_rho{i,j}=strrep(SFlux_rho{i,j},ch_rho{a},ch_Rho{a});
  endfor
 endfor
endfor

fprintf(fid,"********** MODEL NAME\n");
fprintf(fid,"%s\n",model_name);
fprintf(fid,"********** MODEL NOTES\n");
fprintf(fid,"exponential parameters\n");
fprintf(fid,"********** MODEL STATES\n");
for i=1:n
 fprintf(fid,"d/dt(%s) = %s\n",to_char(x{i}),to_char(f_literal{i}));
endfor
fprintf(fid,"\n");
for j=1:m
 for i=1:n
  fprintf(fid,"d/dt(%s) = %s\n",to_char(S{i,j}),SFlux_rho{i,j});
 endfor
endfor
fprintf(fid,"\n");
% initial conditions for state variables, sensitivities and second order sensitivities
for i=1:n
 fprintf(fid,"%s(0) = %s\n",to_char(x{i}),x0{i});
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
for i=1:m
 fprintf(fid,"%s = %s\n",ch_Rho{i},ch_rho{i});
endfor
fprintf(fid,"********** MODEL REACTIONS\n\n");
for i=1:nr
 fprintf(fid,"%s=%s\n",to_char(Flux{i}),Flux_rho{i});
endfor
fprintf(fid,"********** MODEL FUNCTIONS\n\n");
fprintf(fid,"********** MODEL EVENTS\n\n");
fprintf(fid,"********** MODEL MATLAB FUNCTIONS\n\n");

if (fid>2)
 fclose(fid);
endif

