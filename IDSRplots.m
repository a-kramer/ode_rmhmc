Sample=cell(1,3);
n=zeros(1,3);
load('Results/ODE_RMHMC_NRMifa_1_1_7DPS_735531_21  10  43_6235.mat')
[LL,dLL,ddLL,tau,dtau]=UWerr_np(LLHistory);
n(1)=ceil(tau+dtau);
Sample{1}=ParaHistory;

load('Results/ODE_RMHMC_NRMma_1_1_7DPS_735544_16  25   4_6450.mat')
[LL,dLL,ddLL,tau,dtau]=UWerr_np(LLHistory);
n(2)=ceil(tau+dtau);
Sample{2}=ParaHistory;

load('Results/ODE_RMHMC_SBodeInsulinFast_1_1_6DPS_735536_21  10   3_5495')
[LL,dLL,ddLL,tau,dtau]=UWerr_np(LLHistory);
n(3)=ceil(tau+dtau);
Sample{3}=ParaHistory;

IM=cell(1,2);
make_Mifa_model
IM{1}=F;
make_Mma_model
IM{2}=F;
IM{3}=F;

measurements=[18.3851,29.1153,26.3649,28.2924,42.3078,105.0286,94.9714];
sd=[20.3221,29.8399,31.7904,32.9355,16.5959,23.3895,23.3895];
U={0,0.01,0.1,0.3,1,10,100};
nU=length(U);

z_file="IDR.octave";
z=cell(1,3);
if exist(z_file,"file")
 load(z_file);
else
 for k=1:3
  N=rows(Sample{k})
  I=1:n(k):N;
  nI=length(I);
  z{k}=zeros(nI,nU);
  for i=1:nI
   for j=1:nU
    theta=Sample{k}(I(i),:);
    f=@(x) IM{k}.f(x,theta,U{j});
    Jf=@(x) IM{k}.Jf(x,theta,U{j});
    M={f,Jf};
    xs=lsode(M,IM{k}.x0,[0,800])(:,2);
    z{k}(i,j)=F.C*xs;  
   end%for
  end%for
 end%for
 save("-binary",z_file,"z");
endif