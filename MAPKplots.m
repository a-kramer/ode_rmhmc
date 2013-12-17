addpath(genpath('.'))
files={'Results/ODE_RMHMC_SBodeMAPK_1_1_9DPS_735532_2  21  29_4538.mat';
       'Results/ODE_RMHMC_expMAPK_1_1_9DPS_735531_19  36  31_2742.mat';       
       'Results/NRsensRefHMC_NRrefHMC_MAPK_1_1_Obs9_LFsteps_10_FPsteps_10_date_19-Nov-2013_time_16_5.mat'};


NumOfSamples=length(files);
Sample=cell(1,NumOfSamples);
Speed=cell(1,NumOfSamples);
n=zeros(1,NumOfSamples);

for k=1:NumOfSamples
 load(files{k});
 N=length(LLHistory);
 [LL,dLL,ddLL,tau,dtau]=UWerr_np(LLHistory);
 n(k)=ceil(tau+dtau);
 Sample{k}=ParaHistory;
 Speed{k}=[N/(2*tau*PosteriorTime), dtau*N/(2*tau^2*PosteriorTime)];
endfor

make_expMAPK_model;

U=fliplr({0.92,0.751,0.633,0.359,0.389,0.256,0.197,0.194,0.097});
NumOfObs=length(U);
Data=cell(2,NumOfObs);
measurements=fliplr({0.993,0.789,0.646,0.737,0.798,0.487,0.6,0.689,0.043});
[Data{1,:}]=deal(measurements{:});
[Data{2,:}]=deal(0.2);

nU=length(U);

I=1:ceil(mean(n)):N;
nI=length(I);

nu=64;
u=linspace(0,1,nu);


z_file="MAPKoutputs.octave";
if exist(z_file,"file")
 load(z_file);
else
 z=zeros(nI,nu);
 for i=1:nI
  for j=1:nu
   theta=ParaHistory(I(i),:);
   f=@(x) F.f(x,theta,u(j));
   Jf=@(x) F.Jf(x,theta,u(j));
   M={f,Jf};
   xs=lsode(M,F.x0,[0,800])(2,:)';
   z(i,j)=F.C*xs;  
  end%for
 end%for
 save("-binary",z_file,"z");
endif