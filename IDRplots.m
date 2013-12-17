addpath(genpath('.'))
paper_folder="~/Dropbox/OxfJrnl_Bioinf_steady_state_paper";

files={'Results/ODE_RMHMC_NRMifa_1_1_Obs7_LFsteps_20_FPsteps_20_date_20-Nov-2013_time_4_2.mat';
       'Results/NRsensRefHMC_NRrefHMC_Mifa_1_1_Obs7_LFsteps_20_FPsteps_20_date_19-Nov-2013_time_19_12.mat';
       'Results/ODE_RMHMC_NRMma_1_1_Obs7_LFsteps_28_FPsteps_20_date_18-Nov-2013_time_17_41.mat';
       'Results/NRsensRefHMC_NRrefHMC_Mma_1_1_Obs7_LFsteps_28_FPsteps_20_date_19-Nov-2013_time_16_45.mat';
       'Results/ODE_RMHMC_SBodeMma_1_1_Obs7_LFsteps_28_FPsteps_20_date_25-Nov-2013_time_15_5.mat'};

NumOfSamples=length(files);
Sample=cell(1,NumOfSamples);
Speed=cell(1,NumOfSamples);
n=zeros(1,NumOfSamples);

for k=1:NumOfSamples
 load(files{k})
 N=length(LLHistory);
 [LL,dLL,ddLL,tau,dtau]=UWerr_np(LLHistory);
 n(k)=ceil(tau+dtau);
 Sample{k}=ParaHistory;
 Speed{k}=[N/(2*tau*PosteriorTime), dtau*N/(2*tau^2*PosteriorTime)];
endfor

% set models:
IM=cell(1,NumOfSamples);
for k=1:NumOfSamples
 if strfind(files{k},'Mifa')
  make_Mifa_model
 elseif strfind(files{k},'Mma')
  make_Mma_model
 else
  error("could not determine the used model for sample file: %s\n",files{k});
 endif
 IM{k}=F;
endfor


measurements=[18.3851,29.1153,26.3649,28.2924,42.3078,105.0286,94.9714];
sd=[20.3221,29.8399,31.7904,32.9355,16.5959,23.3895,23.3895];
U={0,0.01,0.1,0.3,1,10,100};
Ulabels=cellfun(@(x) sprintf('%g',x),U,"UniformOutput",false);
nU=length(U);

% relative speeds
r=@(v,u) [v(1)/u(1),v(2)*v(1)/(u(1)*v(1)) + u(2)*v(1)/u(1)^2];

z_file="IDRtrajectories5.octave";
z=cell(1,NumOfSamples);
if exist(z_file,"file")
 load(z_file);
else
 for k=1:NumOfSamples
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
    xs=lsode(M,IM{k}.x0,[0,800])(2,:)';
    z{k}(i,j)=IM{k}.C*xs;  
   end%for
  end%for
 end%for
 save("-binary",z_file,"z");
endif

# figure(1); clf; hold on;
# bar_filename="InsulinBar";
# q=[0.25,0.5,0.75];
# nq=length(q);
# Color=ones(3,3,nq);
# Color(:,:,nq-1)=[0.6, 0.0, 0.0;
#               0.0, 0.6, 0.0;
#               0.6, 0.0, 0.6];
# Color(:,:,nq)=[0.8, 0.0, 0.0;
#               0.0, 0.8, 0.0;
#               0.8, 0.0, 0.8];
# 
# for k=1:3
#  zq{k}=quantile(z{k},q);
#  for i=1:nq
#   Q{k,i}=zq{k}(i,:);
#  endfor
# endfor
# for i=nq:-1:1
#  bh=bar([1:nU],[measurements;cat(1,Q{:,i})]');
#  for k=2:4
#   set(bh(k),'facecolor',Color(k-1,:,i));
#  endfor
# endfor
# TF=2;
# th=text(0.7,measurements(1)*TF,'\bullet\,data');
# set(th,'rotation',45);
# th=text(0.9,zq{1}(3,1)*TF,'\bullet\,Mifa');
# set(th,'rotation',45);
# th=text(1.1,zq{2}(3,1)*TF,'\bullet\,Mma');
# set(th,'rotation',45);
# th=text(1.3,zq{3}(3,1)*TF,'\bullet\,Mma');
# set(th,'rotation',45);
# th=text(2.1,zq{2}(3,2)*TF,'\bullet\,\NRRMHMC');
# set(th,'rotation',45);
# th=text(2.3,zq{3}(3,2)*TF,'RMHMC');
# set(th,'rotation',45);
# ylabel('data and model output $z$ (sample)');
# xlabel('input $u$');
# set(gca,'xtick',[1:nU]);
# set(gca,'xticklabel',Ulabels);
# hold off;
# set(gcf,'papersize',[18,8]/2.54);
# set(gcf,'paperposition',[0,0,18,8]/2.54);
# set(gcf,'paperorientation','landscape');
# print('-depslatex','-color',sprintf("%s.tex",bar_filename));
# system(sprintf("epstopdf %s.eps",bar_filename));
# system(sprintf("mv %s.{pdf,tex} %s",bar_filename,paper_folder));

pkg load statistics
figure(2); clf; hold on;
box_filename="InsulinBox";
xshift=linspace(-0.4,0.4,4);
groupsize=NumOfSamples+2;
to_delete=zeros(3,1);
for k=1:NumOfSamples
 [s,hbp]=boxplot(z{k},0,'.');
 if isfield(hbp,'outliers') 
  set(hbp.outliers,'markersize',3);
 endif
 if isfield(hbp,'outliers2') 
  set(hbp.outliers2,'markersize',3);
 endif
  %struct [value,key]        =
 for [box_handle,box_element]=hbp % ↓ This is all ugly hacking, because boxplot 
                                  % does not support grouping of data, 
                                  % or x-positions of data for that matter.
   for l=1:length(box_handle)     % hbp contains lists of handles
    x=get(box_handle(l),'xdata');
    xi=round(x);
    set(box_handle(l),'xdata',(x-xi)+groupsize*(xi-1)+k);
   endfor
 endfor
endfor
set(gca,'xtick',groupsize*[1:nU]-groupsize+3);
set(gca,'xticklabel',Ulabels);
x_eb=groupsize*[1:nU]-groupsize+NumOfSamples+1;
h=errorbar(x_eb,measurements,sd);
%set(h,'marker','.');
hl=findobj(h, "type", "line");
for l=1:2 % fix errorbar plot (erase line)
 if any(isnan(get(hl(l),'xdata')))
  set(hl(l),'linewidth',2);
 else
  set(hl(l),'linestyle','none');
  set(hl(l),'marker','.');
 endif
endfor
xlabel('input $u$ [ins]');
ylabel('data [IRSP] and model outputs $z$');
axis([0 max(x_eb)+1 -7 300]);
hold off
set(gcf,'papersize',[18,8]/2.54);
set(gcf,'paperposition',[0,0,18,8]/2.54);
set(gcf,'paperorientation','landscape');
print('-depslatex','-mono',sprintf("%s.tex",box_filename));
system(sprintf("epstopdf %s.eps",box_filename));
system(sprintf("mv %s.{pdf,tex} %s",box_filename,paper_folder));
