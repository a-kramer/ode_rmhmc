MAPKplots
%theta prior normal
mu=[3.8713 0.9196];
sigma=2;

%bounds:
ax=[0,10,-4,4];
K=cell(1,NumOfSamples);
for k=1:NumOfSamples
 K{k}=mkde(Sample{k}');
endfor

nx1=128;
nx2=nx1;
x1=linspace(ax(1),ax(2),nx1);
x2=linspace(ax(3),ax(4),nx2);

X=cell(1,2);
[X{:}]=meshgrid(x1,x2);

P=cell(1,NumOfSamples+1);
[P{:}]=deal(zeros(nx1,nx2));

for i=1:nx1
 for j=1:nx2
  theta=[X{1}(i,j);X{2}(i,j)];

  %kde posteriors
  for l=1:NumOfSamples
   P{l}(i,j)+=K{l}(theta);
  endfor

  %exact posterior
  LL=0;
  for k=1:nU
   f=@(x) F.f(x,theta,U{k});
   Jf=@(x) F.Jf(x,theta,U{k});
   M={f,Jf};
   try
    xs=lsode(M,F.x0,[0,800])(:,2);
    LL-=0.5*sumsq((Data{1,k}-F.C*xs)/Data{2,k});
   catch
    LL-=Inf;
   end
  endfor
  %prior
  logprior=sum(log([normpdf(theta(1),mu(1),sigma),normpdf(theta(2),mu(2),sigma)]));
  P{NumOfSamples+1}(i,j)=exp(LL+logprior);
 endfor 
endfor


title('comparison of posteriors')
Title={'RMHMC','\NRRMHMC{}','\NRHMC{}','exact'};
h=zeros(1,3);
filename="MAPKposteriors";
for i=1:NumOfSamples
 figure(10+i);
 h(i)=pcolor(X{:},P{i});
 shading('interp');
% xlabel('$\theta_1$')
% ylabel('$\theta_2$')
 title(Title{i});
 axis(ax);
 set(gca,'ytick',[-4:2:4]);
 if i>1
  YL=cell(1,5);
  [YL{:}]=deal('');
  set(gca,'yticklabel',YL);
 endif
 set(gca,'tickdir','out');
 set(gcf,'papersize',[4,3]/2.54);
 set(gcf,'paperposition',[0,0,4,3]/2.54);
 set(gcf,'paperorientation','landscape');

 print('-depslatex','-color',sprintf("%s%i.tex",filename,i));
 system(sprintf("epstopdf %s%i.eps",filename,i));
 system(sprintf("mv %s%i.{pdf,tex} ~/Dropbox/OxfJrnl_Bioinf_steady_state_paper",filename,i))
endfor


title('comparison of posteriors')
filename="MAPKposteriorsBW";
h=zeros(1,NumOfSamples);
for i=1:NumOfSamples+1
 figure(20+i);
 [~,h(i)]=contour(X{:},P{i},5);
 %xlabel('$\theta_1$')
 %ylabel('$\theta_2$')
 title(Title{i});
 axis(ax);
 set(gca,'ytick',[-4:2:4]);
 %if i>1
  YL=cell(1,5);
  [YL{:}]=deal('');
  set(gca,'yticklabel',YL);
 %endif
 if i<NumOfSamples+1
  set(gca,'xtick',[0:2:8])
 endif
 %set(gca,'tickdir','out');
 set(gcf,'papersize',[5,4]/2.54);
 set(gcf,'paperposition',[0,0,5,4]/2.54);
 set(gcf,'paperorientation','landscape');
 print('-depslatex','-mono','-solid',sprintf("%s%i.tex",filename,i))
 system(sprintf("epstopdf %s%i.eps",filename,i));
 system(sprintf("mv %s%i.{pdf,tex} ~/Dropbox/OxfJrnl_Bioinf_steady_state_paper",filename,i));
endfor

