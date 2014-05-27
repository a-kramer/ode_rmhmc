#!/usr/bin/octave -q

#samplefile={'SB_SMMALA_SBodeMAPKforSMMALA_1_1_Obs9_date_09-May-2014_time_21_37.mat',...
#            'SB_SMMALA_SBodeMmaForSMMALA_1_1_Obs7_date_09-May-2014_time_21_3.mat',...
#            'SB_SMMALA_SBMifaForSMMALA_1_1_Obs7_date_10-May-2014_time_0_0.mat';...
#            'steadystate_SMMALA_NRMAPK_smmala_1_1_Obs9_date_09-May-2014_time_21_39.mat',...
#            'steadystate_SMMALA_Mma_1_1_Obs7_date_09-May-2014_time_21_11.mat',...
#            'steadystate_SMMALA_Mifa_1_1_Obs7_date_10-May-2014_time_0_10.mat'};

samplefile={'SB_SMMALA_SBMifaForSMMALA_1_1_Obs7_date_13-May-2014_time_18_5.mat',...
            'SB_SMMALA_SBodeMmaForSMMALA_1_1_Obs7_date_12-May-2014_time_15_56.mat',...
            'SB_SMMALA_SBodeMAPKforSMMALA_1_1_Obs9_date_12-May-2014_time_15_33.mat';...
            'steadystate_SMMALA_Mifa_1_1_Obs7_date_14-May-2014_time_11_2.mat',...
            'steadystate_SMMALA_Mma_1_1_Obs7_date_12-May-2014_time_15_42.mat',...
            'steadystate_SMMALA_NRMAPK_smmala_1_1_Obs9_date_12-May-2014_time_15_28.mat' ...
};

% relative speed (ratio function error estimation)
r=@(v,u) [v(1)/u(1),v(2)*v(1)/(u(1)*v(1)) + u(2)*v(1)/u(1)^2];

% effective sampling speed
eff=cell(2,3);    
for i=1:2
 for j=1:3
   load(sprintf("./Results/%s",samplefile{i,j}));
   eff{i,j}=EffectiveSpeed(LLHistory,PosteriorTime);
 endfor
endfor

% relative speeds
R=cell(3,1);
for i=1:3
  R{i}=r(eff{2,i}.v,eff{1,i}.v);
endfor
rel_speeds=flipud(cat(1,R{:}));
h=errorbar([2,6,14],rel_speeds(:,1),rel_speeds(:,2))
xlabel('number of partameters');
ylabel('relative sampling speed $v_r$');
set(h,'linestyle','none');
set(h,'marker','.');
set(gca,'xtick',[2,6,14]);
set(gcf,'papersize',[9,8]/2.54);
set(gcf,'paperposition',[0,0,9,8]/2.54);
set(gcf,'paperorientation','landscape');
print -depslatexstandalone -mono figure4.tex
% plot
v={zeros(2,3),zeros(2,3)};
for j=1:3
 for i=1:2
  v{i}(:,j)=eff{i,j}.v';
 endfor
endfor

figure(1); clf;
for i=1:2
 h=errorbar(v{i}(1,:),v{i}(2,:)); hold on;
endfor
hold off
