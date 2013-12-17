function [eff]=EffectiveSpeed(LLHistory,PosteriorTime)
N=size(LLHistory,1)
LL=cell(1,3);
tau_int_LL=cell(1,2);
[LL{:},tau_int_LL{:}]=UWerr(LLHistory);
eff.v=zeros(1,2);
eff.v(1)=N/(2*tau_int_LL{1}*PosteriorTime);
eff.v(2)=tau_int_LL{2}*eff.v(1)/tau_int_LL{1};
eff.tau_int=tau_int_LL;
eff.SampleSize=zeros(1,2);
eff.SampleSize(1)=N/(2*tau_int_LL{1});
eff.SampleSize(2)=tau_int_LL{2}*eff.SampleSize(1)/tau_int_LL{1};
end