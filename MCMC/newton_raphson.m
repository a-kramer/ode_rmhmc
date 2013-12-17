function [xs]=newton_raphson(F,X,rho,U)
%
% [xs]=newton_raphson(F,X,U)
%
% finds the steady states of ode system in F
% one steady state per entry in U{:} (row)
% X{:} (row) provides starting points, one per input
%
% F is a struct containing:
% .f: rhs of ode [function-handle]
%     f(x,rho,u)
% .Jf: Jacobian of f [function-handle]
% .Sf: Sensitivity of steady state [function-handle]
%      dx(rho,u)/drho
%
% returns xs{:} 
%      one steady state xs{i} per input U{i}
%

NumOfObs=length(U);
% tolerances for x and f
rtol=1e-6; 
atol=1e-6;

iterations=0;
maxIterations=40;

%fprintf('log(f(x)/x):\n');
for i=1:NumOfObs
    x=X{i};
    dx=x;
    while (any(abs(dx) > rtol*abs(x)+atol) && iterations < maxIterations)
%      Jf=F.Jf(x,rho,U{i});
%      if log10(rcond(Jf)) < -4
%          Jf=Jf+eye(F.n)*1e-6;
%      end%if
     dx=-F.Jf(x,rho,U{i})\F.f(x,rho,U{i});
     x=x+dx;
     iterations=iterations+1;
    end%while
    % we interpret F.f(x,rho,U{i}) as an error for the steady state.
    if any(abs(F.f(x,rho,U{i})) > rtol*abs(x) + atol)
        xs{i}=x;
        warning('newton_raphson: could not find zero to satisfactory precision. max(abs(F(x{%i})))=%g\n',i,max(abs(F.f(x,rho,U{i}))));
    else
        xs{i}=x;
    end%if
end%for

end%function

