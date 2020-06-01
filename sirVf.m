function [Svf,Ivf,Sdotvf,Idotvf] = sirVf(sFinVf,iFinVf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SIR model for visualisation%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Beta
global Gamma
global conf
global Np

c = conf.co;
lam = Beta*c*iFinVf;
Sdotvf = -lam.*sFinVf/Np;
Idotvf = (lam.*sFinVf - Gamma*iFinVf)/Np;
Svf = sFinVf/Np;
Ivf = iFinVf/Np;
end
