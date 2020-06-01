function [S,I,R,t] = sirOL
%%%%%%%%%%%%%%%%%%%%%%%
% Open-loop SIR model %
%%%%%%%%%%%%%%%%%%%%%%%
global Tmax
global SInit
global IInit
global n

n = 2;
[t,sol] = ode45(@evolMat,[0 Tmax],[SInit; IInit; 0*ones(size(SInit))]);
S = sol(:,1:n);
I = sol(:,1+n:2*n);
R = sol(:,1+2*n:3*n);
end
function dydt = evolMat(t,y)
global Beta
global Gamma
global conf
global n
s = y(1:n); % Susceptible population
inF = y(1+n:2*n); % Infected population
r = y(1+2*n:3*n); % Recovered population
c = conf.co; % Open-loop Contact matrix
lam = (Beta).*(c*inF);
dydt = [-lam.*s; lam.*s - Gamma.*inF; Gamma.*inF];
end