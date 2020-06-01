function daefun(daeh,x,~,u,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIR closed-loop dynamics ODEs %
% Ref: Eqn.1 in arXiv:XXXX.XXXX %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conf = daeh.userdata; % Structure with country specific parameters
cc = conf.cc;
co = conf.co;

c = co - u.F.*cc; % Closed-loop Control matrix

lam = (p.Beta).*(c*[x.I1; x.I2]);
wtI = lam.*[x.S1; x.S2] - (p.Gamma)*[x.I1; x.I2];

% SIR equation RHS expressions
sRHS = -lam.*[x.S1; x.S2];
iRHS = -sRHS - (p.Gamma)*[x.I1; x.I2];
rRHS = (p.Gamma)*[x.I1; x.I2];

daeh.setODE('S1', sRHS(1)); % Evolution of Susceptible population
daeh.setODE('S2', sRHS(2)); % Evolution of Susceptible population
daeh.setODE('I1', iRHS(1)); % Evolution of Infected population
daeh.setODE('I2', iRHS(2)); % Evolution of Infected population
daeh.setODE('R1', rRHS(1)); % Evolution of Recovered population
daeh.setODE('R2', rRHS(2)); % Evolution of Recovered population
daeh.setODE('It', (p.Py)*wtI(1) + (p.Po)*wtI(2)); % Evolution of probability weighted total infected
daeh.setODE('Time', 1);
end