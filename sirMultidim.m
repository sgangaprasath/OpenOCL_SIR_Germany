function [solution,times,problem] = sirMultidim()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to evolve SIR Optimal control dynamics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Beta % Beta parameter in SIR model
global Gamma % Recovery rate in SIR model
global Tmax % Total time of simulation
global Np % Total population in millions
global SInit % Initial susceptible population
global IInit % Initial infected population
global conf % Structure with country specific parameters
global uM % Upper bound for value of controller
global Py % Probability of young infected requiring ICU
global Po % Probability of old infected requiring ICU
global AlpL % Importnace of Life cost
global AlpE % Importnace of Economic cost
global AlpS % Importnace of Social cost

problem = ocl.Problem([], @varsfun, @daefun, @pathcosts, 'N', 150, 'userdata', conf);
% @daefun - SIR closed-loop dynamics ODEs
% @varsfun - Definition of variables
% @pathcosts - Objective function definition
problem.setParameter('Beta', Beta);
problem.setParameter('Gamma', Gamma);
problem.setParameter('Py', Py);
problem.setParameter('Po', Po);
problem.setParameter('uM', uM);

problem.setParameter('Np', Np);
problem.setParameter('AlpL', AlpL);
problem.setParameter('AlpE', AlpE);
problem.setParameter('AlpS', AlpS);

% intial state bounds
problem.setInitialBounds('S1',SInit(1));
problem.setInitialBounds('S2',SInit(2));
problem.setInitialBounds('I1',IInit(1));
problem.setInitialBounds('I2',IInit(2));
problem.setInitialBounds('R1',0);
problem.setInitialBounds('R2',0);
problem.setInitialBounds('It',Py*IInit(1)+Po*IInit(2)); 
problem.setInitialBounds('Time',0);

problem.setEndBounds('Time', Tmax);

%   % Get and set initial guess
initialGuess = problem.getInitialGuess();

%   % Run solver to obtain solution
[solution,times] = problem.solve(initialGuess);
end