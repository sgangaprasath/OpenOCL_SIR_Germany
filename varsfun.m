function varsfun(svh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to define variables for OpenOCL %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Np
global Icus
global uM

svh.addState('S1','lb', 0); % Susceptible population
svh.addState('S2','lb', 0); % Susceptible population
svh.addState('I1','lb', 0); % Infected population
svh.addState('I2','lb', 0); % Infected population
svh.addState('R1','lb', 0); % Recovered population
svh.addState('R2','lb', 0); % Recovered population
svh.addState('It', 'lb', 0, 'ub', Icus*Np); % Probability weighted infected population with upper bound

% Scalar u: 0 <= F <= uub
svh.addControl('F', 'lb', 0, 'ub', uM); % Controller with upper bound
svh.addParameter('Beta');
svh.addParameter('Gamma');
svh.addParameter('Py');
svh.addParameter('Po');
svh.addParameter('uM');
svh.addParameter('Icus');
svh.addState('Time');

svh.addParameter('Np');
svh.addParameter('AlpL');
svh.addParameter('AlpE');
svh.addParameter('AlpS');
svh.addParameter('Disc');
end