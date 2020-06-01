function pathcosts(ch,x,~,u,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Objective function definition %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conf = ch.userdata;
Acost = conf.Acost;
Bcost = conf.Bcost;
% Sum of life cost, economic cost, social cost
ch.add((p.AlpL*sum(Acost*([x.I1; x.I2]))/p.Np + p.AlpE- p.AlpE*(1/p.Np)*((1-u.F)*sum(Bcost*([x.S1; x.S2] + [x.R1; x.R2]))) + p.AlpS*(u.F/p.uM).^2));
end