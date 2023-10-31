%This code is made to see how biomass growth and D-lactate production
%changes with CO2 and acetate
clc
clear all

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;
model = changeObjective(model, {'BIOMASS_Ec_SynMixo','LDH_D'});
c = model.c; % Objective Coeficients

Aeq = (model.S);
beq = (model.b)';
%num_var = length(model.c);
num_var = 863;

rng default
func = @(v) -objective_infile(v, c);
options = optimoptions('ga','ConstraintTolerance',1e-6,'PlotFcn', @gaplotbestf,...
                       'CrossoverFrac',0.2,'PopulationSize',10,'StallGen',12,'Generations',3);
                        % 'CrossoverFrac',0.5,'PopulationSize',100,'StallGen',125,'Generations',100);
%options = optimoptions('ga','PlotFcns', @gaplotbestfun, 'PlotInterval', 10, 'PopInitRange', [-10 ; 10]);

x=0;
for i = -0.1:-0.1:-1;
    x=x+1;
    y=0;
    for j = -0.2:-0.2:-2;
        y=y+1;
    
        model = changeRxnBounds(model, {'CO2tex'}, i, 'b');
        model = changeRxnBounds(model, {'ACtex'}, j, 'b');

        % Call Genetic Algorithm
        FBAsolution = ga(func,num_var,[],[],Aeq,beq,[],[],[], options);

        biomass(x,y)=FBAsolution(687);
        lactate(x,y)=FBAsolution(323);
    end
end

CO2tex = (0.1:0.1:1);
Acetate = (0.2:0.2:2);
lactate;

subplot (1,2,1)
surf (CO2tex,Acetate,biomass)
title ('6803 autotrophic biomass per acetate and CO2 injection')
xlabel('CO2 consumption rate','fontweight','bold','fontsize',11)
ylabel('Acetate consumption rate','fontweight','bold','fontsize',11)
zlabel('Biomass production rate','fontweight','bold','fontsize',11)

subplot (1,2,2)
surf (CO2tex,Acetate,lactate)
title ('6803 D-lactate production per Acetate and CO2 injection')
xlabel('CO2 consumption rate','fontweight','bold','fontsize',11)
ylabel('Acetate consumption rate','fontweight','bold','fontsize',11)
zlabel('D-lactate production rate','fontweight','bold','fontsize',11)



function f = objective_infile(v, c)
    f = zeros(1,size(v,1));
    for i = 1:size(v,1)
        f(i) = v(i,:)*c;
    end

end

