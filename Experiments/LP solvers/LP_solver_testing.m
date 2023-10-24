%Example_1
%This is for adding a model to system and see how it works
clc
clear

%Load Model (Cyano model)
load('iJN678.mat')
model = iJN678;


% using mosek solver
% test running solver ten times
elapsed_time_mosek = 0;
changeCobraSolver('mosek', 'LP', 1, 1);

for e = 1:10
    mosek_obj_values = zeros(10);
    x=0;
    
    tic;
    for i = -0.1:-0.1:-1
        x = x + 1;
        y = 0;
        for j = -0.2:-0.2:-2
            y = y + 1;
    
            model = changeRxnBounds(model, {'CO2tex'}, i, 'b');
            model = changeRxnBounds(model, {'ACtex'}, j, 'b');
    
            model = changeObjective(model, {'BIOMASS_Ec_SynMixo', 'LDH_D'});
            FBAsolution_mosek = optimizeCbModel(model, 'max');
            %mosek_obj_values(x, y) = FBAsolution.f;
    
            biomass_mosek(x, y) = FBAsolution_mosek.x(687);
            lactate_mosek(x, y) = FBAsolution_mosek.f;
    
        end
    end
    elapsed_time_mosek = elapsed_time_mosek + toc;
end

% using glpk solver
% test running solver ten times
elapsed_time_glpk = 0;
changeCobraSolver('glpk', 'LP', 1, 1);

for e = 1:10
    glpk_obj_values = zeros(10);
    x=0;
    
    tic;
    for i = -0.1:-0.1:-1
        x = x + 1;
        y = 0;
        for j = -0.2:-0.2:-2
            y = y + 1;
    
            model = changeRxnBounds(model, {'CO2tex'}, i, 'b');
            model = changeRxnBounds(model, {'ACtex'}, j, 'b');
    
            model = changeObjective(model, {'BIOMASS_Ec_SynMixo', 'LDH_D'});
            FBAsolution_glpk = optimizeCbModel(model, 'max');
            %glpk_obj_values(x, y) = FBAsolution.f;
    
            biomass_glpk(x, y) = FBAsolution_glpk.x(687);
            lactate_glpk(x, y) = FBAsolution_glpk.f;
    
        end
    end
    elapsed_time_glpk = elapsed_time_glpk + toc;
end
disp("mosek avg runtime: " + elapsed_time_mosek/10);
disp("glpk avg runtime: " + elapsed_time_glpk/10);


% graph displays
CO2tex = (0.1:0.1:1);
Acetate = (0.2:0.2:2);
lactate_mosek;
biomass_mosek;

subplot(1,2,1)
surf (CO2tex, Acetate, biomass_mosek)
title('6803 autotrophic biomass per Acetate and CO2 injection')
xlabel('CO2 consumption rate', 'fontweight', 'bold', 'fontsize', 11)
xlabel('Acetate consumption rate', 'fontweight', 'bold', 'fontsize', 11)
xlabel('biomass production rate', 'fontweight', 'bold', 'fontsize', 11)

subplot(1,2,2)
surf (CO2tex, Acetate, lactate_mosek)
title('6803 D-lactate production per Acetate and CO2 injection')
xlabel('CO2 consumption rate', 'fontweight', 'bold', 'fontsize', 11)
xlabel('Acetate consumption rate', 'fontweight', 'bold', 'fontsize', 11)
xlabel('D-lactate production rate', 'fontweight', 'bold', 'fontsize', 11)

