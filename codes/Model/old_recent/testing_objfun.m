%Setting parameters, intial conditions and optmization boundaries
J_stack = [];
c = 0;
while c <8
[pars, init] = load_global_easy();

type_of_input = 5;
dt = 0.01;
settling_time = 11;
pars('dt') = dt;
pars('settling_time') = settling_time;
pars('type_of_input') = type_of_input;

%Importing experimental data
[time_exp, y_exp] = data_preprocessing(1);
t_exp = time_exp;

%Loading initial conditions
load('../simulations_saves/90sec_simulation.mat', 'x_vars');
init_values_loaded = x_vars(:, end);
for i = 1:length(init.keys)
    if init.keys{i} ~= "vO2" && init.keys{i} ~= "PAO2" && init.keys{i} ~= "P_1O2" && init.keys{i} ~= "P_2O2" && init.keys{i} ~= "P_3O2" && init.keys{i} ~= "P_4O2" && init.keys{i} ~= "P_5O2" && init.keys{i} ~= "MRtO2"
        init(init.keys{i})= init_values_loaded(i);
    end
end

list_of_pars = {"C2", "G_R_e_p", "KpCO2", "T0", "MRbCO2", "lambda1"};
        

for i = 1:length(list_of_pars)
    upper_boundry = pars(list_of_pars(i))*1.7;
    lower_boundry = pars(list_of_pars(i))*0.7;
    pars(list_of_pars(i)) = lower_boundry + (upper_boundry - lower_boundry) * rand();
    
    
end

J = obj_fun(pars,t_exp,y_exp, init, dt, settling_time, x_keys);
J_stack = [J_stack, J];
c = c +1;
end