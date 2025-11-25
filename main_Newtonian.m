% Simple example script to run a single simulation with the Newtonian
% approach

% Inputs
routing = 2; % 3 = single segment variable (helix)
tau_tot = @(t) [5*sin(t); 3*cos(t)]; % 时变肌腱张力函数
loading_steps = []; % 动力学无需加载步，移除或设为1

% to further modify the robot geometry and material properties, refer to
% functions/files (tendons.m and materials.m)

% to modify the numerical parameters of the simulation, refer to the
% function/file Newtonian.m itself.

% Run the simulation
Newtonian(routing, [], tau_tot); % 传入时变张