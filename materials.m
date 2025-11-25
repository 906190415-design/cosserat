function [E, Gj, R_b, rho] = materials()
    E = 1e9; % [Pa] Young modulus (硅胶级别)
    nu = 0.3; % [-] Poisson's ratio
    Gj = E / (2 * (1 + nu)); % [Pa] shear modulus
    R_b = 1e-3; % [m] backbone radius (增至 1 mm)
    rho = 1000; % [kg/m^3] equivalent mass density
    % materials.m 中的定义（已存在）
    I = pi*R_b^4/4;   % 截面惯性矩（弯曲相关）
    J = 2*I;          % 极惯性矩（扭转相关）
end