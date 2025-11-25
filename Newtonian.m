function flag = Newtonian(routing, ls, tau_tot, save_true)
    % 牛顿法求解 tendon 驱动连续体机器人的运动学模型
    % 输入参数：
    %   routing -  tendon 布线方式
    %   ls - 加载步数
    %   tau_tot - 总张力
    %   save_true - 是否保存结果（默认 true）
    
    if nargin < 4
        save_true = true;
    end

    % 添加工具函数路径
    addpath('./tools')

    % 独立参数定义
    [E, Gj, R_b, rho] = materials();
    grav = [0; 0; -9.81]; % 重力加速度
    [~, ~, ~, n_segments, n_tendons, l_j] = tendons(routing);

    % 边界条件
    r0 = zeros(n_segments, 3);
    Q0 = [rot2quat(eye(3))'; zeros(n_segments - 1, 4)];

    % 派生参数计算
    area = pi * R_b^2;
    I = pi * R_b^4 / 4;
    J = 2 * I;
    Kse = diag([Gj * area, Gj * area, E * area]);
    Kbt = diag([E * I, E * I, Gj * J]);

    % 封装共享参数为结构体（传递给子函数）
    Const = struct();
    Const.E = E;
    Const.Gj = Gj;
    Const.R_b = R_b;
    Const.rho = rho;
    Const.grav = grav;
    Const.routing = routing;
    Const.n_segments = n_segments;
    Const.n_tendons = n_tendons;
    Const.l_j = l_j;
    Const.area = area;
    Const.I = I;
    Const.J = J;
    Const.Kse = Kse;
    Const.Kbt = Kbt;
    Const.r0 = r0;
    Const.Q0 = Q0;

    % 打靶法求解初始化
    init_guess = zeros(6, 1);
    Y = cell(n_segments, 1); % 存储各段的 ODE 解
    flag = 1;
    t = 1;
    % 设置 fsolve 选项
    opt = optimoptions('fsolve', ...
        'Algorithm', 'levenberg-marquardt', ...
        'Display', 'off');

    tic;
    % 逐步加载张力并求解
    while flag > 0 && t <= ls
        Const.tau = t * tau_tot / ls; % 当前步张力
        % 调用外部子函数 shooting_fun，传递参数 Const 和 Y
        [init_guess, ~, flag] = fsolve(@(g) shooting_fun(g, Const, Y), init_guess, opt);
        t = t + 1;
    end
    time = toc;

    % 结果处理
    if flag <= 0
        warning('求解器未正常退出');
    else
        % 可视化
        l = sum(l_j);
        n_nodes = size(Y{1}, 1);
        Y = cell2mat(Y)';
        rX = Y(1:3, :);
        plot_TACR(rX, n_segments, l, n_nodes, 'Newtonian');
        
        % 保存结果
        if save_true
            QX = Y(4:7, :);
            rl_j = Y(1:3, end);
            Ql_j = Y(4:7, end);
            xi0 = Y(8:13, 1);
            save_TACR(routing, rl_j, Ql_j, rX, QX, xi0, time, Const.tau, ls);
        end
    end

end