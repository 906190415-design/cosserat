function dy = TACR_ODE(X, y, Const)  % 新增Const参数传递物理常数
    % 解包状态变量（扩展后）
    r = y(1:3);         % 位置
    Q = y(4:7);         % 姿态四元数
    K = y(8:10);        % 曲率
    G = y(11:13);       % 应变
    dr_dt = y(14:16);   % 线速度
    dQ_dt = y(17:20);   % 四元数一阶导数

    % 物理参数（从Const获取）
    R = quat2rot(Q);                % 旋转矩阵
    mass = Const.rho * Const.area * Const.l_j_seg;  % 当前段质量（段长×线密度）
    I_tensor = diag([Const.I, Const.I, Const.J]);   % 惯性张量（已在Newtonian.m中计算）

    % 1. 计算角速度ω（从四元数一阶导数转换）
    omega = quaternion_omega_from_q_dot(Q, dQ_dt);  % 工具函数：dQ/dt = 0.5·ω×Q → 反解ω

    % 2. 静力学载荷计算（复用原逻辑：肌腱力、重力等）
    a = zeros(3, 1);
    b = zeros(3, 1);
    A = zeros(3, 3);
    O = zeros(3, 3);
    H = zeros(3, 3);
    for jj = Const.j:n_segments  % Const.j为当前段索引
        for i = 1 : Const.n_tendons(jj)
            it = i + sum(Const.n_tendons(1:jj - 1));
            [D, dD, ddD] = tendons(Const.routing, Const.j, it, X);
            G_i = cross(K, D) + dD + G;
            G_i_norm = norm(G_i);
            A_i = -hat(G_i)^2 * (Const.tau(jj, i) / G_i_norm^3);
            O_i = -A_i * hat(D);
            a_i = A_i * (cross(K, G_i + dD) + ddD);
            a = a + a_i;
            b = b + cross(D, a_i);
            A = A + A_i;
            O = O + O_i;
            H = H + hat(D) * O_i;
        end
    end

    % 3. 静力学合力/合力矩（原rhs）
    N = Const.Kse*(G - [0; 0; 1]);
    C = Const.Kbt*K;
    static_force = -cross(K, C) - cross(G, N) - b;       % 静力合力
    static_torque = -cross(K, N) - R.'*Const.rho*Const.area*Const.grav - a;  % 静力合力矩

    % 4. 动力学平衡方程（加入惯性项）
    % 线加速度：静力合力 = 惯性力（mass·d_dr_dt）
    d_dr_dt = static_force / mass;  % 线加速度 = 合力 / 质量

    % 角加速度：静力合力矩 = 惯性力矩（I·dω/dt + ω×(I·ω)）
    inertia_coriolis = cross(omega, I_tensor * omega);  % 科氏项
    d_omega = inv(I_tensor) * (static_torque - inertia_coriolis);  % 角加速度

    % 5. 四元数二阶导数（d_dQ_dt）：从角加速度转换
    d_dQ_dt = 0.5 * (quaternion_dot_omega(d_omega) * Q + quaternion_dot_omega(omega) * dQ_dt);

    % 6. 状态导数打包（对应20维状态变量）
    dr = dr_dt;                  % 线速度（r的一阶导数）
    dQ = dQ_dt;                  % 四元数一阶导数（Q的一阶导数）
    dxi = [Const.Mat\rhs(1:3); Const.Mat\rhs(4:6)];  % 曲率/应变导数（复用原计算）
    dy = [dr; dQ; dxi; d_dr_dt; d_dQ_dt];  % 完整导数向量
end