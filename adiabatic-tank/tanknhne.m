clear; clc; close all;

%load dat
load n2o.mat;

T_data     = Temperature_C_ + 273.15;
Psat_data  = Pressure_bar_ * 1e5;
rho_l_data = Density_l_Kg_m3_;
rho_v_data = Density_v_Kg_m3_;
h_l_data   = Enthalpy_l_KJ_kg_ * 1e3;
h_v_data   = Enthalpy_v_KJ_kg_ * 1e3;
s_l_data   = Entropy_l_J_g_K_ * 1e3;
s_v_data   = Entropy_v_J_g_K_ * 1e3;
Cv_l_data  = Cv_l_J_g_K_ * 1e3;
Cv_v_data  = Cv_v_J_g_K_ * 1e3;
Lv_data    = h_v_data - h_l_data;

T_fit_min = min(T_data);
T_fit_max = max(T_data);

% Polynomial degree for density fits (used in ODE derivatives)
n_rho = 7;

% Fit polynomials to rho_l and rho_v for smooth derivatives
[p_rho_l_ode, ~, mu_rho_l_ode] = polyfit(T_data, rho_l_data, n_rho);
[p_rho_v_ode, ~, mu_rho_v_ode] = polyfit(T_data, rho_v_data, n_rho);

% Analytic derivative coefficients
% Use polyval with mu for evaluation; derivative requires chain rule: d/dT = d/dT_norm * (1/sigma)
sigma_rho_l = mu_rho_l_ode(2);
sigma_rho_v = mu_rho_v_ode(2);
dp_rho_l    = polyder(p_rho_l_ode);
dp_rho_v    = polyder(p_rho_v_ode);

rho_l_fn  = @(T) polyval(p_rho_l_ode, T, [], mu_rho_l_ode);
rho_v_fn  = @(T) polyval(p_rho_v_ode, T, [], mu_rho_v_ode);
drho_l_fn = @(T) polyval(dp_rho_l, (T - mu_rho_l_ode(1)) / sigma_rho_l) / sigma_rho_l;
drho_v_fn = @(T) polyval(dp_rho_v, (T - mu_rho_v_ode(1)) / sigma_rho_v) / sigma_rho_v;

% pchip interpolation for all other thermodynamic properties
Cv_l_fn  = @(T) interp1(T_data, Cv_l_data,  T, 'pchip');
Cv_v_fn  = @(T) interp1(T_data, Cv_v_data,  T, 'pchip');
Psat_fn  = @(T) interp1(T_data, Psat_data,  T, 'pchip');
Lv_fn    = @(T) interp1(T_data, Lv_data,    T, 'pchip');
h_l_fn   = @(T) interp1(T_data, h_l_data,   T, 'pchip');
h_v_fn   = @(T) interp1(T_data, h_v_data,   T, 'pchip');
s_l_fn   = @(T) interp1(T_data, s_l_data,   T, 'pchip');
s_v_fn   = @(T) interp1(T_data, s_v_data,   T, 'pchip');

T_check = 280.15;
fprintf('=== Fit at %.1f K (%.1f°C) ===\n', T_check, T_check-273.15);
fprintf('rho_l = %.2f kg/m^3\n', rho_l_fn(T_check));
fprintf('rho_v = %.2f kg/m^3\n', rho_v_fn(T_check));
fprintf('Psat  = %.3f bar\n',    Psat_fn(T_check)/1e5);
fprintf('h_l   = %.2f kJ/kg\n',  h_l_fn(T_check)/1e3);
fprintf('s_l   = %.4f J/kg/K\n', s_l_fn(T_check));

D_tank = 5.25  * 0.0254;            % m — inner diameter
H_tank = 17.64 * 0.0254;            % m — height
V0     = pi/4 * D_tank^2 * H_tank;  % m^3 — internal volume
C_d       = 0.65;
kappa     = 2;
alpha     = 1/(1+kappa);
beta      = kappa/(1+kappa);
A_inj     = 16e-6;
T0        = 10 + 273.15;
m_N2O     = 3.0;
P_chamber = 1.013e5;
t_burn    = 6;

P_tank0 = Psat_fn(T0);
rl0     = rho_l_fn(T0);
rv0     = rho_v_fn(T0);
m_l0    = (V0 - m_N2O/rv0) / (1/rl0 - 1/rv0);

T_star   = T0;
rho_star = rho_l_fn(T_star);
P_star   = Psat_fn(T_star);
M_star   = m_N2O;
h_star   = h_l_fn(T_star);

% Set P_chamber and compute T_chamber
T_sat_chamber = fzero(@(T) Psat_fn(T) - P_chamber, [T_fit_min, T0]);

DeltaP0   = max(P_tank0 - P_chamber, 0);
mdot_SPI0 = C_d * A_inj * sqrt(2 * rho_l_fn(T0) * DeltaP0);

s_tank0 = s_l_fn(T0);
[mdot_HEM0, x_exit0, rho_2ph0, delta_h0] = hem_model(T0, ...
    s_tank0, C_d, A_inj, T_sat_chamber, ...
    s_l_fn, s_v_fn, rho_l_fn, rho_v_fn, h_l_fn, h_v_fn);

mdot_nom = alpha * mdot_SPI0 + beta * mdot_HEM0;
t_star   = M_star / mdot_nom;
V_star   = V0 * rho_star / M_star;
tau_burn = t_burn / t_star;

fprintf('\n=== Initial Conditions ===\n');
fprintf('T0 = %.1f°C\n',         T0-273.15);
fprintf('P_tank0 = %.2f bar\n',  P_tank0/1e5);
fprintf('P_chamber = %.2f bar\n',P_chamber/1e5);
fprintf('ΔP0 = %.2f bar\n',      DeltaP0/1e5);
fprintf('T_sat_chamber = %.2f K (%.1f°C)\n', T_sat_chamber, T_sat_chamber-273.15);
fprintf('SPI0 = %.4f kg/s\n',    mdot_SPI0);
fprintf('HEM0 = %.4f kg/s (x=%.3f, ρ_2ph=%.1f kg/m³, Δh=%.2f kJ/kg)\n', ...
    mdot_HEM0, x_exit0, rho_2ph0, delta_h0/1000);
fprintf('NHNE mdot0 = %.4f kg/s\n', mdot_nom);
fprintf('t* = %.2f s\n',         t_star);
fprintf('tau_burn = %.4f\n',     tau_burn);

params = struct();
params.V_star         = V_star;
params.P_chamber      = P_chamber;
params.P_star         = P_star;
params.T_star         = T_star;
params.rho_star       = rho_star;
params.h_star         = h_star;
params.M_star         = M_star;
params.t_star         = t_star;
params.tau_burn       = tau_burn;
params.T_sat_chamber  = T_sat_chamber;  
params.varrho_l   = @(Th) rho_l_fn(Th*T_star) / rho_star;
params.varrho_v   = @(Th) rho_v_fn(Th*T_star) / rho_star;
params.dvarrho_l  = @(Th) drho_l_fn(Th*T_star) * T_star / rho_star;
params.dvarrho_v  = @(Th) drho_v_fn(Th*T_star) * T_star / rho_star;
params.Omega_fn   = @(Th) Cv_v_fn(Th*T_star) ./ Cv_l_fn(Th*T_star);
params.Lambda_fn  = @(Th) Lv_fn(Th*T_star) ./ (Cv_l_fn(Th*T_star) .* T_star);
params.Psi_fn     = @(Th) P_star ./ (Cv_l_fn(Th*T_star) .* T_star .* rho_star);
params.Pi_fn      = @(Th) Psat_fn(Th*T_star) / P_star;
params.A_inj      = A_inj;
params.C_d        = C_d;
params.alpha      = alpha;
params.beta       = beta;
params.mu_l_min   = 0.02;
params.rho_l_fn   = rho_l_fn;
params.rho_v_fn   = rho_v_fn;
params.h_l_fn     = h_l_fn;
params.h_v_fn     = h_v_fn;
params.s_l_fn     = s_l_fn;
params.s_v_fn     = s_v_fn;
params.Psat_fn    = Psat_fn;

xi0  = [m_l0/M_star; T0/T_star];
opts = odeset('Events', @(tau,xi) stop_event(tau, xi, params), ...
              'RelTol', 1e-7, 'AbsTol', 1e-9);

[tau_sol, xi_sol] = ode45(@(tau,xi) ode_rhs(tau, xi, params), ...
                          [0, 10*tau_burn], xi0, opts);

mu_l  = xi_sol(:,1);
Theta = xi_sol(:,2);
t     = tau_sol * t_star;
T_sol = Theta * T_star;
P_sol = arrayfun(@(Th) Psat_fn(Th*T_star)/1e5, Theta);

varrho_l = params.varrho_l;
varrho_v = params.varrho_v;
mu_v     = arrayfun(@(ml,Th) varrho_v(Th)*(V_star - ml/varrho_l(Th)), mu_l, Theta);
m_l_sol  = mu_l * M_star;
m_v_sol  = mu_v * M_star;
m_total  = m_l_sol + m_v_sol;

[mdot_out, mdot_SPI_sol, mdot_HEM_sol] = arrayfun(@(Th) ...
    compute_mdot(Th*T_star, params), Theta);

% Evaporation rate: re-evaluate RHS at each solution point to avoid
% finite-difference noise on the irregularly-spaced ode45 output grid.
mdot_evap = arrayfun(@(ml, Th) get_mdot_evap(ml, Th, params), mu_l, Theta);

fprintf('\n=== Results ===\n');
fprintf('Burn time = %.2f s\n',       t(end));
fprintf('Final mass = %.3f kg\n',     m_total(end));
fprintf('Mass expelled = %.3f kg\n',  m_N2O - m_total(end));
fprintf('Average mdot = %.4f kg/s\n', (m_N2O - m_total(end))/t(end));
fprintf('Final T = %.1f°C\n',         T_sol(end)-273.15);
fprintf('Final P = %.2f bar\n',       P_sol(end));

%% ========== PLOTS ==========

% ---- Figure 1 ----
n = 5;
% Reuse degree-7 poly fits already computed for ODE density derivatives
p_rho_l  = p_rho_l_ode;  mu_rho_l = mu_rho_l_ode;
p_rho_v  = p_rho_v_ode;  mu_rho_v = mu_rho_v_ode;
[p_Cv_l,  ~, mu_Cv_l]  = polyfit(T_data, Cv_l_data,  n);
[p_Cv_v,  ~, mu_Cv_v]  = polyfit(T_data, Cv_v_data,  n);
[p_Psat,  ~, mu_Psat]  = polyfit(T_data, Psat_data,  n);
[p_Lv,    ~, mu_Lv]    = polyfit(T_data, Lv_data,    n);

T_plt = linspace(min(T_data), max(T_data), 300)';

figure('Name','NIST vs Polynomial Fits','Position',[50 50 1020 580]);
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

nexttile; hold on;
plot(T_data, rho_l_data, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_rho_l, T_plt, [], mu_rho_l), 'b-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('\rho_l (kg/m^3)'); title('Liquid density');
grid on; box on; legend('NIST','Fit','Location','best');

nexttile; hold on;
plot(T_data, rho_v_data, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_rho_v, T_plt, [], mu_rho_v), 'r-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('\rho_v (kg/m^3)'); title('Vapour density');
grid on; box on; legend('NIST','Fit','Location','best');

nexttile; hold on;
plot(T_data, Cv_l_data, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_Cv_l, T_plt, [], mu_Cv_l), 'b-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('c_{v,l} (J/kg/K)'); title('Liquid c_v');
grid on; box on; legend('NIST','Fit','Location','best');

nexttile; hold on;
plot(T_data, Cv_v_data, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_Cv_v, T_plt, [], mu_Cv_v), 'r-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('c_{v,v} (J/kg/K)'); title('Vapour c_v');
grid on; box on; legend('NIST','Fit','Location','best');

nexttile; hold on;
plot(T_data, Psat_data/1e5, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_Psat, T_plt, [], mu_Psat)/1e5, 'g-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('P_{sat} (bar)'); title('Saturation pressure');
grid on; box on; legend('NIST','Fit','Location','best');

nexttile; hold on;
plot(T_data, Lv_data/1e3, 'o', 'Color',[0.5 0.5 0.5], 'MarkerSize',3);
plot(T_plt, polyval(p_Lv, T_plt, [], mu_Lv)/1e3, 'm-', 'LineWidth',1.5);
xlabel('T (K)'); ylabel('L_v (kJ/kg)'); title('Latent heat');
grid on; box on; legend('NIST','Fit','Location','best');

% ---- Figure 2 ----
% Compute non-dimensional groups along solution
Omega_sol  = arrayfun(params.Omega_fn,  Theta);
Lambda_sol = arrayfun(params.Lambda_fn, Theta);
Psi_sol    = arrayfun(params.Psi_fn,    Theta);
Pi_sol     = arrayfun(params.Pi_fn,     Theta);

% Normalised groups (divide by initial value)
Omega_norm  = Omega_sol  / Omega_sol(1);
Lambda_norm = Lambda_sol / Lambda_sol(1);
Psi_norm    = Psi_sol    / Psi_sol(1);

figure('Name','Additional Diagnostics','Position',[100 100 1200 800]);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

% (a) P-T trajectory
nexttile; hold on;
plot(T_sol-273.15, P_sol, 'b-', 'LineWidth',1.5);
plot(T_data-273.15, Psat_data/1e5, 'k:', 'LineWidth',1.2);
xlabel('Temperature (°C)'); ylabel('Pressure (bar)');
title('P-T'); grid on; box on;
legend('Model','NIST P_{sat}','Location','best');

% (b) Normalised non-dimensional groups
nexttile; hold on;
plot(t, Lambda_norm, 'b-', 'LineWidth',1.5);
plot(t, Omega_norm,  'r--', 'LineWidth',1.5);
plot(t, Psi_norm,    'g-.', 'LineWidth',1.5);
xlabel('Time (s)'); ylabel('X(\tau)/X(0)');
title('Non-dimensionalised groups (normalised)'); grid on; box on;
legend('\Lambda/\Lambda_0','\Omega/\Omega_0','\Psi/\Psi_0','Location','best');

% (c) Non-dimensional state vector
nexttile; hold on;
plot(tau_sol, mu_l,  'b-', 'LineWidth',1.5);
plot(tau_sol, Theta, 'r--', 'LineWidth',1.5);
xline(1,'k:','LineWidth',1);
xlabel('\tau = t/t^*'); ylabel('State');
title('Non-dimensional state vector'); grid on; box on;
legend('\mu_l','\Theta','Location','best');

% (d) Non-dimensional pressure Pi
nexttile;
plot(tau_sol, Pi_sol, 'm-', 'LineWidth',1.5);
xlabel('\tau = t/t^*'); ylabel('\Pi = P/P^*');
title('Non-dimensional pressure \Pi(\Theta)'); grid on; box on;

% ---- Figure 3 ----
figure('Name','Blowdown','Position',[100 100 1400 900]);
subplot(2,3,1);
plot(t, P_sol, 'b-', 'LineWidth', 1.5); hold on;
xlim([0 7]);
yline(P_chamber/1e5, 'r--');
xlabel('Time (s)'); ylabel('Pressure (bar)');
title('Tank Pressure'); grid on;
legend('Tank','Chamber');

subplot(2,3,2);
plot(t, T_sol-273.15, 'r-', 'LineWidth', 1.5);
xlim([0 7]);
xlabel('Time (s)'); ylabel('Temperature (°C)');
title('Fluid Temperature'); grid on;

subplot(2,3,3);
plot(t, [m_l_sol, m_v_sol, m_total], 'LineWidth', 1.5);
xlim([0 7]);
xlabel('Time (s)'); ylabel('Mass (kg)');
legend('Liquid','Vapour','Total');
title('Phase Masses'); grid on;

subplot(2,3,4);
fill_level = (m_l_sol./rho_l_fn(T_sol))/V0*100;
plot(t, fill_level, 'g-', 'LineWidth', 1.5);
xlim([0 7]);
xlabel('Time (s)'); ylabel('Liquid Fill (%)');
title('Liquid Volume Fraction'); grid on;

subplot(2,3,5);
plot(t, [mdot_out, mdot_SPI_sol, mdot_HEM_sol], 'LineWidth', 1.5);
xlim([0 7]);
xlabel('Time (s)'); ylabel('Mass Flow Rate (kg/s)');
legend('NHNE','SPI','HEM');
title('Mass Flow Rates'); grid on;

subplot(2,3,6);
plot(t, mdot_evap, 'm-', 'LineWidth', 1.5);
xlim([0 7]);
xlabel('Time (s)'); ylabel('Evaporation Rate (kg/s)');
title('Evaporation Rate'); grid on;

%% Functions

function [mdot_HEM, x_exit, rho_2ph, delta_h] = hem_model(T_tank, ...
        s_tank, C_d, A_inj, T_sat_chamber, ...
        s_l_fn, s_v_fn, rho_l_fn, rho_v_fn, h_l_fn, h_v_fn)
    s_l_exit = s_l_fn(T_sat_chamber);
    s_v_exit = s_v_fn(T_sat_chamber);
    x_exit   = (s_tank - s_l_exit) / (s_v_exit - s_l_exit);
    x_exit   = max(0, min(1, x_exit));
    rho_2ph  = 1 / (x_exit/rho_v_fn(T_sat_chamber) + (1-x_exit)/rho_l_fn(T_sat_chamber));
    h_exit   = h_l_fn(T_sat_chamber) + x_exit*(h_v_fn(T_sat_chamber) - h_l_fn(T_sat_chamber));
    delta_h  = max(0, h_l_fn(T_tank) - h_exit);
    mdot_HEM = C_d * A_inj * rho_2ph * sqrt(2 * delta_h);
end

function [mdot_out, mdot_SPI, mdot_HEM] = compute_mdot(T_tank, p)
    rho_l  = p.rho_l_fn(T_tank);
    P_tank = p.Psat_fn(T_tank);
    DeltaP = max(P_tank - p.P_chamber, 0);
    mdot_SPI = p.C_d * p.A_inj * sqrt(2 * rho_l * DeltaP);
    s_tank = p.s_l_fn(T_tank);
    [mdot_HEM, ~, ~, ~] = hem_model(T_tank, s_tank, ...
        p.C_d, p.A_inj, p.T_sat_chamber, ...
        p.s_l_fn, p.s_v_fn, ...
        p.rho_l_fn, p.rho_v_fn, p.h_l_fn, p.h_v_fn);
    mdot_out = max(p.alpha * mdot_SPI + p.beta * mdot_HEM, 0);
end

function dxi = ode_rhs(~, xi, p)
    mu_l  = xi(1);
    Theta = xi(2);
    if mu_l <= p.mu_l_min
        dxi = [0; 0];
        return;
    end
    vrl  = p.varrho_l(Theta);
    vrv  = p.varrho_v(Theta);
    mu_v = vrv * (p.V_star - mu_l/vrl);
    T_tank  = Theta * p.T_star;
    [mdot_out, ~, ~] = compute_mdot(T_tank, p);
    mdot_nd = mdot_out * p.t_star / p.M_star;
    dvrl = p.dvarrho_l(Theta);
    dvrv = p.dvarrho_v(Theta);
    ahat_l = (mu_l/vrl^2) * dvrl;
    ahat_v = (mu_v/vrv^2) * dvrv;
    ghat   = 1/vrv - 1/vrl;
    Omega  = p.Omega_fn(Theta);
    Lambda = p.Lambda_fn(Theta);
    Psi    = p.Psi_fn(Theta);
    Pi     = p.Pi_fn(Theta);
    num       = -(Lambda/(vrl*ghat)) * mdot_nd;
    den       = mu_l + Omega*mu_v + (ahat_l+ahat_v)*(Lambda/ghat - Psi*Pi);
    Theta_dot = num / den;
    mdot_evap_nd = ((ahat_l+ahat_v)*Theta_dot + mdot_nd/vrl) / ghat;
    mu_l_dot     = -mdot_nd - mdot_evap_nd;
    dxi = [mu_l_dot; Theta_dot];
end

function mdot_evap = get_mdot_evap(mu_l, Theta, p)
    if mu_l <= p.mu_l_min
        mdot_evap = 0;
        return;
    end
    vrl  = p.varrho_l(Theta);
    vrv  = p.varrho_v(Theta);
    mu_v = vrv * (p.V_star - mu_l/vrl);
    T_tank  = Theta * p.T_star;
    [mdot_out, ~, ~] = compute_mdot(T_tank, p);
    mdot_nd = mdot_out * p.t_star / p.M_star;
    dvrl = p.dvarrho_l(Theta);
    dvrv = p.dvarrho_v(Theta);
    ahat_l = (mu_l/vrl^2) * dvrl;
    ahat_v = (mu_v/vrv^2) * dvrv;
    ghat   = 1/vrv - 1/vrl;
    Omega  = p.Omega_fn(Theta);
    Lambda = p.Lambda_fn(Theta);
    Psi    = p.Psi_fn(Theta);
    Pi     = p.Pi_fn(Theta);
    num       = -(Lambda/(vrl*ghat)) * mdot_nd;
    den       = mu_l + Omega*mu_v + (ahat_l+ahat_v)*(Lambda/ghat - Psi*Pi);
    Theta_dot = num / den;
    mdot_evap_nd = ((ahat_l+ahat_v)*Theta_dot + mdot_nd/vrl) / ghat;
    mdot_evap = mdot_evap_nd * p.M_star / p.t_star;
end

function [val, term, dir] = stop_event(tau, xi, p)
    val  = [xi(1) - p.mu_l_min; tau - p.tau_burn];
    term = [1; 1];
    dir  = [-1; -1];
end