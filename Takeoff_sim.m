function [s_to, t_to, v_max, v_to] = simulate_aircraft(c, show_plots)

if nargin < 1
        clc; close all;
        c = 0.25; % Change chord here
end
if nargin < 2
    show_plots = true;
end

% ======== Assumptions ========

% ====== wing
b = 1.1; % Wingspan m
rho_wing = 38.15; %density of foam
avg_thick = 0.08; % as % of chord
max_thick = 0.117;
Cl_max = 1.4; %max lift coef
Cl_roll = 0.6; % assumes rig angle of 2deg

% ====== mass
mass_non_wing = 1.33; % empty mass w/o wing or golf balls
mass_ball = 0.046;
n_ball = 20; % number of balls carried

% ====== Body
Vol = 0.043; % Fuselage Volume
A_w_body = 0.68; % Wetted area w/o wing (m^2)
W_fus = 0.4; %width of fuselage (m)

% ===== Environment
mu = 0.0000182; %Dynamic viscocity
fric = 0.35; %grass friction coeff
rho = 1.225; % Air density (kg/m^3)
g = 9.81; % Gravitational acceleration (m/s^2)

% ========== Chord dependant ==========

S = b*c; % Planfrom aera m^2
AR = b/c;
mass_wing = rho_wing*S*c*avg_thick*1.2; %1.2 to account for extra spar
% Scale 2D airfoil Cl values for 3D wingtip losses:

m = mass_non_wing + n_ball * mass_ball + mass_wing; %total mass

A_f = max_thick*S; %wing front profile area
A_w = A_w_body + 2*sqrt(1 + max_thick^2)*c*(b - W_fus) + W_fus*c; %Wing wetted area
e = 1.78*(1-0.045*AR^0.68)-0.64; %osawald eff factor approximation
K = 1/(pi*AR*e);

v_stall = sqrt((m*g)/(0.5*rho*S*Cl_max));
v_to = 1.2*v_stall;

Components = [ % [Profile Area (m^2), Cd]
    A_f 0.062; % Wing
    0.002 1.4; % Undercarriage
    0.0004*3 0.6; % Wheel*3
    Vol^(2/3) (0.04*Vol^(2/3))/S; % Fuselage
    0.0017 1.4]; % Tail Profile
Cd0_i = sum(prod(Components, 2)) / S; %parasitic drag

% ========== sim setup
time = 20; % (simulation time (s))
dt = 0.1;
v = 0.01;
t = 0;
s = 0;

s_to = 0; t_to = 0;

takeoff = false;

n = time/dt;
output = zeros(n, 4);

% =============== Simulation Loop ===============
for i = 1:n

    % ========== chord dependant

    Re = (rho*v*c)/mu;
    C_f = 0.455*(log10(Re))^-2.58; %Cd for skin friction, from myplace doc
    Cd0 = Cd0_i + A_w*C_f/S; % add skin friction

    T = -0.003036 *v^2 - 0.4814 *v + 18.8;


    if ~takeoff
        if v >= v_to
            takeoff = true;
            s_to = s; t_to = t; Re_to = Re;
            Cl = (m * g) / (0.5 * rho * v^2 * S); % find c_l required
        else
            Cl = Cl_roll;
        end
    else
        % Find Cl required for level flight
        Cl = (m * g) / (0.5 * rho * v^2 * S);
    end

    L = 0.5 * rho * v^2 * S * Cl;
    N = (m*g - L)*(~takeoff); %
    D = 0.5 * rho * v^2 * S * (Cd0 + K * Cl^2);
    
    F_x = T - D - fric*N;
    a_x = F_x/m;

    s = s + v*dt + 0.5*a_x*dt^2;
    v = v + a_x*dt;
    t = t + dt;

    output(i, :) = [t v s Cl];
end

v_max = max(output(:, 2));

% =============== Plotting & Console Summary ===============
if show_plots
    fprintf('================ SIMULATION RESULTS ================\n');
    fprintf('Chord (c):          %.3f m  (AR = %.2f, S = %.3f m^2)\n', c, AR, S);
    fprintf('Total Mass (m):     %.3f kg (Wing: %.3f kg)\n', m, mass_wing);
    fprintf('Stall Speed:        %.2f m/s\n', v_stall);
    fprintf('Takeoff Speed:      %.2f m/s\n',  v_to);
    fprintf('Takeoff Time:       %.2f s\n', t_to);
    fprintf('Takeoff Distance:   %.2f m\n', s_to);
    fprintf('Max Speed (v_max):  %.2f m/s\n [Note: this is unrealistic]', v_max);
    fprintf('Takeoff Reynolds    %.2f\n', Re_to);  
    fprintf('Max Reynolds:       %.2f\n', Re);
    fprintf('====================================================\n');

    figure('Name', sprintf('Aircraft Simulation (c = %.2f m)', c), 'Position', [100 100 1000 700]);
    tiledlayout(2, 2);

    % --- 1. Velocity vs Time ---
    nexttile
    plot(output(:,1), output(:,2), 'LineWidth', 1.5); hold on;
    if ~isnan(t_to)
        xline(t_to, '--r', sprintf('Takeoff (%.2f s)', t_to));
        yline( v_to, ':w', sprintf('V_{to} = %.1f m/s',  v_to));
    end
    yline(v_max, '-.g', sprintf('V_{max} = %.1f m/s', v_max));
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    title('Velocity vs. Time');
    grid on;

    % --- 2. Velocity vs Distance ---
    nexttile
    plot(output(:,3), output(:,2), 'LineWidth', 1.5); hold on;
    if ~isnan(s_to)
        xline(s_to, '--r', sprintf('s_{to} = %.2f m', s_to));
        yline( v_to, ':w', 'V_{takeoff}');
        xlim([0, max(20, s_to * 3)]); % Zoom to show takeoff roll + initial climb clearly
    end
    yline(v_max, '-.g', 'V_{max}');
    xlabel('Distance (m)');
    ylabel('Velocity (m/s)');
    title('Velocity vs. Distance');
    grid on;

    % --- 3. Distance vs Time ---
    nexttile
    plot(output(:,1), output(:,3), 'LineWidth', 1.5); hold on;
    if ~isnan(t_to)
        xline(t_to, '--r', 'Takeoff');
        yline(s_to, ':w', sprintf('%.2f m', s_to));
        xlim([0, max(5, t_to * 2.5)]);
        ylim([0, max(20, s_to * 3)]);
    end
    xlabel('Time (s)');
    ylabel('Distance (m)');
    title('Distance vs. Time');
    grid on;

    % --- 4. Acceleration vs Time (Bonus diagnostic tile) ---
    nexttile
    plot(output(:,2), output(:,4), 'LineWidth', 1.5); hold on;
    if ~isnan(v_to)
        xline(v_to, '--r', 'Takeoff');
    end
    yline(0, ':w');
    xlabel('Velocity (m/s)');
    ylabel('Cl');
    title('Cl vs. Vel');
    grid on;
end

end
