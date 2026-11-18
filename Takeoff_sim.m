% ------------------- Variables -------------------------

% ------ Plane ------
m = 2.324; %mass kg
c = 0.25; %chord m
b = 1.1; % Wingspan m
S = b*c; %Planfrom wing aera m^2
Cl = 1.35; %lift coeff
Re = 316181; %Reynolds Number
Vol = 0.0175; % Fuselage Volume

Cf = 0.455*(log10(Re))^-2.58; %Cd for skin friction, from myplace doc

Components = [ % [Profile Area (m^2), Cd]
    0.00915 0.062; % Wing
    0.002 1.4; % Undercarriage
    0.0004 0.6; % Wheel
    0.02475 Cf; % Skin Friction (use wetted area)
    Vol (0.0273*Vol^0.666666)/0.275; % Fuselage
    0.0017 1.4]; % Tail Profile

% ------ Environmental ------
mu = 0.35; %grass friction
g = 9.81;
rho = 1.225; % kg/m^3

% ------ Derived ------
AR = b/c; %aspect ratio
Cd0 = sum(prod(Components, 2)) / S; %parasitic drag

e = 1.78*(1-0.045*AR^0.68)-0.64; %osawald eff factor approximation
K = 1/(pi*AR*e);

C_induced = K*Cl^2; % Induced Drag

V_stall = sqrt((m*g)/(0.5*rho*S*Cl));
V_takeoff = 1.2*V_stall;

% ------ Thrust & Velocity ------
dt = 0.05; % time step

T = @(V) -0.003036 * V.^2 - 0.4814 * V + 18.8; %A, B, C from prop crve

L = @(V) 0.5 * rho * V.^2 * S * Cl;
D = @(V) 0.5 * rho * V.^2 * S * (Cd0 + C_induced);
N = @(V) max(0, m*g - L(V)); % set normal force to 0 after lift off

delV = @(V) ( (T(V) - mu*N(V) - D(V)) / m ) * dt;

% ------------------- Simulation -------------------------

%initial conditions

V = 0; 
s = 0; %distance
t = 0;

n = ceil(10/dt); %number of iteraitons

output = zeros(n, 3);


for i = 2:n
    dV = delV(V);
    V = V + dV;
    s = s + V*dt;
    t = t + dt;
    output(i, :) = [t V s];
end

figure
tiledlayout(2,2)

idx_to = find(output(:,2) >= V_takeoff, 1, 'first');
t_to   = output(idx_to, 1); % Time at takeoff [s]
s_to   = output(idx_to, 3); % Distance at takeoff [m]

% --- Velocity vs Time ---
nexttile
plot(output(:,1), output(:,2)); hold on
xline(t_to, '--r', 'Takeoff');
yline(V_takeoff, ':k', 'V_{takeoff}');           
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity vs. Time');
grid on;

% --- Velocity vs Distance ---
nexttile
plot(output(:,3), output(:,2), '-'); hold on
xline(s_to, '--r', 'Takeoff');
yline(V_takeoff, ':k', 'V_{takeoff}');
xlabel('Distance (m)');
ylabel('Velocity (m/s)');
title('Velocity vs. Distance');
grid on;

% --- Distance vs Time ---
nexttile
plot(output(:,1), output(:,3), '-'); hold on
xline(t_to, '--r', 'Takeoff');
yline(s_to, ':k', sprintf('%.1f m', s_to));      
xlabel('Time (s)');
ylabel('Distance (m)');
title('Distance vs. Time');
grid on;
