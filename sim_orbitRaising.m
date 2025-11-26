%% Configure Matlab Instance

clc; % Clear the command window
clear; % Clear workspace
close all; % Close all figures

% Setup plotting

set(0, 'DefaultFigureWindowStyle', 'normal');

% ODE integration tolerances

odeopts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9);

%% Define constants

rEarth = 6378e3; % Earth radius @ equator [m]
muEarth = 3.986e14; % Earth gravitational parameter [m3s-2]

%% Define initial conditions

h0 = 500e3; % Initial orbit altitude [m]
t0 = datetime('now'); % Initial instant

inc = 45; % Orbit inclination [º]
ecc = 0; % Orbit eccentricity [-]
RAAN = 0; % Orbit RAAN [º]
argp = 0; % Argument of the perigee [º] (Not relevant for circular orbits)
nu0 = 0; % Initial true anomaly [º]

%% Define target orbit altitude

% hf = 1000e3; % Target orbit altitude [m]
hf = 35786e3; % GSO altitude [m]

%% Calculate initial state vector

[r0, v0] = keplerian2ijk( ... 
    h0 + rEarth, ... 
    ecc, ... 
    inc, ... 
    RAAN, ... 
    argp, ... 
    nu0, ... 
    arglat = 0);

x0 = [r0, v0]; % Initial ECI state vector

%% Calculate transfer orbit velocity @ perigee

% The transfer orbit will have its perigee altitude equal to the initial
% orbit, and its apogee altitude equal to the target orbit

% To calculate the velocity at perigee of the transfer orbit, we use the
% Vis-Viva equation

aTransfer = ((h0 + rEarth) + (hf + rEarth)) / 2; % Semi major axis of the transfer orbit [m]

vTransferPer = sqrt(muEarth * ((2 / (h0 + rEarth)) - (1 / aTransfer))); % Velocity module at perigee of transfer orbit [m/s]

%% Calculate velocity module at target orbit

vF = sqrt(muEarth / (hf + rEarth)); % Velocity module at target orbit [m/s]

%% Propagate initial orbit

torbit = 2 * pi * sqrt((norm(r0) ^ 3) / muEarth); % Initial orbital period [s]

tprop = 0.5 * torbit; % Propagate the initial conditios for 50 % of an orbital period [s]

[tinitorb, xinitorb] = ode45( ... 
    @(t, x) freeFlightDynamics(t, x, muEarth), ... 
    [0, tprop], ... 
    x0, ... 
    odeopts);

%% Apply apogee-raising burn

% To obtain the velocity vector after the first, apogee raising burn, we
% calculate the unitary vector of the previous velocity to maintain the
% direction: v0t / norm(v0t)

% We then apply the desired velocity norm by multiplying element wise by
% the previously calculated velocity at transfer arc perigee

r0t = xinitorb(end, 1:3); % Position vector @ first burn [m]
v0t = xinitorb(end, 4:6); % Velocity vector before applying the first burn [m/s]

vTransferVec = (v0t ./ norm(v0t)) .* vTransferPer; % Velocity vector @ first burn [m/s]

x0t = [r0t, vTransferVec]; % Initial ECI state vector at transfer arc perigee

% Print Delta-V of the maneuver

dV1 = norm(vTransferVec) - norm(v0t); % Delta-V of the burn
fprintf('\nThe Delta-V of the apogee raising burn is: %5.2f\n', dV1); % Print to terminal

%% Propagate transfer arc

% We need to perform the second burn in the transfer apogee, since we begin
% in its perigee, so we need to propagate by exactly half of the orbital
% period of the transfer arc

tTransferOrbit = 2 * pi * sqrt((aTransfer ^ 3) / muEarth);
tTransfer = 0.5 * tTransferOrbit;

tspantransfer = [tinitorb(end), tinitorb(end) + tTransfer]; % Propagation timespan

[ttransferarc, xtransferarc] = ode45( ... 
    @(t, x) freeFlightDynamics(t, x, muEarth), ... 
    tspantransfer, ... 
    x0t, ... 
    odeopts);

%% Apply perigee raising burn

% We follow the same process as before, by keeping the last state vector
% and scaling the velocity by the value we calculated for the target orbit

r0f = xtransferarc(end, 1:3); % Position vector @ first burn [m]
v0f = xtransferarc(end, 4:6); % Velocity vector before applying the second burn [m/s]

vTargetVec = (v0f ./ norm(v0f)) .* vF; % Velocity vector @ second burn [m/s]

x0f = [r0f, vTargetVec]; % Initial ECI state vector at target orbit

% Print Delta-V of the maneuver

dV2 = norm(vTargetVec) - norm(v0f); % Delta-V of the burn
fprintf('\nThe Delta-V of the apogee raising burn is: %5.2f\n', dV2); % Print to terminal

%% Propagate target orbit

tTargetOrbit = 2 * pi * sqrt((norm(r0f) ^ 3) / muEarth); % Initial orbital period [s]

tTargetprop = 1 * tTargetOrbit; % Propagate the initial conditios for 50 % of an orbital period [s]

tspantarget = [ttransferarc(end), ttransferarc(end) + tTargetprop]; % Propagation timespan

[ttarget, xtarget] = ode45( ... 
    @(t, x) freeFlightDynamics(t, x, muEarth), ... 
    tspantarget, ... 
    x0f, ... 
    odeopts);

%% Join results

tres = [tinitorb; ttransferarc; ttarget];
xres = [xinitorb; xtransferarc; xtarget]; % Combine state vectors for plotting

%% Show results

plotResults(tres, xres, t0);
