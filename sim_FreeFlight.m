%% Configure Matlab Instance

clc; % Clear the command window
clear; % Clear workspace
close all; % Close all figures

% Setup plotting

set(0, 'DefaultFigureWindowStyle', 'normal');

%% Define constants

rEarth = 6378e3; % Earth radius @ equator [m]
muEarth = 3.986e14; % Earth gravitational parameter [m3s-2]

%% Define initial conditions

h0 = 500e3; % Initial orbit altitude [m]
t0 = datetime('now'); % Initial instant

%% Calculate initial state vector

r0 = [h0 + rEarth, 0, 0]; % Initial ECI position vector [m]
v0 = [0, sqrt(muEarth / norm(r0)), 0]; % Initial ECI velocity vector [m]
x0 = [r0, v0]; % Initial ECI state vector

%% Orbital period

torbit = 2 * pi * sqrt((norm(r0) ^ 3) / muEarth);

%% Propagate initial orbit

tprop = 0.95 * torbit; % Propagate the initial conditios for 95 % of an orbital period [s]

odeopts = odeset('AbsTol', 1e-9, 'RelTol', 1e-9); % Ode integration tolerances

[tres, xres] = ode45(@(t, x) freeFlightDynamics(t, x, muEarth), [0, tprop], x0, odeopts);

%% Show results

plotResults(tres, xres, t0);
