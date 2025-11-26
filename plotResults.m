%
%           [] = plotResults(tres, xres)
%
% This function plots the results of the orbital maneuvers
%
% Inputs:
%   - tres (n x 1 vector): Time signal of the simulation [s]
%   - xres (n x 6 matrix): ECI states during the simulation [SI units]
%   - utc (1 x 6 vector): UTC time corresponding to the first element of
%   tres [MATLAB datetime object]
%

function [] = plotResults(tres, xres, utc)


%% Process inputs

npoints = length(tres); % Number of propagation points

if size(xres, 1) == size(xres, 2)
    warning('Please ensure that the xres matrix has the ECI state elements arranged by columns');
end

if size(xres, 2) == npoints
    xres = xres'; % Transform to major row if not the case
end

%% Transform ECI states into orbital elements, ECEF and spherical

% Preallocate transformations

orbel = zeros(npoints, 6);
xecef = zeros(npoints, 6);

% Orbital elements

[orbel(:, 1), ... 
    orbel(:, 2), ... 
    orbel(:, 3), ... 
    orbel(:, 4), ... 
    orbel(:, 5), ... 
    orbel(:, 6)] = ... 
    ijk2keplerian( ... 
    xres(:, 1:3)', ... 
    xres(:, 4:6)');

% ECEF

for ii = 1:npoints % The eci2ecef from MATLAB's Aerospace Toolbox is not vectorized

    [xecef(ii, 1:3), ... 
        xecef(ii, 4:6)] = ... 
        eci2ecef( ... 
        utc + seconds(tres(ii)), ... 
        xres(ii, 1:3), ... 
        xres(ii, 4:6));

end

% Spherical

[sph] = ... 
    ecef2lla( ... 
    xecef(:, 1:3));


%% Plot ECI 3D trajectory

figure('Name', '3D trajectory - ECI');

hold on;

% Draw the Earth as a sphere

[xsph, ysph, zsph] = sphere(50); % Create a sphere for the Earth (sphere outputs a radius 1 sphere)

xsph = xsph .* 6378e3; % Scale X coordinates
ysph = ysph .* 6378e3; % Scale Y coordinates
zsph = zsph .* 6378e3; % Scale Z coordinates

surf( ... 
    xsph, ... 
    ysph, ... 
    zsph, ... 
    'FaceColor', 'b', ... 
    'EdgeColor', 'none', ... 
    'HandleVisibility', 'off'); % Draw the Earth

% Plot trajectory

plot3( ... 
    xres(:, 1), ... 
    xres(:, 2), ... 
    xres(:, 3), ... 
    'LineWidth', 2, ... 
    'Color', 'red', ... 
    'DisplayName', 'Spacecraft trajectory');

% Draw Earth's rotation axis

plot3( ... 
    [0, 0], ... 
    [0, 0], ... 
    [-1.1 * 6378e3, 1.1 * 6378e3], ... 
    'LineWidth', 1.5, ... 
    'LineStyle', ':', ... 
    'Color', 'm', ... 
    'DisplayName', 'Rotation axis of the Earth');

% Configure plot

title('3D ECI trajectory');

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

legend;

grid off;

axis equal; % Set equal scaling

hold off;

% set(gca, 'Color', 'White'); % Set background color of the plot

%% Plot ECEF 3D trajectory

figure('Name', '3D trajectory - ECEF');

hold on;

% Draw the Earth as a sphere

[xsph, ysph, zsph] = sphere(50); % Create a sphere for the Earth (sphere outputs a radius 1 sphere)

xsph = xsph .* 6378e3; % Scale X coordinates
ysph = ysph .* 6378e3; % Scale Y coordinates
zsph = zsph .* 6378e3; % Scale Z coordinates

surf( ... 
    xsph, ... 
    ysph, ... 
    zsph, ... 
    'FaceColor', 'b', ... 
    'EdgeColor', 'none', ... 
    'HandleVisibility', 'off'); % Draw the Earth

% Plot trajectory

plot3( ... 
    xecef(:, 1), ... 
    xecef(:, 2), ... 
    xres(:, 3), ... 
    'LineWidth', 2, ... 
    'Color', 'red', ... 
    'DisplayName', 'Spacecraft trajectory');

% Draw Earth's rotation axis

plot3( ... 
    [0, 0], ... 
    [0, 0], ... 
    [-1.1 * 6378e3, 1.1 * 6378e3], ... 
    'LineWidth', 1.5, ... 
    'LineStyle', ':', ... 
    'Color', 'm', ... 
    'DisplayName', 'Rotation axis of the Earth');

% Configure plot

title('3D ECEF trajectory');

xlabel('X [m]');
ylabel('Y [m]');
zlabel('Z [m]');

legend;

grid off;

axis equal; % Set equal scaling

hold off;

% set(gca, 'Color', 'White'); % Set background color of the plot

%% Plot state vector evolution

figure('Name', 'State vector evolution');

% Plot altitude evolution

subplot(6, 2, [1, 3, 5]);

alt = vecnorm(xres(:, 1:3)') - 6378e3; % Altitude vector [m]
alt = alt' .* 1e-3; % Convert altitude to row vector and to km

plot( ... 
    tres, ... 
    alt, ... 
    'LineWidth', 2, ... 
    'Color', 'b');

title('Altitude vs time');

xlabel('Time [s]');
ylabel('Altitude [km]');

grid on;
grid minor;

hold off;

% set(gca, 'Color', 'White'); % Set background color of the plot

% Plot velocity evolution

subplot(6, 2, [7, 9, 11]);

vel = vecnorm(xres(:, 4:6)')'; % Altitude vector [m]

plot( ... 
    tres, ... 
    vel, ... 
    'LineWidth', 2, ... 
    'Color', 'b');

title('Inertial velocity vs time');

xlabel('Time [s]');
ylabel('Velocity [m/s]');

grid on;
grid minor;

hold off;

% set(gca, 'Color', 'White'); % Set background color of the plot

% Plot state vector elements evolution

titles = { ... 
    'X coordinate', ... 
    'Y coordinate', ... 
    'Z coordinate', ... 
    'X velocity', ... 
    'Y velocity', ... 
    'Z velocity'};

ylabels = { ... 
    'X [m]', ... 
    'Y [m]', ... 
    'Z [m]', ... 
    'X velocity [m/s]', ... 
    'Y velocity [m/s]', ... 
    'Z velocity [m/s]'};

for ii = 1:6

    subplot(6, 2, ii * 2);

    plot( ... 
        tres, ... 
        xres(:, ii), ... 
        'LineWidth', 2, ... 
        'Color', 'b');

    title(strcat(titles{ii}, ' vs time'));
    
    xlabel('Time [s]');
    ylabel(ylabels{ii});
    
    grid on;
    grid minor;
    
    hold off;
    
    % set(gca, 'Color', 'White'); % Set background color of the plot

end

% Final plot setting

% set(gcf, 'Color', 'White'); % Set background color of the plot

%% Plot orbit ground track

figure('Name', 'Trajectory ground track');

geoplot( ... 
    sph(:, 1), ... 
    sph(:, 2), ... 
    'LineWidth', 2, ... 
    'Color', 'g');

geobasemap 'streets-dark'

end
