%
%                 [dx] = freeFlightDynamics(t, x, mu)
%
% This function computes the state vector derivative for the non-linear
% dynamics of a spacecraft flying in a two-body system
%
% Inputs:
%   - t (scalar): Not used internally, needed for integration
%   - x (6 x 1 vector): ECI state vector [SI units]
%   - mu (scalar): Gravitational parameter of the central body [m3s-2]
%
% Output:
%   - dx (6 x 1 vector): Time derivative of the ECI state vector [SI units]
%

function [dx] = freeFlightDynamics(t, x, mu)

dx = zeros(size(x)); % Preallocate

dx(1:3) = x(4:6); % Instant derivative of the position vector

dx(4:6) = ((-1 * mu) / (norm(x(1:3)) ^ 3)) .* x(1:3); % Instant derivative of the velocity vector

end