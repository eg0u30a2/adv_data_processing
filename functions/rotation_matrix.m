
function F = rotation_matrix(angles,u_i_bar,DoR)

% Purpose: compute the rotation matrix for v and z velocity components
% 
% Arguments IN:
% Mean velocities (u_i_bar)
% Degrees of rotation (DoR)
% Angles = variables to be solved
% 
% Arguments OUT:
% Rotation coefficients when applied to mean velocities (F)
%
% Version control:
% 30/06/2023 - Initial release
%
% NOTES:
% Code adapted from Mujal, et al. UPC (2015) https://github.com/vicente-medina/TurbulenceToolbox

if DoR == 2

    A = [
        sin(angles(2))*cos(angles(1)), cos(angles(2)), -sin(angles(2))*sin(angles(1)) ;
        sin(angles(1)),                 0,              cos(angles(1))];

elseif DoR == 1

    A = [
        sin(angles(2))  , cos(angles(2)), 0;
        0               , 0             , 1];
elseif DoR == 0
    A = [
        0, 1, 0;
        0, 0, 1]; 
end

% Apply rotation to mean velocity
F = A*u_i_bar';

end
