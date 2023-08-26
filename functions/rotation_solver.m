
function A = rotation_solver(u_i_bar,DoR,p_bar,y_bar)

% Purpose - Apply rotation to velocity components
%
% Arguments IN:
% Mean velocities (u_i_bar)
% Degrees of rotation (DoR)
% mean pitch (p_bar)
% mean yaw (y_bar)
% Manually override solver and rely on calculated pitch and yaw
% 
% Arguments OUT:
% Rotation matrix (A)
%
% Version control:
% 14/05/2023 - Release
% 
% References:
% http://planning.cs.uiuc.edu/node102.html
% Code adapted from Mujal, et al. UPC (2015) https://github.com/vicente-medina/TurbulenceToolbox

% Initial guess at angular misalignment [pitch yaw]
x0 = [p_bar y_bar];

% Optimisation options (structure), where TolX = Termination tolerance on x, the current point, and TolFun = Termination tolerance on function value
options = optimset('TolX',1.1*10^(-5),'TolFun',1.1*10^(-7));

% Create a handle to control function parameters
f_handle = @(angles) rotation_matrix(angles,u_i_bar,DoR); % anonymous function with dummy variable 'angles'

% Call optimisation (solve system of nonlinear equations)
x = fsolve(f_handle,x0,options);

% Rotation matrix
if DoR == 2

    % Negative pitch value rotates the velocity vector towards the YZ plane
    x(1) = -p_bar;

    fprintf(['\n>> Pitch correction = %.3f rad (%.2f deg) \n>> Yaw correction = %.3f rad (%.2f deg) \n\n'],x(1),x(1)*180/pi,x(2),x(2)*180/pi)  
    A =     [cos(x(2))*cos(x(1)),   -sin(x(2)),      cos(x(2))*sin(x(1));
            sin(x(2))*cos(x(1)) ,    cos(x(2)) ,     sin(x(2))*sin(x(1));
            -sin(x(1))          ,    0         ,     cos(x(1))];

elseif DoR == 1
    fprintf('\n>> Yaw correction = %.3f rad (%.2f deg) \n\n',x(2),x(2)*180/pi)
    A = [cos(x(2))  ,   -sin(x(2))  ,   0;
        sin(x(2))   ,   cos(x(2))   ,   0;
        0           ,   0           ,   1];

elseif DoR == 0
    fprintf('\n>> Do not apply rotation in non-uniform flows \n\n')
    A = [1  ,   0  ,   0;
         0  ,   1  ,   0;
         0  ,   0  ,   1];
end
disp('------------------------------------------------------')

return
