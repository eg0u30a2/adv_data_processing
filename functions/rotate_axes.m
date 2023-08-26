
function [u_i_post_rot] = rotate_axes(u_i_pre_rot,fs,t,DoR)

% Purpose: rotate the axes to obtain zero Y velocity. Z velocity rotated according to prescribed correction.
%
% Arguments IN:
% DoR = number of degrees of rotation (see below for description)
% fs = sampling freq. (Hz)
% t = time vector (s)
% u_i_pre_rot = input 2D velocity array (pre-rotated) (m/s)
% 
% Arguments OUT:
% u_i_post_rot = rotated velocity 2D array acc. to the number of DoR (m/s)
%
% VERSION CONTROL:
% 17/03/2023 - Removal of static pitch offset for rotational correction
%
% NOTES:
% Use the yaw (alpha), pitch (beta), roll (gamma) rotation angles.
% http://planning.cs.uiuc.edu/node102.html
%   
% 2 Degrees of Rotation (DoR)
% Assume that the roll (gamma) is 0.
% Rotation matrix for yaw Y & pitch Z (sin(alpha)=sa, cos(alpha)=ca, sin(beta)=sinb, cos(beta)=cb)
%
%       [ca·cb	-sa	 ca·sb]
%       [sa·cb	 ca	 sa·sb]
%       [-sb	 0	 cb   ]
%
% Rotate the axes until the vertical and transversal velocities are 0
% 
% 1 Degree of Rotation
% Assumes roll is 0 and no pitch correction is to be applied, thus b = 0.
% Rotation matrix for yaw Y (sin(alpha)=sa, cos(alpha)=ca)
%
%       [ca	-sa	 0]
%       [sa	 ca	 0]
%       [0	 0	 1]
%
% Applies only when the data corresponds to a straight channel (e.g. a flume) with negligible secondary currents


% Calculate mean values
u_i_bar = mean(u_i_pre_rot);

% Return average pitch and yaw angles
[y_bar,p_bar] = pitch_yaw(u_i_pre_rot);

if DoR == 0
    u_i_post_rot = u_i_pre_rot;
    return
else
end

% Optimise yaw and pitch to define the rotation matrix A
A = rotation_solver(u_i_bar,DoR,-p_bar,y_bar);
       
% Apply the rotation
u_i_post_rot = u_i_pre_rot * (A');

end