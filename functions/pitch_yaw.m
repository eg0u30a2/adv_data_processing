
% Purpose - calculate pitch and yaw. Roll is rotation of the vector and does not affect the position in 3D space.

function [y_bar,p_bar] = pitch_yaw(u_i_pre_rot,plot)

% Arguments OUT:
% y_bar = mean yaw angle for the dataset
% p_bar = mean pitch angle for the dataset 
%
% Arguments IN:
% u_i_pre_rot = table of velocity components (m/s)
% plot - flag to indicate whether histograms of pitch and yaw are to be plotted
%
% VERSION CONTROL:
% 14/05/2023 - Updated release
%
% NOTES:
% https://stackoverflow.com/questions/58469297/how-do-i-calculate-the-yaw-pitch-and-roll-of-a-point-in-3d

tf = istable(u_i_pre_rot);

if tf ~= 1
    u_i_pre_rot = array2table(u_i_pre_rot(:,1:3));
    u_i_pre_rot.Properties.VariableNames = {'Vx' 'Vy' 'Vz'};
else
end

if nargin == 1
    u_i_bar = mean(table2array(u_i_pre_rot));
    plot = 0;
else
end

% Yaw calculation
y = atan(u_i_pre_rot.Vy ./ u_i_pre_rot.Vx);
y_bar = mean(y,'omitnan');

% Pitch calculation
p = atan(u_i_pre_rot.Vz ./ (((u_i_pre_rot.Vx.^2)+(u_i_pre_rot.Vy.^2)).^0.5));
p_bar = mean(p);

% Plot
if plot
    % Add plot.. 
else
end

end
