function [fo, ip] = func_despike_phasespace3d_wahl( fi, i_plot, i_opt)
%======================================================================
%
% Version 1.2
%
% This subroutine excludes spike noise from Acoustic Doppler
% Velocimetry (ADV) data using phase-space method, using
% modified Goring and Nikora (2002) method by Nobuhito Mori (2007).
% Further modified by Joseph Ulanowski to remove offset in output (2014).
%
% Further modified by Simon Spagnol
%   - nansum for atan2 calculation in case pass data with NaN
%   - control printing
%   - remove offset in output (2018).
%
%======================================================================
%
% Input
%   fi     : input data with dimension (n,1)
%   i_plot : = 9 plot results (optional)
%   i_opt  : = 0 or not specified ; return spike noise as NaN
%            = 1            ; remove spike noise and variable becomes shorter than input length
%            = 2            ; interpolate NaN using cubic polynomial
%
% Output
%   fo     : output (filtered) data
%   ip     : excluded array element number in fi
%
% Example:
%   [fo, ip] = func_despike_phasespace3d( fi, 9 );
%     or
%   [fo, ip] = func_despike_phasespace3d( fi, 9, 2 );
%
%
%======================================================================
% Terms:
%
%       Distributed under the terms of the terms of the BSD License
%
% Copyright:
%
%       Nobuhito Mori
%           Disaster Prevention Research Institue
%           Kyoto University
%           mori@oceanwave.jp
%
%========================================================================
%
% Update:
% x.x : 2019/05/08 : nansum for atan2 calculation, control printing
%       1.2     2014/03/18 Offset removed for non-zero mean data [J.U.]
%       1.11    2009/06/09 Interpolation has been added.
%       1.01    2009/06/09 Minor bug fixed
%       1.00    2005/01/12 Nobuhito Mori
%
%========================================================================

nvar = nargin;
if nvar == 1
    i_opt  = 0;
    i_plot = 0;
elseif nvar == 2
    i_opt  = 0;
end

%% Initial setup

% Maximum number of iterations
n_iter = 20;
n_out  = 999;

n      = length(fi);
f_median = 0.0;

% Scale parameter of threshold - Median of Absolute Deviations (MAD) rather than std deviation 
isMAD = false(size(fi)); % false returns an array of logical zeros
lambda = sqrt(2*log(n)); % Universal threshold
eps = 1e-6;
f_median = nanmedian(fi);
isMAD = abs(fi-f_median) > (1.483*lambda*nanmedian(abs(fi-f_median)) + eps); % returns logical 1 where abs diff exceeds MAD
f      = fi;
f(isMAD) = NaN;

f_median = 0.0; % reset median value
n = length(~isnan(f));
lambda = sqrt(2*log(n));

%% Loop

% Initialise iterations required to update mean values
% At each iteration, identification of the spikes reduces the standard deviations and thus the size of the ellipsoid reduces
n_loop = 1;

while (n_out ~= 0) & (n_loop <= n_iter)
    % step 0
    isSpike = false(size(f));
    f_median = f_median + nanmedian(f); % accumulate offset value at each step [J.U.]
    f = f - nanmedian(f); % apply offset (remove median) to input data
    % fprintf('Standard deviation = %.6e m/s \n',nanstd(f))
    
    % step 1: 1st and 2nd partial derivatives (note: no division by time step)
    f_t  = gradient(f);
    f_tt = gradient(f_t);
    
    % step 2: estimate angle between f and f_tt axis, i.e. rotation angle of principal axis using cross-correlation
    % four-quadrant inverse tangent, atan2(Y,X), returns values (radians) in the closed interval [-pi,pi] based on the values of Y and X.
    if n_loop==1
        theta = atan2( nansum(f.*f_tt), nansum(f.^2) );
    end
    
    % step 3: checking outlier in the 3D phase space
    [isSpike, coef] = func_excludeoutlier_ellipsoid3dv2a(f,f_t,f_tt,theta);

    % Exclude data
    n_nan_1 = size(find(isnan(f)==1),1);
    f(isSpike) = NaN; % replace indices of ellipsoid outliers with NaN in the velocity data 
    n_nan_2 = size(find(isnan(f)==1),1);
    n_out = n_nan_2 - n_nan_1; % converges to zero once (provided max n_iter not reached)
    n_outliers = sum(isnan(f)); % converges once no new NaNs are identified
   
    % --- end of loop
    n_loop = n_loop + 1;
end

%% Post-process

go = f + f_median;    % reapply offset of array mean
ip = find(isnan(go)); % returns indices of outliers

if n_loop < n_iter
%     disp(['>> Number of outliers = ', num2str(size(find(isnan(f)==1),1))])
%     disp(['>> Number of iterations = ', num2str(n_loop-1)])
else
%     disp(['>> Number of outliers = ', num2str(size(find(isnan(f)==1),1))])
%     disp(['>> Number of iterations = ', num2str(n_loop-1)])
%     disp(['>> Error: iteration limit reached!'])
end

%% Interpolation or shorten NaN data

if abs(i_opt) >= 1
    % remove NaN from data
    isSpike = find(~isnan(go));
    fo = go(isSpike);
    % interpolate NaN data
    if abs(i_opt) == 2
        x   = find(~isnan(go));
        y   = go(x);
        xi  = 1:max(length(fi));
        fo = interp1(x, y, xi, 'cubic')';
    end
else
    % output despiked value as NaN
    fo = go;
end

%% Check and plot

if i_plot == 9
    
    %theta/pi*180
    F    = fi - f_median; % remove offset (mean) from unfiltered input array
    F_t  = gradient(F); % 1st derivative
    F_tt = gradient(F_t); % 2nd derivative
    
    % Rotation matrices (opposing directions)
    RF = [ cos(theta) 0  sin(theta); 0 1 0 ; -sin(theta) 0 cos(theta)];
    RB = [ cos(theta) 0 -sin(theta); 0 1 0 ;  sin(theta) 0 cos(theta)];
    
    % making ellipsoid data
    a = coef(1);
    b = coef(2);
    c = coef(3);
    ne  = 32; % no. elements
    dt  = 2*pi/ne;
    dp  = pi/ne; % angular increment
    t   = 0:dt:2*pi;
    p   = 0:dp:pi;
    n_t = max(size(t));
    n_p = max(size(p));
    
    % making ellipsoid
    for it = 1:n_t
        for is = 1:n_p
            xe(n_p*(it-1)+is) = a*sin(p(is))*cos(t(it));
            ye(n_p*(it-1)+is) = b*sin(p(is))*sin(t(it));
            ze(n_p*(it-1)+is) = c*cos(p(is));
        end
    end
    xer = xe*RB(1,1) + ye*RB(1,2) + ze*RB(1,3);
    yer = xe*RB(2,1) + ye*RB(2,2) + ze*RB(2,3);
    zer = xe*RB(3,1) + ye*RB(3,2) + ze*RB(3,3);
    
end
