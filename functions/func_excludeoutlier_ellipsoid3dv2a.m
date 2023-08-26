function [iOutlier, coef] = func_excludeoutlier_ellipsoid3dv2a(xi,yi,zi,theta)
%======================================================================
%
% Version 1.01
%
% This program excludes the points outside of ellipsoid in two-
% dimensional domain
%
% Input
%   xi : input x data
%   yi : input y data
%   zi : input z data
%   theta  : angle between xi and zi
%
% Output
%   iOutlier : logical mask of excluded xi and yi
%   coef : coefficients for ellipsoid
%
% Example:
%   [xp,yp,zp,ip,coef] = func_excludeoutlier_ellipsoid3d(f,f_t,f_tt,theta);
%
%
%======================================================================
% Terms:
%
%       Distributed under the terms of the terms of the BSD License
%
% Copyright:
%
%       Nobuhito Mori, Kyoto University
%
%========================================================================
%
% Update:
% 1.xx 2019/05/08 Simon Spagnol : code simplifications
% 1.01    2009/06/09 Nobuhito Mori
% 1.00    2005/01/12 Nobuhito Mori
%
%========================================================================

%% Initial setup

n = length(~isnan(xi));
lambda = sqrt(2*log(n));

%% Rotate data

%theta = atan2( nansum(xi.*zi), nansum(xi.^2) );

if abs(theta) < 1e-6
    X = xi;
    Y = yi;
    Z = zi;
else % apply rotation matrix yielding X, Y and Z projections in phase space
    R = [ cos(theta) 0  sin(theta); 0 1 0 ; -sin(theta) 0 cos(theta)]; % Pitch = counterclockwise rotation of theta about y axis
    X = xi*R(1,1) + yi*R(1,2) + zi*R(1,3);
    Y = xi*R(2,1) + yi*R(2,2) + zi*R(2,3);
    Z = xi*R(3,1) + yi*R(3,2) + zi*R(3,3);
end

%% Pre-process

% Expected maxima using Univeral Criterion
a = lambda*nanstd(X);
b = lambda*nanstd(Y);
c = lambda*nanstd(Z);

%% Main

coef = zeros(3,1);
x2 = zeros(size(X));
y2 = zeros(size(Y));
zt = zeros(size(Z));
z2 = zeros(size(Z));

% point on the ellipsoid
x2 = a*b*c*X ./ sqrt((a*c*Y).^2 + b^2*(c^2*X.^2 + a^2*Z.^2));
y2 = a*b*c*Y ./ sqrt((a*c*Y).^2 + b^2*(c^2*X.^2 + a^2*Z.^2));
zt = c^2* ( 1 - (x2./a).^2 - (y2./b).^2 );
z2 = sign(Z) .* sqrt(zt);

% check outlier from ellipsoid
dis = (x2.^2 + y2.^2 + z2.^2) - (X.^2 + Y.^2 + Z.^2);

iOutlier = dis < 0; ; % index the time series of the outlier

% Ellipsoid coefficients, i.e. expected maxima 
coef(1) = a;
coef(2) = b;
coef(3) = c;

end