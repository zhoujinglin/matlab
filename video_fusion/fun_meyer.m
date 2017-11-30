function y = fun_meyer(x, param)
% FUN_MEYER   Return 1-D meyer function from x sample 
%
%	y = fun_meyer(x, param)
%
% Input:
%   x     : Sampling points of meyer function  
%   param : Parameter of the meyer function, length 4tion step
%
% Output:
%   y     : Return of meyer function sampling at x
%   
% Example 
%   plot(0:0.1:5, fun_meyer(0:0.1:5, 1:4))
% See also 
% size(x)
p = [-20 70 -84 35 0 0 0 0];

% y = polyval(p,x);

y = 0*x;
y( x <= param(1)) = 0;
y( (x >= param(1)) & (x <= param(2))) = polyval(p, ...
    (x((x >= param(1)) & (x <= param(2))) - param(1))./(param(2)-param(1)));
y( (x > param(2)) & (x <= param(3))) = 1;

y( (x >= param(3)) & (x <= param(4))) = polyval(p, ...
    (x((x >= param(3)) & (x <= param(4))) - param(4))./(param(3)-param(4)));
y( x > param(4)) = 0;
