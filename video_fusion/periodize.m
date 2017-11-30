function f = periodize(x, s, opt)
% PERIODIZE   Make the 2-D signal x periodize by s
%
%	f = periodize(x, s)
%
% Input:
%   x : 2-D signal
%   s : desired periodicity, typically 2^n
%   opt : control parameter
%           's'  thethe image is periodize by delay -1:1
%           'r'  image is periodize not by delay, running faster
%
% Output:
%   f  : periodic signal by s : f(p+s) = f(p)
%   
% Example 
%    x = rand(5)
%    f = periodize(x, 3)
%
% See also DELAY, CIRCSHIFT
%

if length(s) == 1
    s = [s s];
end

sz = size(x);
if ~exist('opt','var')
    in = sz./s;
    f = zeros(s);
    for in1 = 1:in(1)
        for in2 = 1:in(2)
            f = f+x((in1-1)*s(1)+1:in1*s(1), (in2-1)*s(2)+1:in2*s(2) );
        end
    end
    
else
    switch opt
        case {'s'}
            in = [1 1];
            f = 0*x;
            for in1 = -in(1):in(1)
                for in2 = -in(2):in(2)
                    f = f+delay(x,[in1, in2].*s);
                end
            end
        case {'r'}
            
            f = x(1:s(1),1:s(2));
            f(:,1:s(2)/2+1) = f(:,1:s(2)/2+1) + x(1:s(1),s(2)+1:end);
            f(1:s(1)/2+1,:) = f(1:s(1)/2+1,:) + x(s(1)+1:end,1:s(2));
            f(1:s(1)/2+1,1:s(2)/2+1) = f(1:s(1)/2+1,1:s(2)/2+1) + x(s(1)+1:end,s(2)+1:end);

            f = circshift(f,-[s(1)/4 s(2)/4]);
    end

end

% f = x+delay(x, [s(1) 0])+delay(x, [0 s(2)])+delay(x, [-s(1) 0])+ ...
%     delay(x, [0 -s(2)]) + delay(x, [s(1) s(2)])+delay(x, -[s(1) s(2)])+ ...
%     delay(x, [-s(1) s(2)])+delay(x, [s(1) -s(2)]);

