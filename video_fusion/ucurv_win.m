function [F, ind, cf] = ucurv_win(Sz, Cf, r, alpha, fo)
% UCURV_WIN   Generate the curvelet windows that used in uniform curvelet
% toolbox
%
%	F = ucurv_win(Sz, Cf, r, alpha, [fo])
%
% Input:
%   Sz : size of the generated window
%   Cf : number of directional curvelet. In uniform curvelet, n is 3*2^n
%   r : parameter of the meyer window used in diameter direction
%   alpha : paramter for meyer window used in angle function
%   fo : [optional] format option, if yes then F is in compressed format.
%
% Output:
%   F : cell of size n containing 2-D matrices of size s*s
%   in case of compressed format
%   F : a column of curvelet window value that are diff. from zero
%   ind : index of position correspond to value stored in F
%   cf  : configuration paramter, specify ending point of window value in F
%   and ind.
%   
% Example 
%
% S = 512;
% N = [6 12];
% alpha = 0.15;
% r = pi*[0.3 0.5 0.85 1.15];
% F = ucurv_win(S, N, r, alpha);
%
% See also FUN_MEYER
%

if ~exist('fo','var')
    fo = 0;
end

% the multiresolution length
Lt = length(Cf);

r = [r(1:2); r(3:4)];
for in = 1: Lt-1
    r = [0.5*r(1,:); r ];
end

% Size of window
if max(size(Sz)) == 1
    Sz = [Sz Sz];
end

cs = Sz;
% create the grid
S1 = -1.5*pi:pi/(cs(1)/2):0.5*pi;
S2 = -1.5*pi:pi/(cs(2)/2):0.5*pi;
[x1, x2] = meshgrid(S2,S1);

% scale the grid approximate the tan theta function ------------------
% creat two scale grid for mostly horizontal and vertical direction
% firt grid
t1 = zeros(size(x1));
ind = and(x1~=0, abs(x2) <= abs(x1));
t1(ind) = -x2(ind)./x1(ind);

t2 = zeros(size(x1));
ind = and(x2~=0, abs(x1) < abs(x2));
t2(ind) = x1(ind)./x2(ind);
t3 = t2;
t3(t2<0) = t2(t2<0)+2;
t3(t2>0) = t2(t2>0)-2;

M1 = t1+t3;
M1(x1>=0) = -2;

% second grid
t1 = zeros(size(x1));
ind = and(x2~=0, abs(x1) <= abs(x2));
t1(ind) = -x1(ind)./x2(ind);

t2 = zeros(size(x1));
ind = and(x1~=0, abs(x2) < abs(x1));
t2(ind) = x2(ind)./x1(ind);
t3 = t2;
t3(t2<0) = t2(t2<0)+2;
t3(t2>0) = t2(t2>0)-2;

M2 = t1+t3;
M2(x2>=0) = -2;
clear t1 t2 t3;

% for multiresolution window
rd = max(abs(x1),abs(x2));

% for each resolution estimation
for inres = 2:Lt+1
    
    
    % parameter for the resolution
    % current size of image at resoltuion
    % cs = Sz./2^(Lt-inres);
    n = Cf(inres-1);
    
    if or((mod(n,6) ~=0), (2^nextpow2(n/6)  ~= n/6))
        disp('Error: number of direction must be in the form 3*2^n');
        F = 0;
        return;
    end
        
    % multiresolution window ---------------------------------------------
    % Low pass window
    flrow = fun_meyer(abs(S1),[-2 -1 r(inres-1,:)]);
    flcol = fun_meyer(abs(S2),[-2 -1 r(inres-1,:)]);
    FL =  flrow'*flcol;
    fhrow = fun_meyer(abs(S1),[-2 -1 r(inres,:)]);
    fhcol = fun_meyer(abs(S2),[-2 -1 r(inres,:)]);
    FH =  fhrow'*fhcol;
    % high pass window
    fr = FH-FL;
    
    % lowpass band
    if inres == 2
        if fo
            sbind = 1;
            tmp = sqrt(circshift(FL(1:cs(1),1:cs(2)), 0.25*cs ));
            ind = find(tmp);
            F = tmp(ind);
            cf = length(ind);
        else
        F{1}{1} = sqrt(circshift(FL(1:cs(1),1:cs(2)), 0.25*cs ));
        end 
    end
    
    
    % angle meyer window
    angd = 4/n;
    ang = angd*[-alpha alpha 1-alpha 1+alpha];

    if (n == 6)
        in = 1;
        
        % band 1 and 3
        ang2 = -1+(in-1)*angd+ang;
        fang =  fun_meyer(M1,ang2);
        f = fang.*fr;
        % f = fliplr(f);
        if fo
            tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, ...
                circshift(flipud(tmp), [1 0]) );
        else
	    F{inres}{n/2+1-in} = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            F{inres}{in} = circshift(flipud(F{inres}{n/2+1-in}), [1 0]);
        end

        % band 4 and 6
        fang =  fun_meyer(M2,ang2);
        f = fang.*fr;
        % f = flipud(f);
        
        if fo
            tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, ...
                circshift(fliplr(tmp), [0 1]) );
        else
            F{inres}{in+n/2} = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            F{inres}{n+1-in} = circshift(fliplr(F{inres}{in+n/2}), [0 1]);
        end

        in = 2;
        % band 2 and 5
        ang2 = -1+(in-1)*angd+ang;
        fang = fun_meyer(M1,ang2);
        f = fang.*fr;
        % f = fliplr(f);
        
        if fo
            tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
        else
            F{inres}{in} = sqrt(circshift(f(1:cs(1),1:cs(2)),  [0.25*cs(1), 0.25*cs(2)] ));
        end
        
        fang = fun_meyer(M2,ang2);
        f = fang.*fr;
        % f = flipud(f);
        if fo
            tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
            [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
        else
            F{inres}{n+1-in} = sqrt(circshift(f(1:cs(1),1:cs(2)),  [0.25*cs(1), 0.25*cs(2)] ));
        end


    else
        % estimate the window for the first half of all filters
        for in = 1:n/4
            ang2 = -1+(in-1)*angd+ang;
            fang =  fun_meyer(M1,ang2);

            f = fang.*fr;
            % f = fliplr(f);
            % normalization and fit the window
            % then shift and square root
            if fo
                tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
                [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
                [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, ...
                    circshift(flipud(tmp), [1 0]) );
            else
                F{inres}{n/2+1-in} = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
                F{inres}{in} = circshift(flipud(F{inres}{n/2+1-in}), [1 0]);
            end

        end
        
        % estimate the window for the second half of all filters
        for in = 1:n/4
            ang2 = -1+(in-1)*angd+ang;
            fang =  fun_meyer(M2,ang2);

            f = fang.*fr;
            % f = flipud(f);
            % normalization and fit the window
            % then shift and square root
            if fo
                tmp = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
                [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp);
                [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, ...
                    circshift(fliplr(tmp), [0 1]) );
            else
                F{inres}{in+n/2} = sqrt(circshift(f(1:cs(1),1:cs(2)), [0.25*cs(1), 0.25*cs(2)] ));
                F{inres}{n+1-in} = circshift(fliplr(F{inres}{in+n/2}), [0 1]);
            end

        end
    end
    
end

if ~fo
    cf = [];
    ind = [];
end

end


% internal function
    function [sbind, F, ind, cf] = addwindow(sbind, F, ind, cf, tmp)
        sbind = sbind+1;
        ind2 = find(abs(tmp)>10^(-9));
        ind = [ind;ind2];
        
        F = [F;tmp(ind2)];
        cf = [cf;length(F)];
        
        %         tmp2 = 0*tmp;
        %         tmp2(ind2) = 1;
        %         imagesc(tmp2);
        %         pause
    end


