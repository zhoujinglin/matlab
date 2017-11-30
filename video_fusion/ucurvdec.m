function y = ucurvdec(im, Cf, F)
% 
% UCURVDEC   Uniform Curvelet Decomposition
%
%	y = ucurvdec(im, Cf, F)
%
% Input:
%   Im : Input image
%   Cf : number of directional curvelet. In uniform curvelet, all elements 
%       are 3*2^n
%   F : Window cell matrices of same leng as F
%
% Output:
%   y : Cell structure contain coeff
%   
% Example 
%
% S = 512;
% N = [6 12];
% alpha = 0.15;
% r = pi*[0.3 0.5 0.85 1.15];q
% F = ucurv_win(S, N, r, alpha);
% y = ucurvdec(im, N, F);
% 
% See also UCURV_WIN, UCURV_REC
%

% the multiresolution length
Lt = length(Cf);
Sz = size(im);
imf = fft2(im);

% low resolution path
decim = 2^(Lt);
sz = Sz/decim;

imf2 = (1/prod(Sz./sz))*periodize(imf.*F{1}{1}, sz);
y{1} = decim*real(ifft2(imf2(1:sz(1),1:sz(2))));

for inres = 2:Lt+1
    % number of direction
    N = Cf(inres-1);
    
    % size of image at current resolution
    % S = Sz/2^(Lt-inres);
    
    % imf = fft2(im);

    % sz = S./2;
    % dm = [2 2];
    indec = Lt+1-inres;
    decim = 2^(indec);

    % iml2 = ifft2(imf.*fl);
    % imf2 = (1/prod(dm))*periodize(imf.*fl, sz);
    % lowpass multiply by 2 
    % y{1} = 2*real(ifft2(imf2(1:sz(1),1:sz(2))));

    for in = 1:N
        if in > N/2
            dm = decim*[2 N/3];
        else
            dm = decim*[N/3 2];
        end
        
        sz = Sz./dm;

% disp('window in')
% F{inres}{in}
        imf2 = imf.*F{inres}{in};
        
        % the nominator is the scaling of complex band
        % the denominator compensate for periodizing
%	disp('periodization ')
	periodize(imf2, sz);
        imf2 = ((sqrt(2*prod(dm)))/prod(dm))*periodize(imf2, sz);

% disp('one band');
y{inres}{in} = ifft2(imf2(1:sz(1),1:sz(2)));
y{inres}{in};

    end
    
%     im = y{1};
end
