function im2 = ucurvrec(y, Cf, F)
% 
% UCURVREC   Uniform Curvelet Reconstruction
%
%	im2 = ucurvrec(y, Cf, F)
%
% Input:
%   y : Cell structure contain coeff
%   Cf : number of directional curvelet. In uniform curvelet, all elements 
%       are 3*2^n
%   F : Window cell matrices of same leng as F
%
% Output:
%   Im : Reconstructed image
%   
% Example 
%
% S = 512;
% N = [6 12];
% alpha = 0.15;
% r = pi*[0.3 0.5 0.85 1.15];
% F = ucurv_win(S, N, r, alpha);
% y = ucurvdec(im, N, F);
% im2 = ucurvrec(y, N, F);
% 
% See also UCURV_WIN, UCURV_DEC
%

% the multiresolution length
Lt = length(Cf);
% low resolution image
im2 = y{1};

decim = 2^Lt;
sz = size(im2);
Sz= sz*decim;
imhf = zeros(Sz);

iml2f = kron(ones(decim),fft2(im2));
imhf = full(decim*iml2f.*F{1}{1});

for inres = 2:Lt+1

    % number of direction
    N = Cf(inres-1);
    % size of image at current resolution
    S = 2*size(im2);

    indec = Lt+1-inres;
    decim = 2^(indec);

    % imhf = zeros(S);

    for in = 1:N
        if in > N/2
            dm = decim*[2 N/3];
        else
            dm = decim*[N/3 2];
        end
        % sz = S./dm;

        % fft the subband, usample in freq domain
        imh2f = kron(ones(dm),fft2(y{inres}{in}));
        
        % filtering in freq. domain by multiplication
        imhf = imhf + sqrt(2*prod(dm))*imh2f.*F{inres}{in};
    end

end

im2 = real(ifft2(imhf) );



