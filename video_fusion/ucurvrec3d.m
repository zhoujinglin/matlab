function im = ucurvrec3d(ydec, Cf, F)
% UCURVREC3D  3-d ucurvelet reconstruction using normal 3-D window
%
%	im = ucurvrec3d_s(ydec, Cf, F)
%
% Input:
%   Sz : size of the generated window
%   Cf : number of directional curvelet. In uniform curvelet, n is 3*2^n
%   r : parameter of the meyer window used in diameter direction
%   alpha : paramter for meyer window used in angle function
%
% Output:
%   FL: 2-D window of low pass function for lower resolution curvelet
%   F : cell of size n containing 2-D matrices of size s*s
%
% Example
%
% S = 512;
% Cf = [6 6];
% alpha = 0.15;
% r = pi*[0.3 0.5 0.85 1.15];
% Sz = 64;Lt = 2;
% F = ucurv_win(S, N, r, alpha);
%
% See also FUN_MEYER
%
% tic
% the multiresolution length
if size(Cf,1) == 1
    Cf = kron(ones(3,1), Cf);
end

% the multiresolution length
Lt = size(Cf,2);

% Size of image at lowest resolution
SzL = size(ydec{1});
Sz = 2^Lt*SzL;
% locate im
imhf = zeros(SzL.*2^(Lt));

decim = 2^Lt*[1 1 1];
dir = 1;
Dec = sqrt(prod(decim));

% low resolution image
imhf = (1./sqrt(2))*upsamp_filtfft(ydec{1}, F{1}{1}, decim);

% for each resolution estimation
for inres = 2:(Lt+1)
    % parameter for the resolution ---------------------------------------
    % current size of image at resoltuion
    cs = SzL.*2^(inres-1);
    ncf = Cf(:,inres - 1);
    % tmp_cs_com = zeros(cs+1)+sqrt(-1)*zeros(cs+1);

    for dir = 1:3
        switch dir
            case {1}
                n1 = ncf(2);
                n2 = ncf(3);
                decim2 = 2^(Lt+1-inres)*[2, 2*n1/3, 2*n2/3];
            case {2}
                n1 = ncf(1);
                n2 = ncf(3);
                decim2 = 2^(Lt+1-inres)*[2*n1/3, 2, 2*n2/3];
            case {3}
                n1 = ncf(1);
                n2 = ncf(2);
                decim2 = 2^(Lt+1-inres)*[2*n1/3, 2*n2/3, 2];
            otherwise
                disp('Error');
        end

        %     switch n
        %         case {12}
        %             decim = 2^(Lt+1-inres)*[2 8 8; 8 2 8; 8 8 2];
        %             Dec = 16;
        %         case {6}
        %             decim = 2^(Lt+1-inres)*[2 4 4; 4 2 4; 4 4 2];
        %             Dec = 8;
        %         case {3}
        %             decim = 2^(Lt+1-inres)*[2 2 2; 2 2 2; 2 2 2];
        %             Dec = 4;
        %     end

        % the main loop ------------------------------------------------------

        for in1 = 1:n1
            for in2 = 1:n2
                Fc = F{inres}{dir}{in1,in2} ;

                imhf = imhf + ...
                    upsamp_filtfft(ydec{inres}{dir}{in1,in2}, Fc, decim2);

            end
        end
    end

    % im = real(im3a+im3b);
end

im = real(ifftn(imhf));

end

% =========================================================================
function interp_fft = upsamp_filtfft(subband, filtfft, decim)
% toc
% size of subband
% disp('processing subband ...');

sz1 = size(subband);

% upsampled in freq domain
tmp = fftn(subband);
interp_fft = repmat(tmp, decim);
clear tmp;

Dec = sqrt(2*prod(decim));

% filtering
interp_fft = Dec*interp_fft.*filtfft;
% tic
end

% =========================================================================


