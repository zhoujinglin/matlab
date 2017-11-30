function im = ucurvrec3d_s(ydec, Cf, F2, ind2, cfind )
% UCURVREC3D_S  3-d ucurvelet reconstruction using sparse 3-D window
%
%	im = ucurvrec3d_s(ydec, Cf, F2, ind2, cfind )
%
% Input:
%   Sz : size of the generated window
%   Cf : number of directional curvelet. In uniform curvelet, n is 3*2^n
%   r : parameter of the meyer window used in diameter direction
%   alpha : paramter for meyer window used in angle function
%
% Output:
%   im : Reconstructed image
%
% Example
%
% Sz = [32 64 128];
% Cf = [3 6];
% r = pi*[0.3 0.5 0.85 1.15];
% alpha = 0.15;
% F2 = ucurvwin3d(Sz, Cf, r, alpha);
% im = rand(Sz);
% [F2, ind, cf] = ucurvwin3d_s(Sz, 6, r, alpha);
% ydec = ucurvdec3d_s(im, 6, F2, ind, cf );
% imr = ucurvrec3d_s(ydec, 6, F2, ind, cf );
%
% See also UCURVDEC3D_S, UCURVWIN3D_S
%

% tic
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
sbind = 1;
Fc = ucurwin_extract(F2, ind2, cfind, sbind, Sz);
imhf = (1./sqrt(2))*upsamp_filtfft(ydec{1}, Fc, decim);

% for each resolution estimation
Lcol = size(cfind, 1);
for incol =  2: Lcol
    Fc = ucurwin_extract(F2, ind2, cfind, incol, Sz);
    % resolution index
    inres = cfind(incol, 2);
    % parameter for the resolution ---------------------------------------
    % current size of image at resoltuion
    cs = SzL.*2^(inres-1);

    ncf = Cf(:,inres - 1);
    in1 = cfind(incol, 3);
    in2 = cfind(incol, 4);
    dir = cfind(incol, 5);

    % decimation ratio
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

    imhf = imhf + ...
        upsamp_filtfft(ydec{inres}{dir}{in1,in2}, Fc, decim2);
end

im = real(ifftn(imhf));

end

% =========================================================================
function interp_fft = upsamp_filtfft(subband, filtfft, decim)
% toc
% size of subband
% disp('processing subband ...');

sz1 = size(subband);
% decimation ratio
% decim2 = decim(dir,:);

% upsampled in freq domain
tmp = fftn(subband);
interp_fft = repmat(tmp, decim);
clear tmp;

Dec = sqrt(2*prod(decim));

% filtering
interp_fft = Dec*interp_fft.*filtfft;
% tic
end

% % =========================================================================
% function ucurwin_extract(F2, ind2, cfind, sbind, Sz)
% tmp_cs_com = single(zeros(Sz));
% if (sbind>0)
%     st = 1;
% else
%     st = cfind(sbind-1, 1)+1;
% end
% 
% ind_tmp = ind2(st:cfind(sbind, 1));
% val_tmp = F2(st:cfind(sbind));
% tmp_cs_com(ind_tmp) = val_tmp;
% end
% % =========================================================================

% =========================================================================
function Fc = ucurwin_extract(F2, ind2, cfind, sbind, Sz)
Fc = zeros(Sz);
if (sbind==1)
    st = 1;
else
    st = cfind(sbind-1, 1)+1;
end
ind_tmp = ind2(st:cfind(sbind, 1));
val_tmp = F2(st:cfind(sbind));

Fc(ind_tmp) = val_tmp;
end



