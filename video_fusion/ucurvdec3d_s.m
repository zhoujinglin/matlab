function [ydec, cf_ydec] = ucurvdec3d_s(im, Cf, F2, ind2, cfind )
% UCURVDEC3D_S   3-d ucurvelet decomposition using sparse 3-D window
%
%	 ydec = ucurvdec3d_s(im, Cf, F2, ind2, cfind )
%
% Input:
%   Sz : size of the generated window
%   Cf : number of directional curvelet. In uniform curvelet, n is 3*2^n
%   r : parameter of the meyer window used in diameter direction
%   alpha : paramter for meyer window used in angle function
%
% Output:
%   FL: 2-D window of low pass function for lower resolution curvelet
%
% Example
%
% See also FUN_MEYER
%
% History
% tic
% 
%
if size(Cf,1) == 1
    Cf = kron(ones(3,1), Cf);
end

% the multiresolution length
Lt = size(Cf,2);
% Size of window
im = single(im);
% 
Sz = size(im);

% image for processing
% cim = im;

fim3 = fftn(im);
clear im;

% low pass band ----------------------------------------------------------
sbind = 1;
% Fc = ucurwin_extract(F2, ind2, cfind, sbind, Sz);
max_ind2 = cfind(end,1);
Fc2 = fim3(ind2(1:max_ind2)).*F2(1:max_ind2);  

clear fim3 F2;
% tmp_cs_com = 
ucurwin_extract(Fc2, ind2, cfind, sbind, Sz);

decim = 2^(Lt)*[1 1 1];

ydec{1} = (1./sqrt(2))*filter_decim(tmp_cs_com, decim);

Lcol = size(cfind, 1);

for incol =  2: Lcol
    % incol;    
    % Fc = ucurwin_extract(F2, ind2, cfind, incol, Sz);
    % resolution index
    inres = cfind(incol, 2);
    % parameter for the resolution ---------------------------------------
    % current size of image at resoltuion
    cs = Sz./2^(inres-1);
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


    % tmp_cs_com = 
    ucurwin_extract(Fc2, ind2, cfind, incol, Sz);

    ydec{inres}{dir}{in1,in2} = filter_decim(tmp_cs_com,  decim2);
%     tmp = filter_decim(tmp_cs_com,  decim, dir);
%     ltmp = prod(size(tmp));
%     cf_ydec(incol,1) = cf_ydec(incol-1,1)+ltmp;
%     
%     % resolution - dir - dir - pyramid
%     cf_ydec(incol,2:5) = cfind(incol+1,2:5);
%     ydec(cf_ydec(incol-1,1)+1:cf_ydec(incol,1)) = tmp(:);
end

% =========================================================================
function ucurwin_extract(F2, ind2, cfind, sbind, Sz)
tmp_cs_com = single(zeros(Sz));
if (sbind==1)
    st = 1;
else
    st = cfind(sbind-1, 1)+1;
end

ind_tmp = ind2(st:cfind(sbind, 1));
val_tmp = F2(st:cfind(sbind,1));
tmp_cs_com(ind_tmp) = val_tmp;
end
% =========================================================================

% % =========================================================================
% function ucurwin_extract(F2, ind2, cfind, sbind, Sz)
% tmp_cs_com = single(zeros(Sz));
% ind_tmp = ind2(cfind(sbind, 1)+1:cfind(sbind+1, 1));
% val_tmp = F2(cfind(sbind)+1:cfind(sbind+1));
% tmp_cs_com(ind_tmp) = val_tmp;
% end
% % =========================================================================

end

% =========================================================================
function subband = filter_decim(tmp, decim)
% toc
% dir

% disp('processing subband ...');

% tmp = datfft.*filtfft;

% size of image
sz1 = size(tmp);
% decimation ratio
% decim2 = decim(dir,:);
% size of subband and shifting step
sz2 = sz1./decim;

% decim in freq domain
tmp2 = zeros(sz2);
for in1 = 0:decim(1)-1
    for in2 = 0:decim(2)-1
        for in3 = 0:decim(3)-1
            p1 = [in1, in2, in3].*sz2+1 ;
            p2 = ([in1, in2, in3]+1).*sz2;
            tmp2 = tmp2 + ...
                    tmp(p1(1):p2(1), p1(2):p2(2), p1(3):p2(3));

        end
    end
end

subband = single(sqrt(2/prod(decim)).*ifftn(tmp2));
% tic
end
% 
% % =========================================================================
% function Fc = ucurwin_extract(F2, ind2, cfind, sbind, Sz)
% Fc = single(zeros(Sz));
% ind_tmp = ind2(cfind(sbind, 1)+1:cfind(sbind+1, 1));
% val_tmp = F2(cfind(sbind)+1:cfind(sbind+1));
% Fc(ind_tmp) = val_tmp;
% end


