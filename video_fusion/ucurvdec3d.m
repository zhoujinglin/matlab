function ydec = ucurvdec3d(im, Cf, F)
% UCURVDEC3D   3-d ucurvelet decomposition using normal 3-D window
%
%	 ydec = ucurvdec3d_s(im, Cf, F)
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
% See also FUN_MEYER
%
% History
% 
% tic
if size(Cf,1) == 1
    Cf = kron(ones(3,1), Cf);
end

% the multiresolution length
Lt = size(Cf,2);
% Size of window
Sz = size(im);

% image for processing
% cim = im;

fim3 = fftn(im);

% low pass band -----------------------------------------------------------
FL = F{1}{1};
tmp_cs_com = fim3.*FL;
decim = 2^(Lt)*[1 1 1];
ydec{1} = (1./sqrt(2))*filter_decim(tmp_cs_com, decim);

clear FL n;
clear tmp_cs_real tmp_cs_com;

% for each resolution estimation
for inres = 2:Lt+1
    % parameter for the resolution ----------------------------------------
    % current size of image at resoltuion
    cs = Sz./2^(inres-1);
    ncf = Cf(:,inres - 1);
    
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

        %         % decimation ratio
        %         switch n
        %             case {12}
        %                 decim = 2^(Lt+1-inres)*[2 8 8; 8 2 8; 8 8 2];
        %                 Dec = 16;
        %             case {6}
        %                 decim = 2^(Lt+1-inres)*[2 4 4; 4 2 4; 4 4 2];
        %                 Dec = 8;
        %             case {3}
        %                 decim = 2^(Lt+1-inres)*[2 2 2; 2 2 2; 2 2 2];
        %                 Dec = 4;
        %         end
        %
        % FFT image at current resolution
        % fim3 = fftn(cim);
        % clear cim;

        % the main loop ---------------------------------------------------
        % in = 0;

        for in1 = 1:n1
            for in2 = 1:n2
                Fc = F{inres}{dir}{in1,in2} ;

                tmp_cs_com = fim3.*Fc;
                ydec{inres}{dir}{in1,in2} = filter_decim(tmp_cs_com,decim2);

            end
        end
    end
end

% ydec{1} = cim;

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
% size of subband and shifting ste
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

% =========================================================================



