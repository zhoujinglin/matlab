function [F2, ind2, cf] = ucurvwin3d_s(Sz, Cf, r, alpha)
% UCURVWIN3D_S   Generate the sparse curvelet windows that used in 3-D 
% uniform curvelet inverse and foward transform
%
%	[F2, ind2, cf] = ucurvwin3d_s(Sz, Cf, r, alpha)
%
% Input:
%   Sz : size of the generated window
%   Cf : number of directional curvelet. In uniform curvelet, n is 3*2^n
%   r : parameter of the meyer window used in diameter direction
%   alpha : paramter for meyer window used in angle function
%
% Output:
%   F2 : a column of curvelet window value that are diff. from zero
%   ind : index of position correspond to value stored in F
%   cf  : configuration paramter, specify ending point of window value in F
%   and ind.
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
% running time 6*64 - 6 sec, 6*128 - 53 sec
%
% See also UCURVDEC3D_S UCURVREC3D_S
%
if size(Cf,1) == 1
    Cf = kron(ones(3,1), Cf);
end

% preallocate memory for F2 and ind2
Lcf = 1.5*prod(Sz);
F2 = single(zeros(Lcf, 1));
ind2 = uint32(zeros(Lcf, 1));
% tic
% the multiresolution length
Lt = size(Cf,2);

r = [r(1:2); r(3:4)];
for in = 1: Lt-1
    r = [0.5*r(1,:); r ];
end

% Size of window
if max(size(Sz)) == 1
    Sz = [Sz Sz Sz];
end

cs = Sz;

% create the grid
S1 = -1.5*pi:pi/(cs(1)/2):(0.5*pi-pi/(cs(1)/2));
S2 = -1.5*pi:pi/(cs(2)/2):(0.5*pi-pi/(cs(2)/2));
S3 = -1.5*pi:pi/(cs(3)/2):(0.5*pi-pi/(cs(3)/2));

% create the mesh for estimate the angle function v(T(theta)) ---------
[M31, M32] =  adapt_grid(S1, S2);
[M12, M13] =  adapt_grid(S2, S3);
[M21, M23] =  adapt_grid(S1, S3);


pack;

% for each resolution estimation
for inres = 2:Lt+1
    % parameter for the resolution ---------------------------------------
    % current size of image at resoltuion
    cs = Sz;
    % inres2 = Lt +2 - inres;
    ncf = Cf(:,inres-1);
    mncf = max(ncf);

    clear tmp_cs_real tmp_cs_com;

    % Low - high pass window ---------------------------------------------
    f1 = single(fun_meyer(abs(S1),[-2 -1 r(inres-1,:)]));
    f2 = single(fun_meyer(abs(S2),[-2 -1 r(inres-1,:)]));
    f3 = single(fun_meyer(abs(S3),[-2 -1 r(inres-1,:)]));
    % 3-d window
    FL3 = repmat(f1(:), [1, cs(2), cs(3)]).* ...
        permute(repmat(f2(:), [1, cs(1), cs(3)]),[2, 1, 3] ).* ...
        permute(repmat(f3(:), [1, cs(1), cs(2)]),[2, 3, 1]);
    
    f1 = single(fun_meyer(abs(S1),[-2 -1 r(inres,:)]));
    f2 = single(fun_meyer(abs(S2),[-2 -1 r(inres,:)]));
    f3 = single(fun_meyer(abs(S3),[-2 -1 r(inres,:)]));
    % 3-d window
    FR3 = repmat(f1(:), [1, cs(2), cs(3)]).* ...
        permute(repmat(f2(:), [1, cs(1), cs(3)]),[2, 1, 3] ).* ...
        permute(repmat(f3(:), [1, cs(1), cs(2)]),[2, 3, 1]);

    % high pass window
    FR3 = single(sqrt(FR3-FL3));
    
    % low pass band ------------------------------------------------------
    %     if inres == 2
    %         FL = sqrt(circshift(FL3(1:cs(1),1:cs(2),1:cs(3)), 0.25*cs ));
    %         F2{1}{inres-1} = FL;
    %         clear FL FL3;
    %     end

    % low pass band ------------------------------------------------------
    if inres == 2
        sbind = 1;
        
        FL = sqrt(circshift(FL3, 0.25*cs ));
        
        ind3 = find(FL);
        lind3 = length(ind3);
        % append index to ind2
        ind2(1: lind3) = uint32(ind3(:));
        % append value F2
        F2(1: lind3) = FL(ind3);
        % append configuration information
        % cf = [cf; cf(end,1)+lind3, inres, in1, in2, dir];
        % ind2 = find(FL);
        % F2 = FL(ind2);
        
        cf(sbind, 1) = lind3;
        cf(sbind, 2:5) = [1 1 1 1];
        % cf = cf(:);
        clear FL FL3;
    end


    % angle meyer window -------------------------------------------------
    %   angd = 4/(2*mncf);
    %   ang = angd*[-alpha alpha 1-alpha 1+alpha];
    %   angex = angd*[-2*alpha 2*alpha 1-2*alpha 1+2*alpha];
    alpha2 = alpha*4/(2*mncf);
    
    for in = 1:3
        n = ncf(in);
        angd = 4/(2*n);
        ang = [-alpha2,alpha2,(angd-alpha2),(angd+alpha2)];
        if (n == 3)
            n2 = 2;
        else
            n2 = n/2;
        end

        for in2 = 1:n2
            ang2 = -1+(in2-1)*angd+ang;

            switch in
                case {1}
                    fang{3}{1}{in2} = sqrt(single(fun_meyer(M21,ang2)));
                    fang{2}{1}{in2} = sqrt(single(fun_meyer(M31,ang2)));
                case {2}
                    fang{1}{1}{in2} = sqrt(single(fun_meyer(M32,ang2)));
                    fang{3}{2}{in2} = sqrt(single(fun_meyer(M12,ang2)));
                case {3}
                    fang{2}{2}{in2} = sqrt(single(fun_meyer(M13,ang2)));
                    fang{1}{2}{in2} = sqrt(single(fun_meyer(M23,ang2)));
            end
        end
    end

    % the very first filter require special handling
    % estimate the index of the point that need to be scale ---------------
    % first 3-D angle
    FA = repmat(fang{1}{1}{1},[1, 1, cs(3)]);
    FB = repmat(fang{1}{2}{1},[1, 1, cs(2)]);
    ang1 = FA.*permute(FB,[1, 3, 2]);
    % second 3-D angle
    FA = repmat(fang{2}{1}{1},[1, 1, cs(3)]);
    FB = repmat(fang{2}{2}{1},[1, 1, cs(1)]);
    ang2 = FA.*permute(FB,[3, 1, 2]);
    % third 3-D angle
    FA = repmat(fang{3}{1}{1},[1, 1, cs(2)]);
    FB = repmat(fang{3}{2}{1},[1, 1, cs(1)]);
    ang3 = permute(FA,[1, 3, 2]).*...
        permute(FB,[3, 1, 2]);
    ind = find(and( and ((ang1>0),  (ang2>0)), (ang3>0) ));
    G1ex = ang1.^2 + ang2.^2 + ang3.^2;
    G1ex = sqrt(G1ex(ind));
    
    clear ang1 ang2 ang3;
    % finish estimate ind and summation.
    
    clear n
    for in = 1:3
        if (ncf(in) == 3)
            ncf2(in) = 2;
        else
            ncf2(in) = ncf(in)/2;
        end
    end
    % the main loop ------------------------------------------------------
    % in = 0;
    for dir = 1:3
        switch dir
            case {1}
                for in1 = 1:ncf2(2)
                    for in2 = 1:ncf2(3)

                        FA = repmat(fang{1}{1}{in1},[1, 1, cs(3)]);
                        FB = permute(repmat(fang{1}{2}{in2},[1, 1, cs(2)]),[1, 3, 2]);
                        tmp_cs_real = FA.*FB.*FR3;
                        % now normalize those points in the angle wedge
                        tmp_cs_real(ind) = tmp_cs_real(ind)./G1ex;
                        clear FA FB;

                        F = circshift((tmp_cs_real(1:cs(1), 1:cs(2), 1:cs(3))), ...
                            [0.25*cs(1), 0.25*cs(2) , 0.25*cs(3)]);

                        % F2{inres}{dir}{in1,in2} = F;
                        addwindow2(F, inres, in1, in2, dir);
                        
                        if (in2~=ncf(3)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, 3);
                        % F2{inres}{dir}{in1, ncf(3)+1-in2} = Fc;
                        addwindow2(Fc, inres, in1, ncf(3)+1-in2, dir);
                        end

                        if (in1~=ncf(2)+1-in1)
                        Fc = rotate_ucurv3d_nest(F, 2);
                        % F2{inres}{dir}{ncf(2)+1-in1,in2} = Fc;
                        addwindow2(Fc, inres, ncf(2)+1-in1, in2, dir);
                        end
                        
                        if and(in1~=ncf(2)+1-in1, in2~=ncf(3)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, [3 ,2]);
                        % F2{inres}{dir}{ncf(2)+1-in1,ncf(3)+1-in2} = Fc;
                        addwindow2(Fc, inres, ncf(2)+1-in1, ncf(3)+1-in2, dir);
                        end

                    end
                end
            case {2}

                for in1 = 1:ncf2(1)
                    for in2 = 1:ncf2(3)

                        FA = repmat(fang{2}{1}{in1},[1, 1, cs(3)]);
                        FB = permute(repmat(fang{2}{2}{in2},[1, 1, cs(1)]),[3, 1, 2]);
                        tmp_cs_real = FA.*FB.*FR3;
                        % now normalize those points in the angle wedge
                        tmp_cs_real(ind) = tmp_cs_real(ind)./G1ex;
                        clear FA FB;

                        F = circshift((tmp_cs_real(1:cs(1), 1:cs(2), 1:cs(3))), ...
                            [0.25*cs(1), 0.25*cs(2) , 0.25*cs(3)]);

                        % F2{inres}{dir}{in1,in2} = F;
                        addwindow2(F, inres, in1, in2, dir);
                        
                        if (in2~=ncf(3)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, 3);
                        % F2{inres}{dir}{in1, ncf(3)+1-in2} = Fc;
                        addwindow2(Fc, inres, in1, ncf(3)+1-in2, dir);
                        end
                        
                        if (in1~=ncf(1)+1-in1)
                        Fc = rotate_ucurv3d_nest(F, 1);
                        % F2{inres}{dir}{ncf(1)+1-in1,in2} = Fc;
                        addwindow2(Fc, inres, ncf(1)+1-in1, in2, dir);
                        end
                        
                        if and(in1~=ncf(1)+1-in1, in2~=ncf(3)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, [3, 1]);
                        % F2{inres}{dir}{ncf(1)+1-in1,ncf(3)+1-in2} = Fc;
                        addwindow2(Fc, inres, ncf(1)+1-in1, ncf(3)+1-in2, dir);
                        end
                    end
                end

            case {3}
                
                for in1 = 1:ncf2(1)
                    for in2 = 1:ncf2(2)
                        FA = permute(repmat(fang{3}{1}{in1},[1, 1, cs(2)]),[ 1 3 2]);
                        FB = permute(repmat(fang{3}{2}{in2},[1, 1, cs(1)]),[3, 1, 2]);
                        tmp_cs_real = FA.*...
                            FB.*FR3;
                        % now normalize those points in the angle wedge
                        tmp_cs_real(ind) = tmp_cs_real(ind)./G1ex;
                        clear FA FB;

                        F = circshift((tmp_cs_real(1:cs(1), 1:cs(2), 1:cs(3))), ...
                            [0.25*cs(1), 0.25*cs(2) , 0.25*cs(3)]);

                        % F2{inres}{dir}{in1,in2} = F;
                        addwindow2(F, inres, in1, in2, dir);
                        
                        if (in2~=ncf(2)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, 2);
                        % F2{inres}{dir}{in1, ncf(2)+1-in2} = Fc;
                        addwindow2(Fc, inres, in1, ncf(2)+1-in2, dir);
                        end
                        
                        if (in1~=ncf(1)+1-in1)
                        Fc = rotate_ucurv3d_nest(F, 1);
                        % F2{inres}{dir}{ncf(1)+1-in1, in2} = Fc;
                        addwindow2(Fc, inres, ncf(1)+1-in1, in2, dir);
                        end
                        
                        if and(in1~=ncf(1)+1-in1, in2~=ncf(2)+1-in2)
                        Fc = rotate_ucurv3d_nest(F, [1, 2]);
                        % F2{inres}{dir}{ncf(1)+1-in1,ncf(2)+1-in2} = Fc;
                        addwindow2(Fc, inres, ncf(1)+1-in1, ncf(2)+1-in2, dir);
                        end
                    end
                end

            otherwise
                disp('Error switch');
        end
    end

end



% addwindow2 --------------------------------------------------------------
% nested function to handle sbind, F2, ind2, cf
    function addwindow2(tmp,  inres, in1, in2, dir)
        sbind = sbind+1;
        % index and length
        ind3 = find(tmp); % 
        lind3 = length(ind3);

        % append index to ind2
        ind2(cf(end,1)+1: cf(end,1)+lind3) = uint32(ind3(:));
        % append value F2
        F2(cf(end,1)+1: cf(end,1)+lind3) = tmp(ind3);
        % append configuration information
        cf = [cf; cf(end,1)+lind3, inres, in1, in2, dir];

    end
% ------------------------------------------------------------------------

end


%--------------------------------------------------------------------------
% utility functions 
%--------------------------------------------------------------------------

% flip 3-D matrix function ------------------------------------------------
function Fc = rotate_ucurv3d_nest(F, para)

Fc = F;

for in = 1: length(para)
    switch para(in)
        case {1}
            Fc = circshift(flipdim(Fc, 1), [1 0 0]);
        case {2}
            Fc = circshift(flipdim(Fc, 2), [0 1 0]);
        case {3}
            Fc = circshift(flipdim(Fc, 3), [0 0 1]);

        otherwise
            disp('Error');

    end
end

end

% create 2-D grid function-------------------------------------------------
function [M1, M2] = adapt_grid(S1, S2)

[x1, x2] = meshgrid(S2,S1);

% scale the grid approximate the tan theta function ------------------
% creat two scale grid for mostly horizontal and vertical direction

% first grid
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

end

% addwindow --------------------------------------------------------------
function [sbind, F2, ind, cf] = addwindow(sbind, F2, ind, cf, tmp,  inres, in1, in2, dir)
sbind = sbind+1;
ind2 = find(abs(tmp)>10^(-9));
ind = [ind;ind2];

F2 = [F2;tmp(ind2)];
cf = [cf;length(F2), inres, in1, in2, dir];
end


