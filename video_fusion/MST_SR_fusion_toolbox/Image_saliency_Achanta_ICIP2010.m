%---------------------------------------------------------
% Copyright (c) 2012 Radhakrishna Achanta [EPFL]
% Contact: firstname.lastname@epfl.ch
%---------------------------------------------------------
% Citation:
% @InProceedings{Achanta_Saliency_ICIP_2010,
%    author      = {Achanta, Radhakrishna and S黶strunk, Sabine},
%    booktitle   = {{IEEE} {I}nternational {C}onference on {I}mage
%                  {P}rocessing},
%    year        = 2010
% }
%---------------------------------------------------------
%
%
%---------------------------------------------------------
% Read image and blur it with a 3x3 or 5x5 Gaussian filter
% 
% The MATLAB version of the code is much slower than the C++
% code provided here:
% http://ivrgwww.epfl.ch/supplementary_material/RK_ICIP2010/index.html
% or here:
% http://ivrg.epfl.ch/research/saliency/MSSS.html
% The MATLAB code should be used for understanding
% purposes. For comparisons one should use the C++ code.
%
% For comments, bugs, and questions, please write to:
% firstname.lastname@epfl.ch, where
% firstname = radhakrishna and lastname = achanta.
%---------------------------------------------------------
function SalientMap=Image_saliency_Achanta_ICIP2010(grayimg)  %grayimg 输入一个图像

% inputfilename = 'filename.jpg';
% img = imread(inputfilename);%Provide input image path

% colorimg= imread('orchid_.bmp'); % color
% colorSaliencyMap = Color_Saliency(colorimg);

graySaliencyMap = Gray_Saliency(grayimg);
% graySaliencyMap = mat2gray(graySaliencyMap,[0 255]);
% figure;imshow(graySaliencyMap,[]);
SalientMap= graySaliencyMap;



function GraySaliencyMap=Gray_Saliency(grayimage)
% img= imread('source24_1.tif');
img=grayimage;
dim = size(img);           %例如256 256，270 260
width = dim(2);height = dim(1);    % height是高，即：行，270
gfrgb = imfilter(img, fspecial('gaussian', 3, 3), 'symmetric', 'conv');   %滤波结果
                     %'symmetric'  Input array values outside the bounds of the array
                     %are computed by mirror-reflecting the array across
                     %the array border.
                     
l = double(img(:,:,1));

%---------------------------------------------------------
% Create integral images
%---------------------------------------------------------
li = cumsum(cumsum(l,2));
 
%---------------------------------------------------------
% Compute Lab average values (note that in the paper this
% averages are found from the unblurred original image, but
% the results are quite similar)
%---------------------------------------------------------
sm = zeros(height, width);
sm2 = zeros(height, width);
for j = 1:height
    yo = min(j, height-j);
    y1 = max(1,j-yo); 
    y2 = min(j+yo,height);
    for k = 1:width
        xo = min(k,width-k);
        x1 = max(1,k-xo); 
        x2 = min(k+xo,width);
        invarea = 1.0/((y2-y1+1)*(x2-x1+1));
        lm = iisum(li,x1,y1,x2,y2)*invarea;
        %---------------------------------------------------------
        % Compute the saliency map
        %---------------------------------------------------------
        sm(j,k) = (l(j,k)-lm).^2 ;
    end
end
% figure;imshow(sm,[]);

GraySaliencyMap=sm;

function ColorSaliencyMap=Color_Saliency(colorimage)
% img= imread('orchid_.bmp');
img=colorimage;
dim = size(img);
width = dim(2);height = dim(1);
gfrgb = imfilter(img, fspecial('gaussian', 3, 3), 'symmetric', 'conv');
%---------------------------------------------------------
% Perform sRGB to CIE Lab color space conversion (using D65)
%---------------------------------------------------------
cform = makecform('srgb2lab');%, 'whitepoint', whitepoint('d65'));
lab = applycform(gfrgb,cform);
l = double(lab(:,:,1));
a = double(lab(:,:,2));
b = double(lab(:,:,3));
%[l a b] = RGB2Lab(gfrgb(:,:,1),gfrgb(:,:,2), gfrgb(:,:,3));
%---------------------------------------------------------
% Create integral images
%---------------------------------------------------------
li = cumsum(cumsum(l,2));
ai = cumsum(cumsum(a,2));
bi = cumsum(cumsum(b,2));     
%---------------------------------------------------------
% Compute Lab average values (note that in the paper this
% averages are found from the unblurred original image, but
% the results are quite similar)
%---------------------------------------------------------
sm = zeros(height, width);
sm2 = zeros(height, width);
for j = 1:height
    yo = min(j, height-j);
    y1 = max(1,j-yo); y2 = min(j+yo,height);
    for k = 1:width
        xo = min(k,width-k);
        x1 = max(1,k-xo); x2 = min(k+xo,width);
        invarea = 1.0/((y2-y1+1)*(x2-x1+1));
        lm = iisum(li,x1,y1,x2,y2)*invarea;
        am = iisum(ai,x1,y1,x2,y2)*invarea;
        bm = iisum(bi,x1,y1,x2,y2)*invarea;
        %---------------------------------------------------------
        % Compute the saliency map
        %---------------------------------------------------------
        sm(j,k) = (l(j,k)-lm).^2 + (a(j,k)-am).^2 + (b(j,k)-bm).^2;
    end
end
% figure;imshow(sm,[]);
% img = (sm-min(sm(:)))/(max(sm(:))-min(sm(:)));
% figure;imshow(img,[]);
% [pathstr,name,ext,versn] = fileparts(inputfilename);
% imwrite(img,strcat(name,'_msss.jpg'));
ColorSaliencyMap=sm;
%---------------------------------------------------------
