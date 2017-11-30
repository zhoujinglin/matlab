function [yind, mark] = ucurv2vec(y)
% PDTDFB2VEC   Convert the output of the PDTDFB into a vector form
%
%       [yind, cfig] = pdtdfb2vec(y)
%
% Input:
%   y:  an output of the PDTDFB
%
% Output:
%   yind :  1-D vector that contains all PDFB coefficients
%   mark :  starting point of each change in band in yind
%
% See also:	PDTDFBDEC, VEC2PDFB, (also PDFB2VEC in Contourlet toolbox)

% take out the directional subband complex amplitude value
tmp = y{1};
yind = tmp(:);
% band index
min = 1;
mark(min, 1) = prod(size(tmp));
mark(min,2) = size(tmp,1);
mark(min,3) = size(tmp,2);
mark(min,4) = 1;
mark(min,5) = 1;

for in = 2:length(y) % for each consider resolution
    
    for d = 1:length(y{in})
        min = min+1;
        tmp = y{in}{d};
        
        % first column is the ending point of the subband
        mark(min,1) = mark(min-1,1)+prod(size(tmp));
        % second column is the row size of the subband
        mark(min,2) = size(tmp,1);
         % third column is the column size of the subband
        mark(min,3) = size(tmp,2);
        % fourth column resolution the subband
        mark(min,4) = in;
        % fifth column direction the subband
        mark(min,5) = d;

        % [inc, inr] = meshgrid(1:Stmp(2), 1:Stmp(1));
        
        % 
        % tmp3 = [(tmp(:));
        
        yind = [yind; tmp(:)];
    end
end
