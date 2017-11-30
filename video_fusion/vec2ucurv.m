function y = vec2ucurv(yind, mark)
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
tmp = yind(1:mark(1,1));

y{1} = reshape(tmp, mark(1,2), mark(1,3));
% band index

for min = 2:size(mark,1) % for each consider band
	%	mark(min-1,1)+1
	%	mark(min,1)
	%	yind(mark(min,1))
	tmp = yind(mark(min-1,1)+1:mark(min,1));

        y{mark(min,4)}{mark(min,5)} = reshape(tmp,mark(min,2),mark(min,3));
        
end

