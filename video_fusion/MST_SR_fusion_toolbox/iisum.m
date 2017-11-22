%---------------------------------------------------------
%---------------------------------------------------------
function sum = iisum(iimg,x1,y1,x2,y2)
if(x1>1 && y1>1)
    sum = iimg(y2,x2)+iimg(y1-1,x1-1)-iimg(y1-1,x2)-iimg(y2,x1-1);
elseif(x1<=1 && y1>1)
    sum = iimg(y2,x2)-iimg(y1-1,x2);
elseif(y1<=1 && x1>1)
    sum = iimg(y2,x2)-iimg(y2,x1-1);
else
    sum = iimg(y2,x2);
end