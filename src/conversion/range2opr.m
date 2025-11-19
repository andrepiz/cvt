function opr = range2opr(range, Rbody, f, muPixel)
% Return the object-to-pixel ratio (OPR) from the range
% The two angles must be the same:
% atan(OPR/2/fpx) = asin(R/d)

if nargin == 3
    opr = 2.*f.*tan(asin(Rbody./range));
elseif nargin == 4
    opr = 2.*(f./muPixel).*tan(asin(Rbody./range));
else
    error('Input the focal length in pixels or the focal length followed by the pixel pitch')
end
    
end