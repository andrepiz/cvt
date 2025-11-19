function range = opr2range(opr, Rbody, f, muPixel)
% Return the object-to-pixel ratio (OPR) from the range
% The two angles must be the same:
% atan(OPR/2/fpx) = asin(R/d)
% --> d = R/(sin(atan(opr/2/(f/muPixel))))

range = Rbody./sin(atan(opr/2./(f./muPixel)));
    
end