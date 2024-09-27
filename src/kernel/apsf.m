function [amtf, apsf] = apsf(spatial_frequency, lambda, fNum, Wpp)
% An aberattion transfer function (ATF) is used in composition with the
% diffraction-limited MTF to compose a aberration PSF

dlmtf_fun = @(F, Fcut) 2/pi*(acos(F./Fcut)-F./Fcut.*sqrt(1-(F./Fcut).^2));  %diffraction-limited MTF
atf_fun = @(F, Fcut, Wpp) 1 - (Wpp/3.5/0.18)^2*(1 - 4*(F./Fcut - 0.5).^2); %aberattion transfer function ATF

Fcut_fun = @(lambda, fNum) 1./(lambda*fNum);

amtf = dlmtf_fun(spatial_frequency, Fcut_fun(lambda, fNum)).*atf_fun(spatial_frequency, Fcut_fun(lambda, fNum), Wpp);

apsf = fft(amtf);

end