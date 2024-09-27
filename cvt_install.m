%% CVT INSTALLATION %%

% MATLAB environment
clear
clc
close all

if isfile("cvt_install.m")
    % SAVE INSTALL TO PATH
    copyfile("cvt_install.m",userpath)
    % SAVE PATH TREE
    writelines(["function cvt_home_path = cvt_home()",join(["cvt_home_path = '",pwd,"';"],''),"end"],'cvt_home.m')
    movefile('cvt_home.m',userpath)
end

% Path
addpath(genpath(cvt_home))
rmpath(genpath(fullfile(cvt_home,'debug')))

% DISPLAY
fprintf(['*** CVT installed. Have fun! ***\n'])