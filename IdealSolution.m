%
% Copyright (c) 2018, Vishal_S
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Title: Phase field modelling
% 
% Developer: Vishal S
% 
% Contact Info: vishalsubbu97@gmail.com
%
temp = 300;
R = 8.314; 
x = 0.001:0.001:0.999;
Gibbs = 0:0.01:1;
Gibbs = R*temp*(x.*log(x) + (1-x).*log(1-x));
plot (x,Gibbs,'p');
