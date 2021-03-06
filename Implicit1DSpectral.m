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
clear 
clf

t = cputime;
delt = 1.0;
N = 128;
D = 1;
c = zeros(N,1);
ctilde = zeros(N,1);
M = 2;
for i = 1:N
  c(i) = 0.5*(1 + sin(2*pi*M*i/N));
endfor
#plot the initial profile 
plot (c , 'r-');
hold on;
halfN = N/2;
delk = 2*pi/N;
for j= 1:20
  for m = 1:500
    ctilde = fft(c) ;
    for i = 1:N
      if(i<halfN)
        k = i*delk;
      endif
      if(i>=halfN)
        k = (i-N)*delk;
      endif
      ctilde(i) = ctilde(i)/(1+(D*k*k*delt));
    endfor
  endfor
  c = real(ifft(ctilde));
  plot(c,'g.');
endfor
print -depsc Implicit1DSpectral.eps

printf("Total CPU time : %f seconds\n",(cputime-t));


    
      

  
