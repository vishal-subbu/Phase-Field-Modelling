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
delt = 0.5;
N = 128;
D = 1;
c = zeros(N,1);
ctilde = zeros(N,1); 
for i = (N/4) + 1: 3*N/4
  c(i) = 1.0
endfor
#plot the initial profile 
plot (c , 'r-');
hold on;
halfN = N/2;
delk = 2*pi/N;
for j= 1:4000
  g = 2.*c.*(1-c).*(1-2.*c);
  ghat = fft(g);
  chat = fft(c) ;
  for i = 1:N
    if((i-1)<=halfN)
      k = (i-1)*delk;
    endif
    if((i-1)>halfN)
      k = (i-1-N)*delk;
    endif
    k2 = k*k;
    k4 = k2*k2;
    chat(i) = (chat(i) - (k2*delt*ghat(i)))/(1+(2*k4*delt));
  endfor
  c = real(ifft(chat));
endfor
plot(c,'g');
chat = fft(c);
for i = 1:N
   if((i-1)<=halfN)
     k = (i-1)*delk;
   endif
   if((i-1)>halfN)
     k = (i-1-N)*delk;
   endif
   chat (i) = complex(0,1) *k*chat(i);
endfor
cprime = real(ifft(chat));
energy1 -0.0;
energy2 = 0.0;
kappa  = 1.0;
A = 1.0;
for i = 1:N
  energy1 = energy1 + A*c(i)*c(i)*(1-c(i))*(1-c(i));
  energy2 = energy2 + kappa*cprime(i)*cprime(i);
  endfor
energy1
energy2
0.5*(energy1 + energy2)

  

 
print -depsc CH1DFT.eps
