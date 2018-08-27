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
clear all
clf
clc
more off;

function y = Xi(e1,e2)
  y = 0.25*e1*e1*e1*e1 - 0.5*e1*e1 + 0.25*e2*e2*e2*e2 - 0.5*e2*e2 + 2*e1*e1*e2*e2 + 0.25;
endfunction

function y = Xiprime(e1,e2)
  y = e1*e1*e1 - e1 + 4*e1*e2*e2;
endfunction

kappac = 1.0;
kappaphi1 = 1.0/3.0;
kappaphi2 = 1.0/3.0;
A = 1.0;
B = 1.0;
Z = 1.0;
N = 64;
dx = 1.0;
c = zeros(N,N);
phi1 = zeros(N,N);
phi2 = zeros(N,N);

for i = 1:N
  for j = 1:N
    if((j<N/4) || (j>3*N/4))
       c(i,j)  =1.0;
       if(i<N/2) phi1(i,j) = 1.0;
       endif
       if(i>N/2) phi2(i,j) = 1.0;
       endif
    endif
  endfor
endfor
mesh(c);
view(2);
pause(0);
halfN = N/2;
delk = 2*pi/N;
dt = 0.5;
for n = 1:50 
  for m = 1:2
    for i = 1:N
      for j = 1:N
        gc(i,j) = 2*A*c(i,j)*(1-c(i,j))*(1-2*c(i,j)) + 2*B*c(i,j)*Xi(phi1(i,j),phi2(i,j));
        gc(i,j) = gc(i,j) - 2*Z*(1-c(i,j))*(phi1(i,j)*phi1(i,j) + phi2(i,j)*phi2(i,j));
        gphi1(i,j) = B*c(i,j)*c(i,j)*Xiprime(phi1(i,j),phi2(i,j)) + 2*Z*(1-c(i,j))*(1-c(i,j))*phi1(i,j);
        gphi2(i,j) = B*c(i,j)*c(i,j)*Xiprime(phi2(i,j),phi1(i,j)) + 2*Z*(1-c(i,j))*(1-c(i,j))*phi2(i,j);
      endfor
    endfor
    gchat = fft2(gc);
    gphi1hat = fft2(gphi1);
    gphi2hat = fft2(gphi2);
    chat = fft2(c);
    phi1hat = fft2(phi1);
    phi2hat = fft2(phi2);
     for  i = 1:N
       if((i-1)<=halfN) kx = (i-1)*delk;
       elseif ((i-1)>halfN) kx = (i-1-N)*delk;
       endif
       for j = 1:N
         if((j-1)<=halfN) ky = (j-1)*delk;
         elseif ((j-1)>halfN) ky= (j-1-N)*delk;
         endif
         k2 = kx*kx + ky*ky;
         k4 = k2*k2;
         chat(i,j) = (chat(i,j) - dt*k2*gchat(i,j))/(1 + 2*kappac*k4*dt);
         phi1hat(i,j) = (phi1hat(i,j) - dt*gphi1hat(i,j))/(1+2*kappaphi1*k2*dt);
         phi2hat(i,j) = (phi2hat(i,j) - dt*gphi2hat(i,j))/(1+2*kappaphi2*k2*dt);
      endfor
    endfor
    c = real(ifft2(chat));
    phi1 = real(ifft2(phi1hat));
    phi2 = real(ifft2(phi2hat));
  endfor
  mesh(c);
  view(2);
  pause(0);
endfor;

       
