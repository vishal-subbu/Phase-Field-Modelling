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

t = cputime;
delt = 0.5;
Nx = 128;Ny = 128;
D = 1;
phi = zeros(Nx,Ny);
phihat = zeros(Nx,Ny);
L = 1;  
for i = 1:Nx
  for j = 1:Ny
    phi(i,j) = 0.5 + 0.1*(0.5 - rand());
  endfor
endfor
%plot the initial profile 
mesh(phi);
view(2);
pause(1);

halfNx= Nx/2;
halfNy = Ny/2;
delkx= 2*pi/Nx;
delky = 2*pi/Ny;

A = 1.0;
M = 1.0;
kappa = 1.0;
L =1.0;
for m = 1:10
  for n = 1:10
    g = 2*A.*phi.*(1-phi).*(1-2.*phi);
    ghat = fft2(g);
    phihat = fft2(phi);
    for i = 1:Nx
      if (i<= halfNx ) kx = (i-1)*delkx;
      endif
      if (i > halfNx) kx = (i-1-Nx)*delkx;
      endif
      for j =1:Ny
        if (j<= halfNy ) ky = (j-1)*delky;
        endif
        if (j  > halfNy ) ky = (j-1-Ny)*delky;
        endif
        k2 = kx*kx + ky*ky;
        phihat(i,j) = (phihat(i,j) - L*delt*k2*ghat(i,j))/(1 + 2*L*kappa*k2*delt);
      endfor
    endfor
    phi = real(ifft2(phihat));     
   endfor
   mesh(phi);
   view(2);
   pause(0);
endfor
print -depsc AC2DFT.eps

printf("Total CPU time : %f seconds\n",(cputime-t));
