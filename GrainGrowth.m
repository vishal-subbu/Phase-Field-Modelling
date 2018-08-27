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
Nx = 64;Ny = 64;
D = 1;
halfNx= Nx/2;
halfNy = Ny/2;
phi = zeros(Nx,Ny);
phihat = zeros(Nx,Ny);
r_2 = zeros (150,1);
r = 15.0;
L = 1; 
a = 1; 
for i = 1:Nx
  for j = 1:Ny
    if (((i-halfNx)*(i-halfNx) + (j-halfNy)*(j-halfNy)) < r*r)
    phi(i,j) = 1.0;
    elseif phi(i,j) = 0.0;
    endif  
  endfor
endfor
#plot the initial profile 
#mesh(phi);
#view(2);
#pause(1);

phiprofile = phi(:,halfNy);
for i = 2:Nx
  if((phiprofile(i-1)<1) && (phiprofile(i) == 1.0))
     r1 = i;
  endif
  if((phiprofile(i-1) == 1.0) && (phiprofile(i) < 1))
     r2 = i-1;
  endif
endfor
r2_(a,1) = (r2-r1)*(r2-r1) /4.0;
a = a+1 ;
delkx= 2*pi/Nx;
delky = 2*pi/Ny;
plot (phiprofile);
hold on;

A = 1.0;
M = 1.0;
kappa = 1.0;
L =1.0;
for m = 1:15
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
    phiprofile = phi(:,halfNy);
    for b = 2:Nx
      if((phiprofile(b-1)<0.5) && (phiprofile(b) >=0.5))
        r1 = b;
      endif
      if((phiprofile(b-1) >=0.5) && (phiprofile(b) < 0.5))
         r2 = b-1;
      endif
    endfor
    r_2(a,1) = (r2-r1)*(r2-r1)/4.0;
    a = a+1 ;   
   endfor
   #mesh(phi);
   #view(2);
   #pause(0);
   plot (phiprofile);
endfor
hold off;
plot (r_2,'r.');
print -depsc GrainGrowth.eps

printf("Total CPU time : %f seconds\n",(cputime-t));
