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
clear all;
clf;

dt = 0.001;
dx = 1.0;
D = 1.0;
K = 1.0;


beta1 = D*dt/(dx*dx);
beta2 = 2*K*D*dt/(dx*dx*dx*dx);

n = 32;
N = 40000;
m = 1;

hold on;

oldc = zeros(n,1);
newc = zeros(n,1);
g = zeros(n,1);
for i = 1:n
  oldc(i) = 0.5*(1 + sin(2*pi*i*m/n));
endfor
plot (oldc,'b');
A = 1;
for k = 1:n
  g(k) = 2*A*(6*oldc(k)*oldc(k) - 6*oldc(k) + 1);
endfor

for j = 1:N
  for i = 1:n
    for k = 1:n
       g(k) = 2*A*(6*oldc(k)*oldc(k) - 6*oldc(k) + 1);
    endfor
    w = i-1;
    ww = i-2;
    e = i+1;
    ee = i+2;
    if(w==0) w = n;
    endif
    if(ww==0) ww = n;
    endif
    if(ww==-1) ww = n-1;
    endif
    if(e == n+1) e = 1;
    endif
    if(ee == n+1) ee = 1;
    endif
    if(ee == n+2) ee = 2;
    endif
    newc(i) = oldc(i) + beta1*g(i)*(oldc(w)+oldc(e) - 2*oldc(i)) - beta2*(oldc(ww) - 4*oldc(w) + 6*oldc(i) - 4*oldc(e) + oldc(ee) );
  endfor
  oldc = newc; 
endfor
plot(oldc,'r');
 
 
    
    
    
  
