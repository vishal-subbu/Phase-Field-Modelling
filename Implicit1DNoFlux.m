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
clf
clear 

t  = cputime ;
%not dimensionalised numbers
delt = 0.1;
delx = 0.5;
D = 1.0;
n = 101;
alpha = D*delt/(delx*delx);

%initial conditions 

%c(1) = C_0 and the remoining is 0.0
c_old = zeros(n,1);
c_new = zeros(n,1);
c_old(1) = 1;
A = zeros(n,n);
%plot the initial profile
plot (c_old, 'r-;Initial profile;');

%populating the matrix
A(1,1) = 1;
for i=2:n
  A(i,i) = 1 + 2*alpha;
endfor
for i = 2:n-1
   A(i,i-1) = -alpha;
endfor
A(n,n-1) = -2*alpha;
for i = 3:n
  A(i-1,i) = -alpha;
endfor
  

%get handle 
ax = gca ;
set(ax, "linewidth",2.0);
axis("square");

5hold the plot 
hold on

% Compostion evolution
for k = 1:20 %this is for plotting 
  for j = 1:500
     c_new = inv(A)*c_old;
     c_old = c_new; 
  endfor # ending j loop
  plot (c_old);
endfor # ending k loop

% Save figure
print -djpg Implicit1D.jpg

printf('Total cpu time: %f seconds\n', cputime-t);
