%
% Copyright (c) 2018, Vishal_S
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Title: Computational Materials Thermodynamics
%
% Developer: Vishal S
%
% Contact Info: vishalsubbu97@gmail.com
%

intensity = 50.0;
pulse_width = 0.05;
pulse_on = 0.02;
n = 50;
length = 200;
a_cos = 1:n;
a_sin = 1:n;
c = pi/pulse_width;
a0 = 0;
t = linspace(0,pulse_width,length);
func = linspace(0,pulse_width,length);

for i = 1:length
    if(t(i)<pulse_on)
        func(i) = intensity;
    else 
        func(i) = 0.0;
    end
end

a0 = (intensity*pulse_on)/(2.0*pulse_width);

for i = 1:n 
    a_cos(i) = (intensity/pulse_width)*(sin(pulse_on*c*i))/(i*c);
    a_sin(i) = (intensity/pulse_width)*(-cos(pulse_on*c*i)+1.0)/(i*c);
end
func_num = 1:length;
for i = 1:length
    func_num(i) = a0;
    for j = 1:n
        func_num(i) = func_num(i) + a_cos(j)*cos(j*c*t(i)) + a_sin(j)*sin(j*c*t(i));
    end
end

plot(t,func);
hold on;
plot(t,func_num);
hold off;


