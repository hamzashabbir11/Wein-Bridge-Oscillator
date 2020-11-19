%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Wien-bridge Oscillator Exploration  %   
%       with Matlab Implementation     %
% Author: M.Sc. Eng. Hristo Zhivomirov %
%               12/22/2012             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear, clc, close all

Cprime = 1;  % uF
Rprime = 1;  % kOhms
R1prime = 1; % kOhms
R2prime = 1.999; % kOhms

C = (Cprime)*0.000001; % Converting positive feedback capacitance to Farads
R = (Rprime)*1000;     % Converting positive feedback resistance to Ohms
R1 = (R1prime)*1000;   % Converting forward gain resistance, R1, to Ohms
R2 = (R2prime)*1000;   % Converting forward gain resistance, R2, to Ohms

% Calculates and displays the oscillating frequency
wo = (1/(C*R));
disp(['Oscillating frequency for the system w0= ', num2str(wo)])

% Calculates and displays the open loop transfer function A*b = tf([(1+(R2/R1)) 0],[C*R 3 (1/(C*R))])
disp('Open loop gain A*b:')
L = tf([(1+(R2/R1)) 0],[C*R 3 (1/(C*R))])
zpk(L) % Displays zero pole gain form of the loop gain transfer function, L

% Calculates and displeys the closed loop transfer function, Af
disp('Closed loop gain Af:')
Af = (1 + R2/R1)/(1-L)
zpk(Af) % Displays zero-pole-gain form of the closed loop transfer function, Af

% Transient step response of the close loop transfer function Af
figure(1)
step(Af)
title(['Oscillator Step Response for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])
grid on

% Transient impulse response of the close loop transfer function Af
figure(2)
impulse(Af)
title(['Oscillator Impulse Response for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])
grid on

% Bode plot of the closed loop transfer function
figure(3)
bode(Af);
title(['Oscillator Bode Plot for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])
grid on

% Nyquist diagram of the open loop transfer function
figure(4)
nyquist(Af);
title(['Nyquist Plot for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])

% Nyquist diagram of the open loop transfer function (polar form)
[mag, phase] = bode(Af);
mag = mag(1, :);
phase = phase(1, :);
phase = phase.*pi/180;
mag = interp(mag, 10);
phase = interp(phase, 10);
figure(5)
polar(phase, mag)
hold on
polar(0, -1, 'xr')
title(['Nyquist Plot for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])

% Pole-Zero plot of the close loop transfer function
figure(6)
pzplot(Af);
title(['Oscillator Pole-Zero Plot for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])

% 3D Pole-Zero plot of the close loop transfer function
figure(7)
pz=[pole(Af); zero(Af)];

rmin = min(real(pz));
rmax = max(real(pz));
imin = min(imag(pz));
imax = max(imag(pz));
res = 100;
a=0.2;
ax=[rmin-(rmax-rmin)*a rmax+(rmax-rmin)*a imin-(imax-imin)*a imax+(imax-imin)*a];
x=linspace(ax(1),ax(2),res);
y=linspace(ax(3),ax(4),res);

[xx,yy] = meshgrid(x,y);
transf_func3D = polyval(Af.num{1},xx+1i*yy)./polyval(Af.den{1},xx+1i*yy);
surfc(x, y, 20*log10(abs(transf_func3D)))
shading('interp')
colorbar
alpha(0.7)
xlabel('Real Axis');
ylabel('Imaginary Axis');
zlabel('Magnitude, dB');
title(['Oscillator Pole-Zero 3D Plot for R1=',num2str(R1prime),'kOhm',' and R2=',num2str(R2prime),'kOhm'])