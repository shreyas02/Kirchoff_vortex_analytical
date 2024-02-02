clear 
clc

%% r = 10
rho = 1.3;
a = 1;
nu = 2 ;
c = 330;
e = 1.25e-3;
mach = 0.11;
t = 0 ;
theta = 0;
r = linspace(10, 300,100);;
o = (mach*c)/(2*a);
k = (2*o)/c;

e2 = (1/pi) * (2/(k*a)^3);
Ha = (besselh(1,1,k*a) - besselh(3,1,k*a))/2  ; %derivative of hankel function


A1 = (4i * rho * a  * e * o * o )/(k * Ha);
A1 = A1 * exp(-1i *(pi/2));

H = besselh(nu, 1 , k.*r);

pa = A1 .* H .* exp(2i.*(theta - o*t));

plot(r,real(pa));
hold on;



