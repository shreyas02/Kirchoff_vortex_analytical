clc 
clear 

%% r = 10
rho = 1.3;
a = 2;
nu = 2 ;
c = 330;
e = 1.25e-3;
mach = 0.076;
t = 0 ;
theta = linspace(0, 2*pi,100);
r = 10;
o = 2*pi;
k = (2*o)/c;

e2 = (1/pi) * (2/(k*a)^3);
Ha = (besselh(1,1,k*a) - besselh(3,1,k*a))/2  ; %derivative of hankel function


A1 = (4i * rho * a  * e * o * o )/(k * Ha);
A1 = A1 * exp(-1i *(pi/2));

H = besselh(nu, 1 , k*r);

pa = A1 * H .* exp(2i.*(theta - o*t));

plot(theta,real(pa));
hold on;

%% r = 2
rho = 1.3;
a = 2;
nu = 2 ;
c = 330;
e = 1.25e-3;
mach = 0.076;
t = 0 ;
theta = linspace(0, 2*pi,100);
r = 2;
o = 2*pi;
k = (2*o)/c;

e2 = (1/pi) * (2/(k*a)^3);
Ha = (besselh(1,1,k*a) - besselh(3,1,k*a))/2  ; %derivative of hankel function


A1 = (4i * rho * a  * e * o * o )/(k * Ha);
A1 = A1 * exp(-1i *(pi/2));

H = besselh(nu, 1 , k*r);

pa = A1 * H .* exp(2i.*(theta - o*t));

plot(theta,real(pa));
hold on;

%% plotting paper data 
n = readtable('r.csv');
x = table2array(n(:,1));
y = table2array(n(:,2));
scatter(x,y);


