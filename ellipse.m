clc 
clear

a = 2.0025; 
b=1.9975;
nu = 0 ;

x = linspace(-2.0025,2.0025);

y = sqrt( b^2.*(1 - ((x.^2)./(a^2))) );
y = [y , -y];
x = [x, x];

% figure(1);
% scatter(x,y);
% xlim([-10 10]);
% ylim([-10 10]);

%declare variables
c = 330;
e = 1.25e-3;
rho = 1.3;
mach = 0.076;
timestep = 0.01;
% x2 = linspace(0.1,10,1000);
x2 = 10;
theta = linspace(0.1, 2*pi,100);
o = (c*mach)/(2*a);
k = (2*o)/c;
H = besselh(1,2,k*x2);
Ha = (besselh(0,2,k*a) - besselh(2,2,k*a))/2  ;

A = (4i * rho * a * e * o^2 * exp(-1i*(pi/2)) )/(k*Ha);

%time loop starts here
for n =1:100
        
        t = n * timestep;
        % rotation of matrix 
        rot = [cos(o*t) , sin(o*t) ; - sin(o*t), cos(o*t)] ;
        k1 = rot* [ x ; y];
        x1 = k1(1,:);
        y1 = k1(2,:);

        figure(2);
        title(t)
        scatter(x1,y1);
        xlim([-10 10]);
        ylim([-10 10]);
        axis equal;

        %calculation of pa along x axis 
        pa = A * H .* exp(2i.*(theta - o*t));
   
        figure(3);
        plot(theta,pa);
        xlim([0 2*pi]);

end


