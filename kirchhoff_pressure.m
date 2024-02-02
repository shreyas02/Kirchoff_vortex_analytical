clc
clear 

% cartesian coordinates

x = linspace(0,5,100);
y = linspace(0,0,100);
xi = linspace(0,5,100);
yi = linspace(0,0,100);
a = 1;
b =  0.5; 
c = sqrt(a^2 - b^2);
 
%data given 
omega = 10;
o = 20/9 ; 
rho = 1.3;
cs = 330;
t = 0;

for n = 0 : 10000

t = 0.0001*n;
%Position in rotational frame
test = [cos(o*t),sin(o*t);-sin(o*t),cos(o*t)] * [x;y];
x = test(1,:);
y = test(2,:);
clear test;

% eliptical coordinates conversion

t = (xi+ 1i*yi) / c;
t = acosh(t);
si =real(t); 
nu = imag(t);

%velocity components inertial frame t = 0

ui = zeros(1,size(x,2));
vi = zeros(1,size(x,2));

for(i = 1 : size(x,2))
    test = (x(i)^2)/a^2 + (y(i)^2)/b^2  -1 ;
    if(test <= 0)
        ui(i) = -omega *  (a/(a+b)) * yi(i);
        vi(i) = omega *  (b/(a+b)) * xi(i);
    end
    if(test > 0)
        J = c.*c.*(sinh(si(i)).^2 + sin(nu(i)).^2);
        J = 1./J;
        temp = 0.5 * omega * a * b * c * J;
        ui(i) = temp * (cosh(si(i)) * sin(nu(i)) * (cos(2*nu(i))*exp(-2*si(i)) -1) + sinh(si(i))*cos(nu(i))*sin(2*nu(i))*exp(-2*si(i)));
        vi(i) = temp * (sinh(si(i)) * cos(nu(i)) * (-cos(2*nu(i))*exp(-2*si(i)) + 1) + cosh(si(i))*sin(nu(i))*sin(2*nu(i))*exp(-2*si(i)));
        clear temp;
    end
end
clear test;

%velocity components in rotational frame at t = 0 
ur = ui + o.*y;
vr = vi - o.*x;

%stream function rotational frame

psir = zeros(1,size(x,2));
for(i = 1 : size(x,2))
    test = (x(i)^2)/a^2 + (y(i)^2)/b^2  -1 ;
    if(test <= 0)
        psir(i) = 0.5*omega.*((b/(a+b))^2.*x(i)^2 + (a/(a+b))^2.*y(i)^2) - 0.5*o*a*b;
    end
end
clear test;

%pressure calculation 
pressure = zeros(1,size(x,2));
for(i = 1 : size(x,2))
    test = (x(i)^2)/a^2 + (y(i)^2)/b^2  -1 ;
    if(test <= 0)
        pressure(i) = rho * (omega * psir(i)  +  0.5 * o *o* (x(i)^2 +y(i)^2) - 0.5*(ur(i)^2 + vr(i)^2));
    end
    if(test > 0)
    pressure(i) = rho*(0.5 * o *o* (x(i)^2 +y(i)^2) - 0.5*(ur(i)^2 + vr(i)^2));
    end
end

%plotting relative vel 
% figure(1)
% plot(x,vr);
% hold on;

%plotting paper data relative vel
% n = readtable('vrel.csv');
% xg = table2array(n(:,1));
% yg = table2array(n(:,2));
% figure(1)
% scatter(xg,yg);

%plotting inertial vel
% figure(2)
% plot(x,vi);
% hold on;

%plotting paper data inertial vel
% n = readtable('vin.csv');
% xg = table2array(n(:,1));
% yg = table2array(n(:,2));
% figure(2)
% scatter(xg,yg);

%plotting pressure 
figure(3)
plot(xi,pressure);
% hold on;

%plotting paper data pressure
% n = readtable('pressure.csv');
% xg = table2array(n(:,1));
% yg = table2array(n(:,2));
% figure(3)
% scatter(xg,yg);
% hold off

end


