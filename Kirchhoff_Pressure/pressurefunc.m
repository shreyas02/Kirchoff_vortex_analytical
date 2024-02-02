function pressure = pressurefunc(x,y,t)

% cartesian coordinates
xi = x;
yi = y;
a = 1.001;
b =  0.999; 
c = sqrt(a^2 - b^2);
 
%data given 

rho = 1;
cs = 330;
o = (0.1*cs)/2 ; 
omega = 2*o;

%Position in rotational frame
test = [cos(o*t),sin(o*t);-sin(o*t),cos(o*t)] * [x;y];
x = test(1,:);
y = test(2,:);
clear test t;

% eliptical coordinates conversion

t = (x+ 1i*y) / c;
t = acosh(t);
si =real(t); 
nu = imag(t);

%velocity components inertial frame t = 0

ui = zeros(1,size(x,2));
vi = zeros(1,size(x,2));

for(i = 1 : size(x,2))
    test = (x(i)^2)/a^2 + (y(i)^2)/b^2  -1 ;
    if(test <= 0)
        ui(i) = -omega *  (a/(a+b)) * y(i);
        vi(i) = omega *  (b/(a+b)) * x(i);
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

%plotting pressure 
% figure(3)
% plot(xi,pressure);
