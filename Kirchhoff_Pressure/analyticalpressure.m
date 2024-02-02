function pa = analyticalpressure(x,t)

rho = 1;
a = 1;
nu = 2 ;
c = 330;
e = 0.001;
mach = 0.1;
theta = 0;
r = x;
o = (0.1*c)/2 ; 
k = (2*o)/c;

e2 = (1/pi) * (2/(k*a)^3);
Ha = (besselh(1,1,k*a) - besselh(3,1,k*a))/2  ; %derivative of hankel function


A1 = (4i * rho * a  * e * o * o )/(k * Ha);
A1 = A1 * exp(-1i *(pi/2));

H = besselh(nu, 1 , k.*r);

% for i = 1:size(r,2)
%         for j = 1:size(theta,2)
%             pa(i,j) = A1 * H(i) * exp(2i*(theta(j) - o*t));
%         end
% end
pa = A1 .* H.* exp(2i*(theta - o*t));

% plotting graph pa
% figure(1)
% plot(r,real(pa)./(330^2));
% ylim([-0.3e-7 0.3e-7]);



