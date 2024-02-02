clc 
clear 

x = linspace(10,300,100);
theta = linspace(0,0,100);
y = linspace(0,0,100);
timestep = 1.515*1e-2;
sum = 0;

for i = 1:1e2
        t = timestep*i;
        p = pressurefunc(x,y,t) ;
        pa = analyticalpressure(x,t);
        pressure = real(p) + real(pa);
        if i==1
            pres0= pressure;
        end
        if i==2
            presPrev = pressure;
        end
        if i==3 
            sum = timestep/2 .* (pressure + pres0 + presPrev);
        end
        if i > 3 
            sum = sum + timestep/2 .* (pressure + presPrev);
            presPrev = pressure;
            presAv = sum./t;
            dp = pressure - presAv;

            %plotting 
            figure(2)
            plot(x,real(dp)./(330^2));
            %ylim ([-1e-7 1e-7])
            ylim tight
            hold on
        end
end

save 1.mat presAv

