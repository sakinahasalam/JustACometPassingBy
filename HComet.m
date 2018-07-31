%% Halley's comet


% // An example to study the time-dependent position and
% // velocity of Halley's comet via the Verlet algorithm

clear all;

%%constant Parameters
global kappa;
global n;
global m;

m =20000;
n = 200000;
kappa = 39.478428;

h = 2.0/(n-1);
h2 = h*h/2;

%% Intial condition

% Two quantities can be assumed as the known
% quantities, the total energy and the angular momentum, as constants
% describe the motion of the comet in the xy plane and choose
% x0 = rmax, vx0 = 0, y0 = 0, and vy0 = vmin

t(1) = 0; %initial time
x(1) = 1.966843; % initial x position of comet
y(1) = 0; %initial y posiyion of
r(1) = x(1);
vx(1) = 0;
vy(1) = 0.851795;
gx(1) = -kappa./(r(1)*r(1));
gy(1) = 0;

Count = m;

%% Function
% Using Verlet Algorithm Method 
% - Formula can be referred from Intro of Comp Physics
% Verlet algorithm yields very high accuracy for the position
% BUT - velocity is evaluated using one time step behind of the position.
for n=1:m
    t(n+1) = h*(n+1);
    x(n+1) = x(n) + h.*vx(n) + h2.*gx(n);
    y(n+1) = y(n) +  h.*vy(n) + h2.*gy(n);
    r2 = x(n+1).*x(n+1) + y(n+1).*y(n+1);
    r(n+1) = sqrt(r2);
    r3 = r2.*r(n+1);
    gx(n+1) = -kappa.*x(n+1)./r3;
    gy(n+1) = -kappa.*y(n+1)./r3;
    vx(n+1) = vx(n) + h.*(gx(n+1)+gx(n))/2;
    vy(n+1) = vy(n) + h.*(gy(n+1)+gy(n))/2; 
    
end

%% Plotting
figure

plot(t,r)
title('Distance between Halley Comet and Sun')
axis tight;
ylabel('r');
xlabel('t');