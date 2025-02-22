%ball - compute the trajectory of a baseball including air resistance using
%euler's method
clear; help ball; %Clear memory and print header
% Set intial position and velocity of baseball
y1 = input ('Enter intial height (meters): ');
r1 = [0, y1]; %initial vector position
speed = input ('Enter initial speed (m/s): ');
theta = input ('Enter initial angle (degrees): ');
v1 = [speed*cos(theta*pi/180), speed*sin(theta*pi/180)];
r = r1; v = v1; %Set initial position and velocity

%* Set physical parameters (mass, Cd, etc. )
Cd = 0.35; %Drag coefficient (dimensionless)
area = 0.6567; %Cross-sectional area of projectile (m^2)
grav = 9.81; %Gravitational acceleration (m/s^2)
mass = 45.3; %Mass of projectile (kg)
airFlag = input ('Air resistance? (Yes: 1, No: 0): ');
if (airFlag == 0)
    rho = 0; %No air resistance
else
    rho = 1.2; %Density of air (kg/m^3)
end
air_const = -0.5*Cd*rho*area/mass; %Air resistance constant

%* Loop until ball hits the ground or max steps completed
tau = input ('Enter time step, tau (sec): '); %(sec)
maxstep = 25; %maximum number of steps
for istep = 1:maxstep


%* Record position (computed and theoretical) for plotting
xplot(istep) = r(1); %Record trajectory for plot
yplot(istep) = r(2);
t = (istep-1)*tau; %Current time
xNoAir(istep) = r1(1) + v1(1)*t;
yNoAir(istep) = r1(2) + v1(2)*t - 0.5*grav*t.^2;

%* Calculate acceleration of ball
accel = air_const*norm(v)*v; %Air resistance
accel(2) = accel(2)-grav; %Gravity

%* Calculate the new position and velocity using Euler's method
r = r + tau*v; %Euler step
v = v + tau*accel;

%* If ball reaches ground (y<0), break out of loop
if (r(2) < 0)
    xplot(istep+1) = r(1); %Record last values computed
    yplot(istep+1) = r(2);
end
end


%*Print maximum range and time of flight
fprintf('Maximum range is %g meters\n', r(1));
fprintf('Time of flight %g seconds\n',istep*tau);

%* Graph the trajectory of the baseball
clf; figure(gcf); %Clear windor and bring it forward
%mark the location of the ground by a straight line
xground = [0, max(xNoAir)]; yground = [0, 0];
% Plot the computed trajectory and parabolic, no-air curve
plot(xplot, yplot, '+', xNoAir, yNoAir, '-', xground, yground, '-');
legend('Euler method', 'Theory (No air)');
xlabel('Range (m)'); ylabel('Height (m)');
title('Projectile Motion');
