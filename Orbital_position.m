%Orbit - Program to compute the orbit of a comet
clear all; help orbit; %Clear memory and print Header

%*Set initial position and velocity of the comet
%r0 = input ('Enter initial radial distance (AU): ');
%v0 = input ('Enter initial tangential velocity (AU/yr): ');
r0 = 2;
v0 = 1/2 * pi;

r = [r0 0]; v = [0 v0];
state = [ r(1) r(2) v(1) v(2)]; %Used by R-K routines

%*Set physical parameters (mass, G*M)
GM = 4*pi^2; %Grav const * Mass of sun (AU^3/yr^2)
mass = 1. ;  %Mass of comet
%adaptErr = 1.e-3; %Error parameter used by adaptive Runge-Katta
time = 0;

%* Loop over desired number of steps using specified numerical method
%nStep = input('Enter number of steps: ');
%tau = input('Enter time step (yr): ');
nStep = 10000;
tau = 0.001;

%NumericalMethod = menu('Choose a numerical method:', 'Euler', 'Euler-Cromer', 'Runge-Kutta', 'Adaptive R-K');
NumericalMethod = 3;
for iStep = 1:nStep
    %* Record position and energy for plotting
    rplot(iStep) = norm(r); % Record position for polar plot
    thplot(iStep) = atan2(r(2),r(1));
    tplot(iStep) = time;
    kinetic(iStep) = 0.5*mass*norm(v)^2; %Record energies
    potential(iStep) = -GM*mass/norm(r);

    %* Calculate new position and velocity using desired method
%    if (NumericalMethod == 1)
%        accel = -GM*r/norm(r)^3;
%        r = r + tau*v; %Euler Step
%        v = v + tau*accel;
%        time = time + tau;
%    elseif (NumericalMethod == 2)
%        accel = -GM*r/norm(r)^3;
%        v = v + tau*accel;
%        r = r + tau*v; %Euler-Cromer step
%        time = time + tau;
    if (NumericalMethod == 3)
        state = rk4(state, time, tau, 'gravrk', GM);
        r = [state(1) state(2)];   %4th order Runge-Kutta
        v = [state(3) state (4)];
        time = time + tau;
%    else
%        [state time tau] = rk4(state, time, tau, adaptErr, 'gravrk', GM);
%        r = [state(1) state(2)]; %Adaptive Runge-Kutta
%        v = [state(3) state(4)];
    end
end

%*Graph the trajectory of the comet
figure(1); clf; %Clear figure 1 window and bring forward 
polarplot(thplot, rplot, '+'); %Use polar plot for graphing orbit
pause(1) %Pause for 1 second before drawing next plot

%*Graph the energy of the comet versus time
figure(2); clf;
totalE = kinetic + potential;
plot(tplot, kinetic, '-', tplot, potential, '--', tplot, totalE, '-')
legend('Kinetic', 'Potential', 'Total');
xlabel('Time (yr)'); ylabel('Energy (M AU^2/yr^2)');
title('Energy v. Time');

function xout = rk4(x, t, tau, derivsRK, param)
half_tau = 0.5*tau;
F1 = feval(derivsRK, x, t, param);
t_half = t + half_tau;
xtemp = x + half_tau*F1;
F2 = feval(derivsRK, xtemp, t_half, param);
xtemp = x + half_tau*F2;
F3 = feval(derivsRK, xtemp, t_half, param);
t_full = t + tau;
xtemp = x + tau*F3;
F4 = feval(derivsRK, xtemp, t_full, param);
xout = x + tau/6 .* (F1 + F4 + 2 .* (F2 + F3));
end

function deriv = gravrk(s, t, GM)
r = [s(1) s(2)];
v = [s(3) s(4)];
accel = -GM*r/norm(r)^3;

deriv = [v(1) v(2) accel(1) accel(2)];
end



