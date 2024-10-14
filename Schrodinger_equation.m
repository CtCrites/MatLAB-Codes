%Schrodinger Equation

clear all; help schro;

%Menu
NUMPOT = menu('Select the Potential Function Form', 'None', ...
    'Quantum Harmonic Oscillator', 'cos^2 Form');

%Constants
i_imag = sqrt(-1);          %Imaginary number i
h_bar = 1;                  %Modified Planck's Constant
mass = 1;                   %Mass
x0 = 0;                     %Sets initial position
velocity = 1;               %Average velocity of the packet
k0 = mass*velocity/h_bar;
tau = 0.025;                    %Size of timestep
omega = 0.1;
lambda = pi/2;

%Grid Stuff
N = 200;                    %Number of grid points
L = 100;                    %System extends from -L/2 to L/2
h = L/(N-1);                %Grid size
x = h * (0:N-1) - L/2;      %Coordinates of grid points

%Hamiltonian Operator Matrix
ham = zeros(N);             %Create an N x N matrix full of 0's
coeff = -(h_bar^2) / (2 * mass * h^2);

for i=2:(N-1)
    ham(i, i - 1) = coeff;  %Sets edge of diagonal
    ham(i, i) = -2 * coeff; %Sets diagonal center
    ham(i, i + 1) = coeff;  %Sets bottom edge of diagonal
end

%First Rows for boundary conditions
ham(1, N) = coeff; 
ham(1, 1) = -2 * coeff;
ham(1, 2) = coeff;

%Last Rows for boundary conditions
ham(N, N - 1) = coeff;
ham(N, N) = -2 * coeff;
ham(N, 1) = coeff;

%Potential Function Stuff
switch NUMPOT
    case 1
        U = 0;
    case 2
        U = 0.5*omega^2*(x'-x0).^2;
    case 3
        U = cos(lambda*omega*(x'-x0)) .^2;
end

potential = diag(U);
ham = ham + potential;
[vec, eval] = eig(ham);
eval = diag(eval)/omega;



%Compute Crank-Nicholson Matrix
dCN = (inv(eye(N) + 0.5 * i_imag * tau / h_bar * ham) * ((eye(N)) - 0.5 * i_imag * tau / h_bar * ham));

%Initializes Wavefunction
sigma0 = L/10;
Norm_psi = 1 / (sqrt(sigma0 * sqrt(pi)));
psi = Norm_psi * exp(i_imag * k0 * x') .* exp(-(x' - x0) .^2 / (2 * sigma0^2));


%Plot the initial wavefunction
figure(1); clf;
plot(x, real(psi), x, imag(psi), '--')
xlabel('X'); ylabel('\psi(x)'); legend('Real', 'Imaginary');
drawnow; pause(1);

%Initialize loop for wave movement
max_iter = 2000;   %L / (velocity * tau);    %Maximum number of iterations
plot_iter = max_iter / 20;         
p_plot(:, 1) = psi .* conj(psi);
iplot = 1;
axisV = [-L/2 L/2 0 max(p_plot)];

%Loop over the desired number of steps
for iter=1:max_iter
    %Compute new wave function using Crank-Nicholson
    psi = dCN * psi;

    %Periodically record new values to plot
    if (rem(iter, plot_iter) < 1)
        iplot = iplot + 1;
        p_plot(:, iplot) = psi .* conj(psi);
        figure(2);
        plot(x, p_plot(:, iplot));
        xlabel('X'); ylabel('P(x,t)');
        title(sprintf('Finished %g of %g iterations', iter, max_iter));
        axis(axisV); %drawnow;
        pause(0.05);

    end
end

%Plot probability density versus position
figure(2);
pFinal = psi .* conj(psi);
plot(x, p_plot(:, 1:3:iplot), x, pFinal);
xlabel('X'); ylabel('P(x,t)');
title('Probability Density at Various Times');

%
figure(3);
plot(x, U);
xlabel('X'); ylabel('U(x)');
title('Potential Energy Plot');












