%%Electron Bonding in Molecular Orbitals
clear all;

%Parameters
h_bar = (6.62607015e-34) / (2*pi);              %J*s

m = 9.1093837e-31;                              %Electron Mass: kg
t = 1;                                          %Hopping/Resonance integral
E0 = 0;                                         %Initial binding energy to center
N = 2;                                          %Number of electrons in chain
max = 100;

%Begin by Constructing the Hamiltonian
%J. G. Analytis notes ham = -(h_bar^2 / 2 * m)(Laplacian) + Va + Vb
%Electron bound in a linear chain.
   
for n=2:max                                 %Loop the program to show multiple N's
    E = 1:n;                                %Create the array E for the eenergy
    ham = zeros(n);                         %Create the empty hamiltonian matrix 
    for i=2:(n)                             %Populate the Hamitonian
        ham(i, i-1) = t;
        ham(i, i) = E0;
        ham(i-1, i) = t;
    end
    for i=1:n                               %Calculate values of E for each N
        E(i) = 2 * abs(t) * cos((pi * i) / (n + 1));
    end
    E = diag(E);
    ham = ham - E;
    [evec, eval] = eig(ham);
    g = diag(eval);
    g = g.';
    figure(1);
    plot(n, g, '_');
    xlabel("N"); ylabel("Orbital Levels"); title("Chain");
    hold on;
end
hold off;
%Electrons bound in a ploygon shape.
for p=2:max

    E2 = 1:p;
    ham2 = zeros(p);
    for i=2:(p)
        ham2(i, i-1) = t;
        ham2(i, i) = E0;
        ham2(i-1, i) = t;
    end
    ham2(1, p) = t;
    ham2(p, 1) = t;
    for i=1:p
        E2(i) = 2 * abs(t) * cos((2 * pi * i) / (p));
    end
    E2 = diag(E2);
    ham2 = ham2 - E2;
    [evec2, eval2] = eig(ham2);
    g2 = diag(eval2);
    g2 = g2.';
    figure(2);
    plot(p, g2, '_');
    xlabel("N"); ylabel("Orbital Levels"); title("Polygon");
    hold on;
end
hold off;

