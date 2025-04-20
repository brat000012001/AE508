%
% Computes the Hamiltonian and adjoins an additional path cost term.
%
function Hamiltonian = hamiltonian_minFuel(t,X,T,c,rho,mu)
    addpath("..");
    Hamiltonian = hamiltonian(t,X,T,c,rho,mu);
    for idx = 1:length(Hamiltonian)
        Hamiltonian(idx) = Hamiltonian(idx) + T/c;
    end
end