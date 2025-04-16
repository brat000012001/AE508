%
% AE508 Optimal Space Trajectories, Spring 2025
% Course Project
% 
% Team: petern4@illinois.edu, pc46@illinois.edu, davisr2@illinois.edu
%

%
% Calculates the Hamiltonian
%
function Hamiltonian = hamiltonian(t,X,T,c,rho,mu)
    Hamiltonian = zeros(length(t),1);
    rvec = X(:,1:3);
    vvec = X(:,4:6);
    mvec = X(:,7);
    rmag = vecnorm(rvec,2,2);
    lam_r = X(:,8:10);
    lam_v = X(:,11:13);
    lam_m = X(:, 14);
    uhat = -lam_v./vecnorm(lam_v,2,2);
    for idx = 1:length(t)
        Hamiltonian(idx) = lam_r(idx,:)*vvec(idx,:)' ...
            - lam_v(idx,:)*rvec(idx,:)'*mu/rmag(idx)^3 ...
            + lam_v(idx,:)*uhat(idx,:)'*T/mvec(idx) - lam_m(idx)*T/c; 
    end
end
