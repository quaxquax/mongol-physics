% md_step.m

% Given solution of MD system at time t,
% compute solution at time t + dt.
% Vex and Aex are given external potentials at time t + dt. 
%
% V and dV are single arrays
% A and dA are 3-item cell arrays for A1, A2, A3
% psi is a 4-item cell array for psi1, psi2, psi3, psi4 
%
% space and fspace are coordinates, as computed by md_domain
% epsilon and delta are scale parameters for asymptotic solution

function [psi_new, V_new, dV_new, A_new, dA_new] = md_step(psi_n , V_n , dV_n , A_n , dA_n , space, fspace, dt, epsilon, delta, Vex, Aex)

    % Fourier transform of psi
    for i=1:4
        psihat_n{i} = fftn(psi_n{i});
    end

    % Forward exponential matrix
    M = md_expmatrix1(epsilon, delta, fspace, dt);
    for i=1:4
        phihat_new{i} = M{i,1}  .*psihat_n{1} + M{i,2}  .*psihat_n{2} + M{i,3}  .*psihat_n{3} + M{i,4}  .*psihat_n{4};
    end 

    % Inverse Fourier transform of phi_{n+1}
    for i=1:4
        phi_new{i} = ifftn(phihat_new{i});
    end

    % Compute current and particle densities
    % 2.21
    % rho_n = |psi_n|^2
    % rho_{n+1} = |phi_{n+1}|^2
    rho_n = md_abs2(psi_n);
    rho_new = md_abs2(phi_new);
    J_n = md_diracprod(delta, psi_n); J_new = md_diracprod(delta, phi_new);

    % FFT
    Vhat_n = fftn(V_n);
    dVhat_n = fftn(dV_n);
    rhohat_n = fftn(rho_n);
    rhohat_new = fftn(rho_new);
    for i=1:3
        Ahat_n{i} = fftn(A_n{i});
        dAhat_n{i} = fftn(dA_n{i});
        %Jhat_n{i} = fftn(J_n{i});
        %Jhat_new{i} = fftn(J_new{i});
        Jhat_sum{i} = fftn(J_n{i} + J_new{i});
    end

    % Update A and V with the Crank-Nicolson method
    % 2.19-2.21

    % Crank-Nicolson update matrix, (2.19)
    % This could be precomputed for fixed dt, but it?s a small
    % part of the total computation time anyway.
    % |xi|^2
    absxi2 = md_abs2(fspace);
    CN_factor = 1 ./ (1 + (dt^2 * absxi2 / (4*delta^2))); CN_M11 = 1 - (dt^2 * absxi2 / (4*delta^2));
    CN_M12 = dt;
    CN_M21 = -dt*absxi2/delta^2;
    CN_M22 = 1 - (dt^2 * absxi2 / (4*delta^2));
    CN_B1 = epsilon * dt^2 / (4*delta^2);
    CN_B2 = epsilon * dt / (2*delta^2);
    % End CN update matrix

    BV = rhohat_n + rhohat_new;
    Vhat_new     = (CN_M11 .*Vhat_n+CN_M12 .*dVhat_n+CN_B1 .*BV) .*CN_factor;
    dVhat_new   = (CN_M21 .*Vhat_n+CN_M22 .*dVhat_n+CN_B2 .*BV) .*CN_factor;
    for i=1:3
        %BA =delta * (Jhat_n{i} + Jhat_new{i});
        BA = delta * Jhat_sum{i};
        Ahat_new{i}     = (CN_M11 .*Ahat_n{i} + CN_M12 .*dAhat_n{i} + CN_B1 .*BA)  .* CN_factor;
        dAhat_new{i}  = (CN_M21 .*Ahat_n{i} + CN_M22 .*dAhat_n{i} + CN_B2 .*BA)  .* CN_factor;
    end

    % Updated V and A and their derivatives
    V_new = ifftn(Vhat_new);
    dV_new = ifftn(dVhat_new);
    for i=1:3
        A_new{i} = ifftn(Ahat_new{i});
        dA_new{i} = ifftn(dAhat_new{i});
    end

    M = md_expmatrix2(epsilon, delta, V_new, Vex, A_new, Aex, dt);
    for i=1:4
        psi_new{i} = M{i,1} .*phi_new{1} + M{i,2} .*phi_new{2} + M{i,3} .*phi_new{3} + M{i,4} .*phi_new{4};
    end

    % Possibly remove imaginary noise
    V_new = real(V_n);
    dV_new = real(dV_n);
    for i=1:3
        A_new{i} = real(A_new{i});
        dA_new{i} = real(dA_new{i});
    end
end


