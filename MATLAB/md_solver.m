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
        phihat_new{i} = M{i,1}.*psihat_n{1} + M{i,2}.*psihat_n{2} + M{i,3}.*psihat_n{3} + M{i,4}.*psihat_n{4};
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
    Vhat_new     = (CN_M11.*Vhat_n+CN_M12.*dVhat_n+CN_B1.*BV).*CN_factor;
    dVhat_new   = (CN_M21.*Vhat_n+CN_M22.*dVhat_n+CN_B2.*BV).*CN_factor;
    for i=1:3
        %BA =delta * (Jhat_n{i} + Jhat_new{i});
        BA = delta * Jhat_sum{i};
        Ahat_new{i}     = (CN_M11.*Ahat_n{i} + CN_M12.*dAhat_n{i} + CN_B1.*BA) .* CN_factor;
        dAhat_new{i}  = (CN_M21.*Ahat_n{i} + CN_M22.*dAhat_n{i} + CN_B2.*BA) .* CN_factor;
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
        psi_new{i} = M{i,1}.*phi_new{1} + M{i,2}.*phi_new{2} + M{i,3}.*phi_new{3} + M{i,4}.*phi_new{4};
    end

    % Possibly remove imaginary noise
    V_new = real(V_n);
    dV_new = real(dV_n);
    for i=1:3
        A_new{i} = real(A_new{i});
        dA_new{i} = real(dA_new{i});
    end
end

% sp_step.m
    function [phi_e,phi_p] = sp_step(phi_e, phi_p, Vex, dt, space, fspace)
        
        x=space{1};
        y=space{2};
        z=space{3};
        px=fspace{1};
        py=fspace{2};
        pz=fspace{3};
        xi2 = px.^2 + py.^2 + pz.^2;
        
        % Solve i * (d/dt) phi_{e/p} = (-/+ Delta/2) phi_{e/p}
        phi_e = ifftn( exp(-0.5i*xi2*dt) .* fftn(phi_e) );
        phi_p = ifftn( exp(0.5i*xi2*dt) .* fftn(phi_p) );
        
        % Poisson?s equation
        % -Delta V = |phi_p|^2 + |phi_e|^2
        Vhat = fftn( abs(phi_e).^2 + abs(phi_p).^2 ) ./ xi2;
        Vhat(1,1,1) = 0;
        V = ifftn(Vhat);
        
        % Solve i * (d/dt) phi_{e/p} = (V + Vex) phi_{e/p}
        phi_e = exp(-1i*dt*(V+Vex)) .* phi_e;
        phi_p = exp(-1i*dt*(V+Vex)) .* phi_p;
    end

% md_expmatrix1.m
% Note: formulas refer to the paper by Huang et al.

    function M = md_expmatrix1(epsilon, delta, xi, dt)
        xi1 = xi{1};
        xi2 = xi{2};
        xi3 = xi{3};
        % |xi|^2
        absxi2 = xi1.*xi1 + xi2.*xi2 + xi3.*xi3;
        % |epsilon * delta * xi|^2
        absepsdeltaxi2 = absxi2 * epsilon^2 * delta^2;
        % 2.16
        lambda = 1i/(epsilon*delta^2)*sqrt(1+epsilon^2*delta^2*absxi2);
        % 2.18
        cl = cos(-1i*lambda*dt);
        sl = sin(-1i*lambda*dt) .* ((1+absepsdeltaxi2) .^ (-0.5));
        M = cell(4,4);
        % 2.17
        M{1,1} = cl - 1i.*sl;
        M{1,2} = 0;
        M{1,3} = -1i .* epsilon .* delta .* sl .* xi3;
        M{1,4} = -epsilon .* delta .* sl .* (xi2+1i.*xi1);
        M{2,1} = 0;
        M{2,2} = cl - 1i*sl;
        M{2,3} = epsilon .* delta .* sl .* (xi2-1i.*xi1);
        M{2,4} = 1i .* epsilon .* delta .* sl .* xi3;
        M{3,1} = -1i .* epsilon .* delta .* sl .* xi3;
        % note: typo in paper
        M{3,2} = -epsilon .* delta .* sl .* (xi2+1i.*xi1);
        M{3,3} = cl + 1i.*sl;
        M{3,4} = 0;
        M{4,1} = epsilon .* delta .* sl .* (xi2-1i.*xi1);
        M{4,2} = 1i .* epsilon .* delta .* sl .* xi3; M{4,3} = 0;
        M{4,4} = cl + 1i.*sl;
    end

% md_expmatrix2.m

% Compute the second exponential matrix
    function M = md_expmatrix2(epsilon, delta, V, Vex, A, Aex, dt)
        V = V + Vex;
        A1 = A{1} + Aex{1};
        A2 = A{2} + Aex{2};
        A3 = A{3} + Aex{3};
        M = cell(4,4);
        Aabs = sqrt(abs(A1).^2 + abs(A2).^2 + abs(A3).^2);
        %RA = 1 ./ Aabs;
        RA = A1 * 0;
        % Note: this is a normalization constant. Where 1/|A| is undefined, % the corresponding component in the matrix will be zero.
        for k=1:numel(Aabs)
            if Aabs(k) == 0
                RA(k) = 0;
            else
                RA(k) = 1./Aabs(k);
            end
        end
        
        B1 = A1 .* RA;
        B2 = A2 .* RA;
        B3 = A3 .* RA;
        u = Aabs .* dt ./ epsilon;
        v = dt .* V ./ epsilon;
        su = sin(u);
        cu = cos(u);
        sv = sin(v);
        cv = cos(v);
        a = cv-1i .*sv;
        b = 1i .*cv + sv;
        M{1,1}=cu .*a; M{1,2}=0;
        M{1,3}=B3 .*su .*b; M{1,4}=(1i .*B1+B2) .*su .*a;
        M{2,1}=0; M{2,2}=M{1,1};
        M{2,3}=(B1+1i .*B2) .*su .*b; M{2,4}=-1i .*B3 .*su .*a;
        M{3,1}=M{1,3}; M{3,2}=M{1,4}; M{3,3}=M{1,1}; M{3,4}=0;
        M{4,1}=M{2,3}; M{4,2}=M{2,4}; M{4,3}=0; M{4,4}=M{1,1};
    end

% md_domain.m

% space--3-item cell array of space cooordinates [x,y,z]
% fspace--3-item cell array of Fourier coordinates [xi1,xi2,xi3] 
% each coordinate array is N x N x N

        function [space, fspace] = md_domain(xd, yd, zd, N)
            
            pos_x1 = xd(1); pos_x2 = xd(2);
            pos_y1 = yd(1); pos_y2 = yd(2);
            pos_z1 = zd(1); pos_z2 = zd(2);
            xwidth = pos_x2 - pos_x1;
            ywidth = pos_y2 - pos_y1;
            zwidth = pos_z2 - pos_z1;
            %volume = xwidth * ywidth * zwidth;
            xs = linspace(pos_x1, pos_x2, N+1);
            ys = linspace(pos_y1, pos_y2, N+1);
            zs = linspace(pos_z1, pos_z2, N+1);
            [xx, yy, zz] = ndgrid(xs, ys, zs);
            xs = xs(1:N); ys = ys(1:N); zs = zs(1:N);
            
            % Coordinates xi in Fourier (momentum) space
            % Note that the Fourier modes range from -N/2+1 to N/2, but
            % the FFT indexing starts at 0, so there is wraparound.
            pv = [];
            for p=0:N-1
                if p < floor(N/2) pv(p+1) = p;
                else
                    pv(p+1) = p - N;
                end
            end
            
            [xi1,xi2,xi3] = ndgrid(2*pi*pv/xwidth,2*pi*pv/ywidth, 2*pi*pv/zwidth);
            
            space = cell(3,1);
            fspace = cell(3,1);
            
            space{1} = xx;
            space{2} = yy;
            space{3} = zz;
            
            fspace{1} = xi1;
            fspace{2} = xi2;
            fspace{3} = xi3;
        end

% md_abs2.m
    function y = md_abs2(x)
        y = 0;
        for i=1:length(x)
            y = y + abs(x{i}).^2;
        end
    end
    
% md_diracprod.m

% Computes {Y_1, Y_2, Y_3} where Y_k = <X, alpha^k X> / c 
% where <x,y> is the C^4 inner product with x conjugated 
% alpha^k is a Dirac matrix

    function Y = md_diracprod(c, X)
        Y = cell(3,1);
        A1 = conj(X{1});
        A2 = conj(X{2});
        A3 = conj(X{3});
        A4 = conj(X{4});
        Y{1} = A1 .*X{4} + A2 .*X{3} + A3 .*X{2} + A4 .*X{1};
        Y{2} = i*(A2 .*X{3} + A4 .*X{1} - A1 .*X{4} - A3 .*X{2});
        Y{3} = A1 .*X{3} + A3 .*X{1} - A2 .*X{4} - A4 .*X{2};
        Y{1} = Y{1} / c;
        Y{2} = Y{2} / c;
        Y{3} = Y{3} / c;
    end
    
% flatten2.m
    function y = flatten2(f, N)
        abspsi = md_abs2(f); y = abspsi(:,:,1);
        for i=2:N
            y = y + abspsi(:,:,i);
        end
        y = y / N;
    end
    
% wvplot.m
    function data = wvplot(x, y, N, f, fnorm, srows, scols, sn, msg)
        subplot(srows,scols,sn);
        Z = flatten2(f, N) ./ fnorm;
        surf(x(:,:,N/2), y(:,:,N/2), Z);
        zlim([0,4]);
        view(-20,30);
        pbaspect([1 1 1]);
        title(msg);
        data = Z;
    end
