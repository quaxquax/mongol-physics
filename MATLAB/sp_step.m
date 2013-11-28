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
