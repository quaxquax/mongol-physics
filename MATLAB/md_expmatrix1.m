% md_expmatrix1.m
% Note: formulas refer to the paper by Huang et al.

    function M = md_expmatrix1(epsilon, delta, xi, dt)
        xi1 = xi{1};
        xi2 = xi{2};
        xi3 = xi{3};
        % |xi|^2
        absxi2 = xi1 .*xi1 + xi2 .*xi2 + xi3 .*xi3;
        % |epsilon * delta * xi|^2
        absepsdeltaxi2 = absxi2 * epsilon^2 * delta^2;
        % 2.16
        lambda = 1i/(epsilon*delta^2)*sqrt(1+epsilon^2*delta^2*absxi2);
        % 2.18
        cl = cos(-1i*lambda*dt);
        sl = sin(-1i*lambda*dt)  .* ((1+absepsdeltaxi2) .^ (-0.5));
        M = cell(4,4);
        % 2.17
        M{1,1} = cl - 1i .*sl;
        M{1,2} = 0;
        M{1,3} = -1i  .* epsilon .* delta .* sl .* xi3;
        M{1,4} = -epsilon .* delta .* sl .* (xi2+1i .*xi1);
        M{2,1} = 0;
        M{2,2} = cl - 1i*sl;
        M{2,3} = epsilon .* delta .* sl .* (xi2-1i .*xi1);
        M{2,4} = 1i .* epsilon .* delta .* sl .* xi3;
        M{3,1} = -1i .* epsilon .* delta .* sl .* xi3;
        % note: typo in paper
        M{3,2} = -epsilon .* delta .* sl  .* (xi2+1i .*xi1);
        M{3,3} = cl + 1i .*sl;
        M{3,4} = 0;
        M{4,1} = epsilon .* delta .* sl .* (xi2-1i .*xi1);
        M{4,2} = 1i .* epsilon .* delta .* sl .* xi3; M{4,3} = 0;
        M{4,4} = cl + 1i .*sl;
    end
