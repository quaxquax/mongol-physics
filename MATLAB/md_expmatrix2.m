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
