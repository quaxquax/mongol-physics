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
