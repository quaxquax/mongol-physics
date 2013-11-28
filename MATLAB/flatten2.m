% flatten2.m
    function y = flatten2(f, N)
        abspsi = md_abs2(f); y = abspsi(:,:,1);
        for i=2:N
            y = y + abspsi(:,:,i);
        end
        y = y / N;
    end
