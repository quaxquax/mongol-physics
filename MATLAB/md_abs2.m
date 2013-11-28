% md_abs2.m
    function y = md_abs2(x)
        y = 0;
        for i=1:length(x)
            y = y + abs(x{i}).^2;
        end
    end