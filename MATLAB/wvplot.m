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