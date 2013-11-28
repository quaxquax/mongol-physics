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
            xs = linspace(pos_x1, pos_x2, N);
            ys = linspace(pos_y1, pos_y2, N);
            zs = linspace(pos_z1, pos_z2, N);
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
