classdef ZpinchArray
%     This class is used to visualize the magnetic field and
%     j x B force on unit j ez for Z pinch cylindical arrays in the 2D x-y plane
%     The current flows in the z-direction
    properties
        I % current, [A]
        N % Number of wires
        R0 % array radius, [m]
        d0 % wire diameter, [m]
        L % Box size
        isExploder % true for exploder arrays
        M % Number of points 
        pos0 % position of wires
    end
    methods
        function obj = ZpinchArray(N,R0,d0,I,isExploder)
            % Creates a Z pinch array 
            % Inputs 
            % N = number of wires
            % R0 = Array radius, [m]
            % d0 = wire diamter, [m]
            % I = current, [A]
            % isExploder = true for inverse wire arrays
            obj.N = N;
            obj.I = I;
            obj.d0 = d0;
            obj.R0 = R0;
            obj.isExploder = isExploder;
            
            % (2) Create wires at positions
            th0 = 0:2*pi/N:2*pi-2*pi/N; 
            x0 = R0 * cos(th0); % wire x-position
            y0 = R0 * sin(th0); % wire z-position
            obj.pos0 = [x0', y0']; % centers of wires

            % (2a) For exploder
            if isExploder
                obj.pos0 = [obj.pos0; 0, 0]; % cathode
            end
        end
        function obj = SquareArray(obj,N,R0,d0,I,isExploder)
            % Creates a square Z pinch array 
            % Inputs 
            % N = number of wires
            % R0 = length of side, [m]
            % d0 = wire diamter, [m]
            % I = current, [A]
            % isExploder = true for inverse wire arrays
                obj.N = N;
                obj.I = I;
                obj.d0 = d0;
                obj.R0 = R0;
                obj.isExploder = isExploder;

                % (2) Create wires at positions - square
                delS = obj.R0 * 4 / obj.N; 
                % R
                yr = -R0:delS:R0; xr = R0 * ones(size(yr)); 
                % L 
                yl = yr; xl =  -1 *  xr;  
                % T 
                yt = xr; xt = yr;
                % B 
                yb = -1 * yt; xb = xt;
                
                obj.pos0 = [[xr,xl,xt,xb]', [yr,yl,yt,yb]']; % centers of wires
                obj.pos0 = unique(obj.pos0,'rows');

                % (2a) For exploder
                if isExploder
                    obj.pos0 = [obj.pos0; 0, 0]; % cathode
                end
            
        end
        function [Bx_t,By_t,B] = getMagField(obj,pts_x, pts_y)
            % Returns the total magnetic field from teh superposition of
            % magentic field around each wire
            % Inputs:
            % pts_x = x-position of grid points
            % pts_y = y-position of grid points
            Bx_t = 0; By_t = 0; 
            for ii = 1:obj.N
                c0 = obj.pos0(ii,:); % current wire
                % Convert to "local" r-theta coordinates
                r = sqrt((pts_x-c0(1)).^2 + (pts_y-c0(2)).^2); 
                th = obj.getth(pts_x-c0(1),pts_y-c0(2));
                % get mag field 
                load physicalConstants-SI mu0
                Bth = mu0 * (obj.I/obj.N) ./ (2 * pi * r); % Field in theta direction
                Bth(r <= obj.d0/2) = 0; 
                Bx = Bth .* -sin(th); % Field in x and y direction 
                By = Bth .* cos(th); 
                % Find total field
                Bx_t = Bx_t + Bx;
                By_t = By_t + By; 
            end
            Bmag = sqrt(Bx_t.^2 + By_t.^2); % Total field magnitude
            if obj.isExploder
                % (3a) Get mag field around cathode for exploder array
                c0 = obj.pos0(end,:); 
                r = sqrt((pts_x-c0(1)).^2 + (pts_y-c0(2)).^2); 
                th = obj.getth(pts_x-c0(1),pts_y-c0(2));
                Bth = mu0 * (-obj.I) ./ (2 * pi * r); % Field in theta direction
                Bth(r <= obj.d0/2) = 0;
                Bx = Bth .* -sin(th); % Field in x and y direction
                By = Bth .* cos(th);
                % Find total field
                Bx_t = Bx_t + Bx;
                By_t = By_t + By;
                Bmag = sqrt(Bx_t.^2 + By_t.^2); % Total field magnitude
                Bmag(r <= obj.d0/2) = 0;
            end
            B = Bmag;
        end
        function  out = getth(~,x,y)
            % Returns azmuthal angle for a point(x,y)
            th = atan(y./x);
            th(y >= 0 & x >= 0) = th(y >= 0 & x >= 0);
            th(y >= 0 & x < 0) = 1 * pi + th(y >= 0 & x < 0);
            th(y < 0 & x < 0) = pi + th(y < 0 & x < 0);
            th(y < 0 & x >= 0) = 2 * pi + th(y < 0 & x >= 0);
            out = th;
        end
        function showBField(obj,M,L)
            % shows the total B field for the grid
            % Inputs:
            % M = number of grid points in one direction
            % L = length of gdif box [m]
            % (1a) Create box 
            x = linspace(-1,1,M) * L; 
            [xx,yy] = meshgrid(x,x); 
            [Bx,By,Bmag] = obj.getMagField(xx, yy);
            % (4a) Plot the field
            imagesc((Bmag)); colormap('hot'); set(gca,'YDir','normal');hold on;
            set(gca,'YTickLabel',round((get(gca,'YTick')*1/M*L-L/2)*1e3))
            set(gca,'XTickLabel',round((get(gca,'XTick')*1/M*L-L/2)*1e3))
            xlabel('$x, mm$','Interpreter','latex'); ylabel('$y, mm$','Interpreter','latex');
            cb = colorbar; ylabel(cb,'$|B|,T$','Interpreter','latex');
            
            % (4b) Show vectors
            smp = M/20;
            qX = 1:M; [qX,qY] = meshgrid(qX,qX); 
            quiver(qX(1:smp:end,1:smp:end), qY(1:smp:end,1:smp:end),...
                Bx(1:smp:end,1:smp:end) ,By(1:smp:end,1:smp:end),'color','b',...
                'Linewidth',2,...
                'DisplayName','B'); 

            % (4c) Show position of wires
            p = (obj.pos0 + (1 * L)/M) ./ (2*L) * M + M/2; % Scale
            d = obj.d0 ./ (2 * L) * M;
            viscircles(p,0.5*d * ones(size(obj.pos0,1),1),'Color','g'); 
            
            
            % (4d) Show force direction
            [fx,fy] = obj.getJxBDirection(Bx,By,Bmag);
            
            quiver(qX(1:smp:end,1:smp:end), qY(1:smp:end,1:smp:end),...
                fx(1:smp:end,1:smp:end) ,fy(1:smp:end,1:smp:end),'color','r',...
                'Linewidth',2,...
                'DisplayName','$\hat{j} \times {\bf B}$'); 
            
            formatPlots(); axis square
            legend('TextColor','w'); 
        end
        function [fx,fy] = getJxBDirection(obj,Bx,By,B)
            % returns the of the j x B force on a unit j 
            % By calculating the force acting on a unit current density 
            % in the z-direction due to the 
            bx = Bx ./ 1; bx = bx(:); % 
            by = By ./ 1; by = by(:); 
            b = [bx, by, zeros(size(bx))]; 
            j = [bx*0, by*0, ones(size(bx))];
            F = cross(j,b,2); % lorentz force
            fx = F(:,1); fy = F(:,2); 
            fx = reshape(fx,[size(B)]); fy = reshape(fy,[size(B)]);
        end
    end
        
end