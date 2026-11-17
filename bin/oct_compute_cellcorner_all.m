function [c1x,c1y,c2x,c2y,c3x,c3y,c4x,c4y,wdth,hgt] = oct_compute_cellcorner_all(x,y)

% -------------------------------------------------------------------------
%| oct_compute_cellcorner_all  (vectorised version of oct_compute_cellcorner)      |
%|                                                                        |
%| Returns, for ALL cells at once, the corner coordinates (c1..c4) and    |
%| the cell width/height that oct_compute_cellcorner(x,y,j,k,Nx,Ny) returns    |
%| one cell at a time. Every output is an (Ny x Nx) array. The arithmetic |
%| is the same (same operations in the same order) so results are         |
%| identical to the scalar routine.                                       |
%|                                                                        |
%| Requires Nx >= 3 and Ny >= 3 (callers fall back to the scalar routine  |
%| otherwise).                                                            |
% -------------------------------------------------------------------------

[Ny,Nx] = size(x);

xp = NaN(Ny+2,Nx+2);  yp = NaN(Ny+2,Nx+2);
xp(2:Ny+1,2:Nx+1) = x;  yp(2:Ny+1,2:Nx+1) = y;

% Corner estimates taken from the diagonal neighbours
[A1x,A1y] = nbcorner(-1,+1);   % neighbour (k-1,j+1)
[A2x,A2y] = nbcorner(+1,+1);   % neighbour (k+1,j+1)
[A3x,A3y] = nbcorner(+1,-1);   % neighbour (k+1,j-1)
[A4x,A4y] = nbcorner(-1,-1);   % neighbour (k-1,j-1)

c1x = zeros(Ny,Nx); c1y = c1x; c2x = c1x; c2y = c1x;
c3x = c1x; c3y = c1x; c4x = c1x; c4y = c1x;

r = 2:Ny-1;   c = 2:Nx-1;

% Interior cells
c1x(r,c) = A1x(r,c); c1y(r,c) = A1y(r,c);
c2x(r,c) = A2x(r,c); c2y(r,c) = A2y(r,c);
c3x(r,c) = A3x(r,c); c3y(r,c) = A3y(r,c);
c4x(r,c) = A4x(r,c); c4y(r,c) = A4y(r,c);

% Left edge (j == 1), k = 2..Ny-1
c1x(r,1) = A1x(r,1); c1y(r,1) = A1y(r,1);
c2x(r,1) = A2x(r,1); c2y(r,1) = A2y(r,1);
c3x(r,1) = 2*x(r,1) - c1x(r,1); c3y(r,1) = 2*y(r,1) - c1y(r,1);
c4x(r,1) = 2*x(r,1) - c2x(r,1); c4y(r,1) = 2*y(r,1) - c2y(r,1);

% Right edge (j == Nx), k = 2..Ny-1
c3x(r,Nx) = A3x(r,Nx); c3y(r,Nx) = A3y(r,Nx);
c4x(r,Nx) = A4x(r,Nx); c4y(r,Nx) = A4y(r,Nx);
c1x(r,Nx) = 2*x(r,Nx) - c3x(r,Nx); c1y(r,Nx) = 2*y(r,Nx) - c3y(r,Nx);
c2x(r,Nx) = 2*x(r,Nx) - c4x(r,Nx); c2y(r,Nx) = 2*y(r,Nx) - c4y(r,Nx);

% Bottom edge (k == 1), j = 2..Nx-1
c2x(1,c) = A2x(1,c); c2y(1,c) = A2y(1,c);
c3x(1,c) = A3x(1,c); c3y(1,c) = A3y(1,c);
c4x(1,c) = 2*x(1,c) - c2x(1,c); c4y(1,c) = 2*y(1,c) - c2y(1,c);
c1x(1,c) = 2*x(1,c) - c3x(1,c); c1y(1,c) = 2*y(1,c) - c3y(1,c);

% Top edge (k == Ny), j = 2..Nx-1
c4x(Ny,c) = A4x(Ny,c); c4y(Ny,c) = A4y(Ny,c);
c1x(Ny,c) = A1x(Ny,c); c1y(Ny,c) = A1y(Ny,c);
c2x(Ny,c) = 2*x(Ny,c) - c4x(Ny,c); c2y(Ny,c) = 2*y(Ny,c) - c4y(Ny,c);
c3x(Ny,c) = 2*x(Ny,c) - c1x(Ny,c); c3y(Ny,c) = 2*y(Ny,c) - c1y(Ny,c);

% Corner (k=1,j=1)
k = 1; j = 1; x0 = x(k,j); y0 = y(k,j);
c2x(k,j) = A2x(k,j); c2y(k,j) = A2y(k,j);
c4x(k,j) = 2*x0 - c2x(k,j); c4y(k,j) = 2*y0 - c2y(k,j);
c3x(k,j) = x0 - (c2y(k,j)-y0); c3y(k,j) = y0 + (c2x(k,j)-x0);
c1x(k,j) = 2*x0 - c3x(k,j); c1y(k,j) = 2*y0 - c3y(k,j);

% Corner (k=Ny,j=1)
k = Ny; j = 1; x0 = x(k,j); y0 = y(k,j);
c1x(k,j) = A1x(k,j); c1y(k,j) = A1y(k,j);
c3x(k,j) = 2*x0 - c1x(k,j); c3y(k,j) = 2*y0 - c1y(k,j);
c2x(k,j) = x0 - (y0-c1y(k,j)); c2y(k,j) = y0 + (x0-c1x(k,j));
c4x(k,j) = 2*x0 - c2x(k,j); c4y(k,j) = 2*y0 - c2y(k,j);

% Corner (k=1,j=Nx)
k = 1; j = Nx; x0 = x(k,j); y0 = y(k,j);
c3x(k,j) = A3x(k,j); c3y(k,j) = A3y(k,j);
c2x(k,j) = x0 - (c3y(k,j)-y0); c2y(k,j) = y0 + (c3x(k,j)-x0);
c1x(k,j) = 2*x0 - c3x(k,j); c1y(k,j) = 2*y0 - c3y(k,j);
c4x(k,j) = 2*x0 - c2x(k,j); c4y(k,j) = 2*y0 - c2y(k,j);

% Corner (k=Ny,j=Nx)
k = Ny; j = Nx; x0 = x(k,j); y0 = y(k,j);
c4x(k,j) = A4x(k,j); c4y(k,j) = A4y(k,j);
c3x(k,j) = x0 - (c4y(k,j)-y0); c3y(k,j) = y0 + (c4x(k,j)-x0);
c1x(k,j) = 2*x0 - c3x(k,j); c1y(k,j) = 2*y0 - c3y(k,j);
c2x(k,j) = 2*x0 - c4x(k,j); c2y(k,j) = 2*y0 - c4y(k,j);

wdth = sqrt((c1x-c4x).^2+(c1y-c4y).^2);
hgt  = sqrt((c2x-c1x).^2+(c2y-c1y).^2);

    function [Ax,Ay] = nbcorner(dk,dj)
        xt = xp(2+dk:Ny+1+dk, 2+dj:Nx+1+dj);
        yt = yp(2+dk:Ny+1+dk, 2+dj:Nx+1+dj);
        d = xt - x;
        m = abs(d) > 270;
        xt(m) = xt(m) - 360*sign(d(m));
        Ax = 0.5*(xt + x);
        Ay = 0.5*(yt + y);
    end

end
