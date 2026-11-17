function mask = oct_clean_mask(x,y,mask,bound_ingrid,lim,offset)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  29-Mar-2013 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| This function checks all the wet cells in a 2D mask array and          |  
%| determines if they lie outside the boundary polygons or not            |
%|                                                                        |
%| mask = oct_clean_mask(x,y,mask,bound_ingrid,lim,offset)                    |
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      |
%|  y            : A 2D array specifying the longitudes of each cell      |
%|  mask         : A 2D mask array of size (Nx,Ny) of initial mask values |
%|                 with wet cells having a flag of 1 and dry cells a flag |
%|                 value of 0. The initial mask can be all set to 1 or for|
%|		   better efficiency be based on the bathymetry data      |
%|  bound_ingrid : Data structure array of boundary polygons that lie     |
%|                 inside the grid                                        |
%|  lim          : Fraction value between 0 and 1 that define the cut-off |
%|                 limit for the proportion of a cell that has to lie     | 
%|                 inside boundary polygons for the cell to be marked dry |
%|  offset       : An additional width around the boundaries to check if  |
%|                 wet cells are being crossed by the boundary. Typically |
%|                 set it to largest grid resolution                      |
%|                                                                        |
%|  OUTPUT                                                                |
%|    mask       : New land/sea mask array that is generated after        |
%|                 checking with the boundary polygons                    |
%|                                                                        |
% -------------------------------------------------------------------------


%@@@ Determine array limits

N1 = length(bound_ingrid);
[Ny,Nx] = size(x);

%@@@ The fast algorithm below needs a rectilinear grid (all rows of x equal
%@@@ and all columns of y equal, increasing). Anything else (curvilinear
%@@@ grids, tiny grids) is handled by the original routine.

x1 = x(1,:);
y1 = y(:,1);

if (~(Nx >= 3 && Ny >= 3 && N1 > 0 && issorted(x1) && issorted(y1) && ...
        isequal(x,repmat(x1,Ny,1)) && isequal(y,repmat(y1,1,Nx))))
    mask = oct_clean_mask_orig(x,y,mask,bound_ingrid,lim,offset);
    return;
end;

%@@@ Corners of every cell (vectorised). Cells that are plain rectangles
%@@@ whose sampling points only depend on their column / row are processed
%@@@ together ("regular" cells). The cells on the grid edges (whose corners
%@@@ are extrapolated) and any non-rectangular cell are "special" and are
%@@@ processed one by one with the original algorithm.

[c1x,c1y,c2x,c2y,c3x,c3y,c4x,c4y] = oct_compute_cellcorner_all(x,y);

box_cell = (c1x==c2x) & (c3x==c4x) & (c1y==c4y) & (c2y==c3y) & ...
           (c1x>c4x) & (c2y>c1y);
xmin_c = min(min(c1x,c2x),min(c3x,c4x));
xmax_c = max(max(c1x,c2x),max(c3x,c4x));
ymin_c = min(min(c1y,c2y),min(c3y,c4y));
ymax_c = max(max(c1y,c2y),max(c3y,c4y));

regular = box_cell & ...
          (xmin_c == repmat(xmin_c(2,:),Ny,1)) & ...
          (xmax_c == repmat(xmax_c(2,:),Ny,1)) & ...
          (ymin_c == repmat(ymin_c(:,2),1,Nx)) & ...
          (ymax_c == repmat(ymax_c(:,2),1,Nx));
regular([1 Ny],:) = false;
regular(:,[1 Nx]) = false;

special = find(~regular & (mask == 1))';

%@@@ Sampling points of the regular cells: 6 x 6 points per cell,
%@@@ equally spaced between the cell limits (same call to linspace as the
%@@@ original code, one column / row at a time)

X6t = linspace(xmin_c(2,:)',xmax_c(2,:)',6)';     % 6 x Nx
Y6t = linspace(ymin_c(:,2),ymax_c(:,2),6)';       % 6 x Ny
jxv = ceil((1:36)'/6);                            % x index of point p
iyv = mod((0:35)',6)+1;                           % y index of point p
Npts = 36;

%@@@ Boundary limits (with the offset), as in the original code

bw = zeros(1,N1); be = bw; bs = bw; bn = bw;
tw = bw; te = bw; ts = bw; tn = bw;
for i = 1:N1
    tw(i) = bound_ingrid(i).west;   te(i) = bound_ingrid(i).east;
    ts(i) = bound_ingrid(i).south;  tn(i) = bound_ingrid(i).north;
end;
bw = tw-offset;  be = te+offset;  bs = ts-offset;  bn = tn+offset;

%@@@ status of each sampling point (inside any of the boundaries so far)

status = false(Npts,Ny*Nx);
tolbb = 1e-9;
itmp = 0;

%@@@ Loop through all the boundaries

for i = 1:N1

    %@@@ Cells (centres) within the boundary limits + offset. The grid is
    %@@@ rectilinear, so these are a block of rows and columns

    c_lo = find(x1 >= bw(i),1,'first');
    c_hi = find(x1 <= be(i),1,'last');
    r_lo = find(y1 >= bs(i),1,'first');
    r_hi = find(y1 <= bn(i),1,'last');

    if (~isempty(c_lo) && ~isempty(c_hi) && ~isempty(r_lo) && ~isempty(r_hi) ...
            && c_hi >= c_lo && r_hi >= r_lo)

        %@@@ Wet regular cells in that block

        blk = (mask(r_lo:r_hi,c_lo:c_hi) == 1) & regular(r_lo:r_hi,c_lo:c_hi);
        [rl,cl] = find(blk);

        if (~isempty(rl))

            rl = rl+r_lo-1;
            cl = cl+c_lo-1;
            idx = rl + (cl-1)*Ny;

            PX = X6t(jxv,cl);                  % 36 x nc
            PY = Y6t(iyv,rl);                  % 36 x nc
            St = status(:,idx);

            %@@@ Only points inside the (padded) bounding box of the boundary
            %@@@ and not yet flagged can be inside the boundary polygon

            test = (PX >= tw(i)-tolbb) & (PX <= te(i)+tolbb) & ...
                   (PY >= ts(i)-tolbb) & (PY <= tn(i)+tolbb) & ~St;

            if (any(test(:)))
                inout = oct_inpolygon_fast(PX(test),PY(test),bound_ingrid(i).x,...
                    bound_ingrid(i).y);
                St(test) = inout;
                status(:,idx) = St;
            end;

            %@@@ Proportion of each cell covered by boundaries; mark the
            %@@@ cell dry when it exceeds the user specified limit

            mask_obstr = round(sum(St,1)/Npts*10)/10;
            mask(idx(mask_obstr >= lim)) = 0;

        end;
    end;

    %@@@ Counter to update proportion of land sea mask clean up

    itmp_prev = itmp;
    itmp = floor(i/N1*100);
    if (mod(itmp,5) == 0 && itmp_prev ~= itmp)
        fprintf(1,'Completed %d per cent of land sea mask clean up\n',itmp);
    end;

end; %@@@ Corresponds to for loop of all the boundaries

%@@@ Special cells: original cell-by-cell algorithm, restricted to the
%@@@ boundaries whose (limits + offset) contain the cell centre

for lin = special

    k = mod(lin-1,Ny)+1;
    j = (lin-k)/Ny+1;

    if (mask(k,j) ~= 1)
        continue;
    end;

    cand = find(x(k,j) >= bw & x(k,j) <= be & y(k,j) >= bs & y(k,j) <= bn);
    if (isempty(cand))
        continue;
    end;

    c1 = [c1x(k,j) c1y(k,j)];  c2 = [c2x(k,j) c2y(k,j)];
    c3 = [c3x(k,j) c3y(k,j)];  c4 = [c4x(k,j) c4y(k,j)];
    xmin = min([c1(1) c2(1) c3(1) c4(1)]);
    xmax = max([c1(1) c2(1) c3(1) c4(1)]);
    ymin = min([c1(2) c2(2) c3(2) c4(2)]);
    ymax = max([c1(2) c2(2) c3(2) c4(2)]);

    xtt = linspace(xmin,xmax,6);
    ytt = linspace(ymin,ymax,6);
    [xtt2,ytt2] = meshgrid(xtt,ytt);
    px1 = [c4(1) c1(1) c2(1) c3(1) c4(1)]';
    py1 = [c4(2) c1(2) c2(2) c3(2) c4(2)]';
    in_cell = inpolygon(xtt2,ytt2,px1,py1);
    loc = find(in_cell > 0);
    xt = xtt2(loc);
    yt = ytt2(loc);
    st = zeros(size(xt));
    Na = numel(xt);

    for i = cand
        loc = find(st == 0);
        if (~isempty(loc))
            inout = oct_inpolygon_fast(xt(loc),yt(loc),bound_ingrid(i).x,...
                bound_ingrid(i).y);
            st(loc) = inout;
        end;
        mo = round(length(find(st > 0))/Na*10)/10;
        if (mo >= lim)
            mask(k,j) = 0;
            break;
        end;
    end;
end;

return;
