function [sx,sy] = oct_create_obstr(x,y,bound,mask,offset_left,offset_right)

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
%| This routine generates the 2D obstruction grid in x and y given a 2D   |
%| mask and set of boundary polygons. Obstructions are only generated for |
%| wet cells and obstructions for cells on either side of a dry cell are  |
%| also set to 0 (to prevent spurious suppression of swell near the coast)|  
%| The routine allows for the possibility of curvilinear coordinates and  |
%| locally rotates the grid to align the coordinates in lat/lon space with|
%| local p/q space which is determined by the 2D x and y matrices         |
%|                                                                        |
%| [sx,sy] = create_obstr_curv(x,y,bound,mask,coords,offset_left,...      |
%|                                                   offset_right)        |
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      | 
%|  y            : A 2D array specifying the lattitudes of each cell      | 
%|  bound        : Data structure array of boundary polygons              |
%|  mask         : 2D array of size (Nx,Ny) that determines land/sea mask |
%|  offset_left  : Flag to determine if neighbor to the left/down in x/y  |
%|                 should be considered. (0/1 = no/yes)                   |
%|  offset_right : Similar for neighbor to the right/up in x/y            |
%|                                                                        |
%| OUTPUT                                                                 |
%|  sx,sy        : 2D obstruction grids of size (Nx,Ny) for obstructions  |
%|                 in x and y. Values range from 0 for no obstruction to 1| 
%|                 for full obstruction                                   | 
%|                                                                        |
% -------------------------------------------------------------------------

%@@@ Initialize variables

[Ny,Nx] = size(x);

N = length(bound);

%@@@ The fast algorithm needs a rectilinear grid (all rows of x equal and all
%@@@ columns of y equal, increasing), more than 2 cells in each direction and
%@@@ boundaries whose stored length agrees with their coordinates. Anything
%@@@ else is handled by the original routine.

x1 = x(1,:);
y1 = y(:,1);
use_fast = (Nx >= 3 && Ny >= 3 && N > 0 && issorted(x1) && issorted(y1) && ...
            isequal(x,repmat(x1,Ny,1)) && isequal(y,repmat(y1,1,Nx)));
if (use_fast)
    use_fast = isequal(arrayfun(@(b) numel(b.x),bound(:)'),[bound.n]);
end;
if (~use_fast)
    [sx,sy] = oct_create_obstr_orig(x,y,bound,mask,offset_left,offset_right);
    return;
end;

sx = zeros(Ny,Nx);
sy = sx;

cell_bnd = mask;

loc_wet = find(mask ~= 0);
N_wet = length(loc_wet);
Nb = Nx*Ny;

fprintf(1,' Total Number of cells = %d\n',Nb);
fprintf(1,'   Number of wet cells = %d\n',N_wet);

%@@@ Cell corners, widths and heights (vectorised)

[c1x,c1y,c2x,c2y,c3x,c3y,c4x,c4y,Wc,Hc] = oct_compute_cellcorner_all(x,y);

%@@@ Preparing the boundaries: all boundary points in two long vectors

fprintf(1,'Preparing the boundaries \n');

%@@@ forced to columns: some boundaries can be stored as rows
bnd_x = cell2mat(cellfun(@(v) v(:),{bound.x},'UniformOutput',false)');
bnd_y = cell2mat(cellfun(@(v) v(:),{bound.y},'UniformOutput',false)');
bnd_indx = repelem((1:N)',[bound.n]');
Mpts = numel(bnd_x);

%@@@ Which cells can contain each boundary point? For a rectilinear grid
%@@@ the cell of a point follows from its column and row. A margin is kept
%@@@ so that points on (or numerically at) a cell edge are candidates for
%@@@ both neighbouring cells; the exact test is done with inpolygon below.

mrg = 1e-7;

xlo = c4x(2,:)';  xhi = c1x(2,:)';
ylo = c4y(:,2);   yhi = c2y(:,2);

[~,jlo] = histc(bnd_x,[-Inf; xhi+mrg; Inf]);       % first column with xhi+mrg > x
[~,jt]  = histc(bnd_x,[-Inf; xlo-mrg; Inf]);       % columns with xlo-mrg <= x : 1..jt-1
jhi = jt-1;
[~,klo] = histc(bnd_y,[-Inf; yhi+mrg; Inf]);
[~,kt]  = histc(bnd_y,[-Inf; ylo-mrg; Inf]);
khi = kt-1;

jlo = min(jlo,Nx+1);   klo = min(klo,Ny+1);
jhi = min(jhi,Nx);     khi = min(khi,Ny);

ok = (jlo <= jhi) & (klo <= khi);
nj = jhi-jlo+1;  nk = khi-klo+1;

pair_cell = [];
pair_pt = [];
for dj = 0:max([0; nj(ok)])-1
    for dk = 0:max([0; nk(ok)])-1
        sel = find(ok & nj > dj & nk > dk);
        pair_cell = [pair_cell; (klo(sel)+dk) + (jlo(sel)+dj-1)*Ny];
        pair_pt = [pair_pt; sel];
    end;
end;
clear jlo jhi klo khi jt kt ok nj nk;

%@@@ Keep wet cells only and group the candidate points by cell

keep = (mask(pair_cell) ~= 0);
pair_cell = pair_cell(keep);
pair_pt = pair_pt(keep);
[pair_cell,ord] = sort(pair_cell);
pair_pt = pair_pt(ord);
clear keep ord;

%@@@ Cells whose polygon is not inside the representative box (the cells
%@@@ in the grid corners when dx ~= dy) can contain points that the binning
%@@@ above misses: they get every boundary point as candidate (as the
%@@@ original code does).

xmin_c = min(min(c1x,c2x),min(c3x,c4x));  xmax_c = max(max(c1x,c2x),max(c3x,c4x));
ymin_c = min(min(c1y,c2y),min(c3y,c4y));  ymax_c = max(max(c1y,c2y),max(c3y,c4y));
outside_box = (xmin_c < repmat(xlo',Ny,1)-0.5*mrg) | (xmax_c > repmat(xhi',Ny,1)+0.5*mrg) | ...
              (ymin_c < repmat(ylo,1,Nx)-0.5*mrg)  | (ymax_c > repmat(yhi,1,Nx)+0.5*mrg);
special_cells = find(outside_box & (mask ~= 0))';
clear xmin_c xmax_c ymin_c ymax_c outside_box;

[ucell,ustart] = unique(pair_cell,'first');
uend = [ustart(2:end)-1; numel(pair_cell)];

special_flag = false(Ny*Nx,1);
special_flag(special_cells) = true;

%@@@ Cells to process: those with candidate points, plus the special cells

extra_special = setdiff(special_cells,ucell(:)');
work_cells = [ucell(:)' extra_special];
n_ucell = numel(ucell);

%@@@ Slot map: cells with boundaries get a slot in a compact struct array.
%@@@ Slot 1 is an empty dummy cell used for all the cells without boundaries.

empty_cell = struct('nx',0,'ny',0,'south_lim',[],'north_lim',[],...
    'east_lim',[],'west_lim',[],'bndx',[],'bndy',[]);
cell = repmat(empty_cell,1,1+numel(work_cells));
sl = ones(Ny,Nx,'int32');

fprintf(1,'Loop through the wet cells to identify boundaries \n');

nslot = 1;
cell_bnd(loc_wet) = 0;

for iw = 1:numel(work_cells)

    lin = work_cells(iw);
    k = mod(lin-1,Ny)+1;
    j = (lin-k)/Ny+1;

    if (special_flag(lin) || iw > n_ucell)
        cand = (1:Mpts)';
    else
        cand = pair_pt(ustart(iw):uend(iw));
    end;

    %@@@ The cell polygon, orientation and size

    px = [c4x(k,j) c1x(k,j) c2x(k,j) c3x(k,j) c4x(k,j)]';
    py = [c4y(k,j) c1y(k,j) c2y(k,j) c3y(k,j) c4y(k,j)]';
    angle = atan2(c1y(k,j)-c4y(k,j),c1x(k,j)-c4x(k,j));
    cwidth = Wc(k,j);
    cheight = Hc(k,j);

    in_box = inpolygon(bnd_x(cand),bnd_y(cand),px,py);
    pin = cand(in_box > 0);
    if (isempty(pin))
        continue;
    end;

    [bnds,~,grp] = unique(bnd_indx(pin));
    Nbnds = length(bnds);
    cell_bnd(lin) = Nbnds;

    %@@@ Rotation Matrix and origin (for each cell the first corner is taken as
    %@@@ the origin about which the rotation takes place)

    RM = [cos(angle) -sin(angle) ; sin(angle) cos(angle)];
    x0 = px(1);
    y0 = py(1);

    nslot = nslot+1;
    c = empty_cell;

    for i = 1:Nbnds

        indx_bnd = bnds(i);
        xt = bnd_x(pin(grp == i));
        yt = bnd_y(pin(grp == i));

        % rotate points within the cell (including the corners)

        tmp = [(xt-x0) (yt-y0)]*RM;
        xt = tmp(:,1);
        yt = tmp(:,2);

        % Boundary segments

        south_limit = min(yt)/cheight;
        north_limit = max(yt)/cheight;
        west_limit = min(xt)/cwidth;
        east_limit = max(xt)/cwidth;

        south_limit = max(0.0,south_limit);
        south_limit = min(1.0,south_limit);
        north_limit = max(0.0,north_limit);
        north_limit = min(1.0,north_limit);
        west_limit = max(0.0,west_limit);
        west_limit = min(1.0,west_limit);
        east_limit = max(0.0,east_limit);
        east_limit = min(1.0,east_limit);

        c.nx = c.nx+1;
        c.south_lim(c.nx) = south_limit;
        c.north_lim(c.nx) = north_limit;
        c.bndx(c.nx) = indx_bnd;
        c.ny = c.ny+1;
        c.east_lim(c.ny) = east_limit;
        c.west_lim(c.ny) = west_limit;
        c.bndy(c.ny) = indx_bnd;
    end;

    cell(nslot) = c;
    sl(k,j) = nslot;
end;
cell = cell(1:nslot);
clear bnd_x bnd_y bnd_indx pair_cell pair_pt;

%@@@ Loop through all the wet cells with boundaries and move boundary segments 
%@@@ that are part of the same boundary and cross neighboring cells. This is 
%@@@ done to prevent double counting in building the obstruction grids. Only the
%@@@  immediate neighbors are checked when moving the segments

loc_bnd = find(cell_bnd ~= 0);
N_bnd = length(loc_bnd);
[row_bnd,column_bnd] = ind2sub(size(cell_bnd),loc_bnd);
fprintf(1,'Number of wet cells enclosing boundaries = %d\n',N_bnd);

fprintf(1,'Computing subgrid obstruction masks for these cells \n');

for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);
    
    %@@@ First check the neighbors in x direction
            
    if (j < Nx)
        jj = j+1;

        %@@@ Check to see that boundary segments are present in both cells

        if (cell(sl(k,j)).nx~=0 && cell(sl(k,jj)).nx~=0)

            %@@@ Save information to temporary variables

            set1 = cell(sl(k,j));
            set2 = cell(sl(k,jj));
            found_common = 0;

            %@@@ Loop through boundary segments and move segments  
            %@@@ of common boundaries to the cell with the larger 
            %@@@ segment

            for l = 1:set1.nx
                for m = 1:set2.nx
                    if (set1.bndx(l) == set2.bndx(m))
                        if ((set1.north_lim(l)-set1.south_lim(l)) >= ...
                                (set2.north_lim(m)-set2.south_lim(m)))
                            set1.north_lim(l) = max([set1.north_lim(l) ...
                                set2.north_lim(m)]);
                            set1.south_lim(l) = min([set1.south_lim(l) ...
                                set2.south_lim(m)]);
                            for n = (m+1):set2.nx
                                set2.bndx(n-1) = set2.bndx(n);
                                set2.north_lim(n-1) = set2.north_lim(n);
                                set2.south_lim(n-1) = set2.south_lim(n);
                            end;
                            set2.nx = set2.nx-1;
                        else
                            set2.north_lim(m) = max([set1.north_lim(l) ...
                                set2.north_lim(m)]);
                            set2.south_lim(m) = min([set1.south_lim(l) ...
                                set2.south_lim(m)]);
                            for n = (l+1):set1.nx
                                set1.bndx(n-1) = set1.bndx(n);
                                set1.north_lim(n-1) = set1.north_lim(n);
                                set1.south_lim(n-1) = set1.south_lim(n);
                            end;
                            set1.nx = set1.nx-1;
                        end;
                        found_common = 1;
                        break;
                    end;
                end;
            end;

            %@@@ Write cell information back from temporary variables 
            %@@@ if common boundaries were found

            if (found_common == 1)
                cell(sl(k,j)).bndx = [];
                cell(sl(k,j)).north_lim = [];
                cell(sl(k,j)).south_lim = [];
                cell(sl(k,jj)).bndx = [];
                cell(sl(k,jj)).north_lim = [];
                cell(sl(k,jj)).south_lim = [];
                        
                cell(sl(k,j)).nx = set1.nx;
                cell(sl(k,j)).bndx = set1.bndx(1:set1.nx);
                cell(sl(k,j)).north_lim = set1.north_lim(1:set1.nx);
                cell(sl(k,j)).south_lim = set1.south_lim(1:set1.nx);
                        
                cell(sl(k,jj)).nx = set2.nx;
                cell(sl(k,jj)).bndx = set2.bndx(1:set2.nx);
                cell(sl(k,jj)).north_lim = set2.north_lim(1:set2.nx);
                cell(sl(k,jj)).south_lim = set2.south_lim(1:set2.nx);
                
            end;

        end;

    end;
            
    %@@@ Now check the neighbors in y direction (comments are omitted
    %@@@ because the operations are same as for x)
            
    if (k < Ny)
        kk = k+1;

        if (cell(sl(k,j)).ny~=0 && cell(sl(kk,j)).ny~=0)
            set1 = cell(sl(k,j));
            set2 = cell(sl(kk,j));
            found_common = 0;

            for l = 1:set1.ny
                for m = 1:set2.ny

                    if (set1.bndy(l) == set2.bndy(m))
                        if ((set1.east_lim(l)-set1.west_lim(l)) >= ...
                                (set2.east_lim(m)-set2.west_lim(m)))
                            set1.east_lim(l) = max([set1.east_lim(l) ...
                                set2.east_lim(m)]);
                            set1.west_lim(l) = min([set1.west_lim(l) ...
                                set2.west_lim(m)]);
                            for n = (m+1):set2.ny
                                set2.bndy(n-1) = set2.bndy(n);
                                set2.east_lim(n-1) = set2.east_lim(n);
                                set2.west_lim(n-1) = set2.west_lim(n);
                            end;
                            set2.ny = set2.ny-1;
                        else
                            set2.east_lim(m) = max([set1.east_lim(l) ...
                                set2.east_lim(m)]);
                            set2.west_lim(m) = min([set1.west_lim(l) ...
                                set2.west_lim(m)]);
                            for n = (l+1):set1.ny
                                set1.bndy(n-1) = set1.bndy(n);
                                set1.east_lim(n-1) = set1.east_lim(n);
                                set1.west_lim(n-1) = set1.west_lim(n);
                            end;
                            set1.ny = set1.ny-1;
                        end;
                        found_common = 1;
                        break;
                    end;
                end;
            end;

            if (found_common == 1)
                cell(sl(k,j)).bndy = [];
                cell(sl(k,j)).east_lim = [];
                cell(sl(k,j)).west_lim = [];
                cell(sl(kk,j)).bndy = [];
                cell(sl(kk,j)).east_lim = [];
                cell(sl(kk,j)).west_lim = [];
            
                cell(sl(k,j)).ny = set1.ny;
                cell(sl(k,j)).bndy = set1.bndy(1:set1.ny);
                cell(sl(k,j)).east_lim = set1.east_lim(1:set1.ny);
                cell(sl(k,j)).west_lim = set1.west_lim(1:set1.ny);
                        
                cell(sl(kk,j)).ny = set2.ny;
                cell(sl(kk,j)).bndy = set2.bndy(1:set2.ny);
                cell(sl(kk,j)).east_lim = set2.east_lim(1:set2.ny);
                cell(sl(kk,j)).west_lim = set2.west_lim(1:set2.ny);                       
            end;

        end;

    end; %@@@ End of if condition to check neighbor in y direction        
    
end;  %@@@ For loop along all the wet cells with boundaries


%@@@ Loop through the cells a second time to reduce overlapping segments 
%@@@ into individual non-overlapping segments

for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);

    %@@@ First for x obstruction

    if (cell(sl(k,j)).nx ~= 0)

    %@@@ Reduction of overlapping segments only done if there are more than
    %@@@ 1 segment

        if (cell(sl(k,j)).nx > 1)
            n_segs = cell(sl(k,j)).nx;
            baseseg_n = cell(sl(k,j)).north_lim;
            baseseg_s = cell(sl(k,j)).south_lim;
            cell(sl(k,j)).north_lim = [];
            cell(sl(k,j)).south_lim = [];
            ind_segs=0;
            indseg_n = [];
            indseg_s = [];

            %@@@ Loop till all the segments in the cell done

            while (n_segs > 0)
                overlap_found = 0;                       
                if (n_segs > 1)
                    for l = 2:n_segs
                        if (baseseg_n(1) >= baseseg_s(l) && ...
                                baseseg_s(1) <= baseseg_n(l))
                            baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                            baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                            overlap_found = 1;
                            if (l == n_segs)
                                n_segs = n_segs-1;
                            else
                                for m = l+1:n_segs
                                    baseseg_n(m-1) = baseseg_n(m);
                                    baseseg_s(m-1) = baseseg_s(m);
                                end;
                                n_segs = n_segs-1;
                            end;
                            break;
                        end;
                    end;
                end;
                        
                if (n_segs == 1)
                    ind_segs = ind_segs+1;
                    indseg_n(ind_segs) = baseseg_n(1);
                    indseg_s(ind_segs) = baseseg_s(1);
                    n_segs = n_segs-1;
                else
                    if (overlap_found == 0)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        for l = 2:n_segs                                    
                            baseseg_n(l-1) = baseseg_n(l);
                            baseseg_s(l-1) = baseseg_s(l);
                        end;
                        n_segs = n_segs-1;
                    end;
                end;
            end;

            %@@@ Store the unique segments back in the cell variables

            cell(sl(k,j)).nx = ind_segs;
            cell(sl(k,j)).north_lim = indseg_n;
            cell(sl(k,j)).south_lim = indseg_s;

        end; %@@@ Corresponds to check for multiple segments in x

   end; %@@@ Corresponds to check for obstruction in x

   %@@@ Now check for overlapping segments in y 

   if (cell(sl(k,j)).ny ~= 0)

        if (cell(sl(k,j)).ny > 1)

            n_segs = cell(sl(k,j)).ny;
            baseseg_n = cell(sl(k,j)).east_lim;
            baseseg_s = cell(sl(k,j)).west_lim;
            cell(sl(k,j)).east_lim = [];
            cell(sl(k,j)).west_lim = [];
            ind_segs=0;
            indseg_n = [];
            indseg_s = [];

            while (n_segs > 0)
                overlap_found = 0;
                if (n_segs > 1)
                    for l = 2:n_segs
                        if (baseseg_n(1) >= baseseg_s(l) && ... 
                                baseseg_s(1) <= baseseg_n(l))
                            baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                            baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                            overlap_found = 1;
                            if (l == n_segs)
                                n_segs = n_segs-1;
                            else
                                for m = l+1:n_segs
                                    baseseg_n(m-1) = baseseg_n(m);
                                    baseseg_s(m-1) = baseseg_s(m);
                                end;
                                n_segs = n_segs-1;
                            end;
                            break;
                        end;
                    end;
                end;
                        
                if (n_segs == 1)
                    ind_segs = ind_segs+1;
                    indseg_n(ind_segs) = baseseg_n(1);
                    indseg_s(ind_segs) = baseseg_s(1);
                    n_segs = n_segs-1;
                else
                    if (overlap_found == 0)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        for l = 2:n_segs                                    
                            baseseg_n(l-1) = baseseg_n(l);
                            baseseg_s(l-1) = baseseg_s(l);
                        end;
                        n_segs = n_segs-1;
                    end;
                end;
            end;

            cell(sl(k,j)).ny = ind_segs;
            cell(sl(k,j)).east_lim = indseg_n;
            cell(sl(k,j)).west_lim = indseg_s;

        end;

   end; %@@@ End of overlapping segment computation in y        
    
end; %@@@ End of for loop along wet cells with boundaries


%@@@ Final loop through the cell rows and columns to construct the 
%@@@ obstruction grid accounting for neighboring cell information


for indx_bnd = 1:N_bnd
    
    j = column_bnd(indx_bnd);
    k = row_bnd(indx_bnd);

    %@@@ Computing x obstruction

    if (cell(sl(k,j)).nx ~=0)

        n_segs = cell(sl(k,j)).nx;
        baseseg_n = cell(sl(k,j)).north_lim;
        baseseg_s = cell(sl(k,j)).south_lim;
    
        no_boundary = 0;    

        %@@@ Compare boundary segments in cell with that of cell 
        %@@@ to left (if applicable)

        for off = 1:offset_left

            jj = j-off;
            if (jj >= 1)
                if (cell(sl(k,jj)).nx ~= 0)
                    set1 = cell(sl(k,jj)); 

                    %@@@ remove segments in shadow of previous cell  
                            
                    shadow_flags = zeros(size(baseseg_n));

                    for m = 1:n_segs
                        for l = 1:set1.nx
                            if (set1.north_lim(l) >= baseseg_n(m) ...
                                    && set1.south_lim(l) <= baseseg_s(m))
                                shadow_flags(m) = 1;
                                break;
                            end;
                        end;                                
                    end;

                    loc = find(shadow_flags == 0);
                    if (isempty(loc))
                        no_boundary = 1;
                        n_segs = 0;
                        baseseg_n = [];
                        baseseg_s = [];
                    elseif (length(loc) < n_segs)
                        tmp_n = baseseg_n;
                        tmp_s = baseseg_s;
                        baseseg_n = [];
                        baseseg_s = [];
                        baseseg_n = tmp_n(loc);
                        baseseg_s = tmp_s(loc);
                        n_segs = length(baseseg_n);
                        clear tmp_n tmp_s;
                    end;
                    clear loc shadow_flags;
                            
                    %@@@ Remove segments in previous cell that are
                    %@@@ shadows of present cell, 

                    if (~no_boundary)

                        shadow_flags = zeros(size(set1.north_lim));

                        for m = 1:set1.nx
                            for l = 1:n_segs
                                if (set1.north_lim(m) <= baseseg_n(l) ...
                                        && set1.south_lim(m) >= baseseg_s(l))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        loc = find(shadow_flags == 0);
                        if (~isempty(loc))
                            if (length(loc) < set1.nx)
                                tmp_n = set1.north_lim;
                                tmp_s = set1.south_lim;
                                set1.north_lim = [];
                                set1.south_lim = [];
                                set1.north_lim = tmp_n(loc);
                                set1.south_lim = tmp_s(loc);
                                set1.nx = length(set1.north_lim);
                                clear tmp_n tmp_s;
                            end;
                            n_segs = n_segs+set1.nx;
                            baseseg_n = [baseseg_n set1.north_lim];
                            baseseg_s = [baseseg_s set1.south_lim];
                        end;
                        clear loc shadow_flags;
                    end;                           
                end;
            end;
        end; %@@@ End of search of segments in the left cell

        %@@@ Repeat for cell to the right (if applicable)

        if (~no_boundary)

            for off = 1:offset_right
                jj = j+off;
                if (jj <= Nx)
                    if (cell(sl(k,jj)).nx ~= 0)
                        set1 = cell(sl(k,jj)); 

                        %@@@ See if segments lie in shadow zone to 
                        %@@@ segments in right cell

                        shadow_flags = zeros(size(baseseg_n));
                        for m = 1:n_segs
                            for l = 1:set1.nx
                                if (set1.north_lim(l) >= baseseg_n(m) ...
                                        && set1.south_lim(l) <= baseseg_s(m))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        %@@@ Set no_boundary flag if all segments removed

                        loc = find(shadow_flags == 0);
                        if (isempty(loc))
                            no_boundary = 1;
                            n_segs = 0;
                            baseseg_n = [];
                            baseseg_s = [];
                        elseif (length(loc) < n_segs)
                            tmp_n = baseseg_n;
                            tmp_s = baseseg_s;
                            baseseg_n = [];
                            baseseg_s = [];
                            baseseg_n = tmp_n(loc);
                            baseseg_s = tmp_s(loc);
                            n_segs = length(baseseg_n);
                            clear tmp_n tmp_s;
                        end;

                        clear loc shadow_flags;

                        %@@@ Remove segments from right cell that 
                        %@@@ are in the shadow zone 
                            
                        if (~no_boundary)
                            shadow_flags = zeros(size(set1.north_lim));
                            for m = 1:set1.nx
                                for l = 1:n_segs
                                    if (set1.north_lim(m) <= baseseg_n(l) ...
                                            && set1.south_lim(m) >= baseseg_s(l))
                                        shadow_flags(m) = 1;
                                        break;
                                    end;
                                end;                                
                            end;

                            loc = find(shadow_flags == 0);
                            if (~isempty(loc))
                                if (length(loc) < set1.nx)
                                    tmp_n = set1.north_lim;
                                    tmp_s = set1.south_lim;
                                    set1.north_lim = [];
                                    set1.south_lim = [];
                                    set1.north_lim = tmp_n(loc);
                                    set1.south_lim = tmp_s(loc);
                                    set1.nx = length(set1.north_lim);
                                    clear tmp_n;
                                    clear tmp_s;
                                end;
                                n_segs = n_segs+set1.nx;
                                baseseg_n = [baseseg_n set1.north_lim];
                                baseseg_s = [baseseg_s set1.south_lim];
                            end;
                            clear loc;
                            clear shadow_flags;
                        end;                           
                    end;
                end;
            end;
        end; %@@@ End of check of cells to the right
                
        %@@@ Now ready to build obstruction grid from the total 
        %@@@ set of segments

        if (~no_boundary)

            %@@@ Obstruction grid straightforward if n_segs == 1

            if (n_segs==1)
                
                sx(k,j) = baseseg_n(1)-baseseg_s(1);

            else    

                %@@@ Remove overlapping segments

                ind_segs=0;
                indseg_n = [];
                indseg_s = [];
                
                while (n_segs > 0)
                    overlap_found = 0;
                    if (n_segs > 1)
                        for l = 2:n_segs
                            if (baseseg_n(1) >= baseseg_s(l) && ...
                                    baseseg_s(1) <= baseseg_n(l))
                                baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                                baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                                overlap_found = 1;
                                if (l == n_segs)
                                    n_segs = n_segs-1;
                                else
                                    for m = l+1:n_segs
                                        baseseg_n(m-1) = baseseg_n(m);
                                        baseseg_s(m-1) = baseseg_s(m);
                                    end;
                                    n_segs = n_segs-1;
                                end;
                                break;
                            end;
                        end;
                    end;
                        
                    if (n_segs == 1)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        n_segs = n_segs-1;
                    else
                        if (overlap_found == 0)
                            ind_segs = ind_segs+1;
                            indseg_n(ind_segs) = baseseg_n(1);
                            indseg_s(ind_segs) = baseseg_s(1);
                            for l = 2:n_segs                                    
                                baseseg_n(l-1) = baseseg_n(l);
                                baseseg_s(l-1) = baseseg_s(l);
                            end;
                            n_segs = n_segs-1;
                        end;
                    end;
                end;
      
                %@@@ Compute the obstruction values from the 
                %@@@ independant segments

                for l = 1:ind_segs
                    sx(k,j) = sx(k,j) + (indseg_n(l)-indseg_s(l));
                end;

                clear baseseg_n baseseg_s baseseg_bound;

            end; %@@@ End of if statement to check number of segments

        end; %@@@ End of if statement to check boundaries in cell
 
    end; %@@@ End of if statement to check obstructions along x
            
    %@@@ Now computing obstruction in y in a similar fashion to x

    if (cell(sl(k,j)).ny~=0)

        n_segs = cell(sl(k,j)).ny;
        baseseg_n = cell(sl(k,j)).east_lim;
        baseseg_s = cell(sl(k,j)).west_lim;
                
        no_boundary = 0;    

        %@@@ First check with boundaries in cell below (if applicable)

        for off = 1:offset_left

            kk = k-off;

            if (kk >= 1)
                if (cell(sl(kk,j)).ny ~= 0)
                    set1 = cell(sl(kk,j)); 
                            
                    %@@@ Like in x obstruction, remove boundaries 
                    %@@@ in shadow zone and combine segments
                           
                    shadow_flags = zeros(size(baseseg_n));
                    for m = 1:n_segs
                        for l = 1:set1.ny
                            if (set1.east_lim(l) >= baseseg_n(m) && ...
                                    set1.west_lim(l) <= baseseg_s(m))
                                shadow_flags(m) = 1;
                                break;
                            end;
                        end;                                
                    end;

                    loc = find(shadow_flags == 0);

                    if (isempty(loc))
                        no_boundary = 1;
                        n_segs = 0;
                        baseseg_n = [];
                        baseseg_s = [];
                    elseif (length(loc) < n_segs)
                        tmp_n = baseseg_n;
                        tmp_s = baseseg_s;
                        baseseg_n = [];
                        baseseg_s = [];
                        baseseg_n = tmp_n(loc);
                        baseseg_s = tmp_s(loc);
                        n_segs = length(baseseg_n);
                        clear tmp_n tmp_s;
                    end;

                    clear loc;
                    clear shadow_flags
                            
                    if (~no_boundary)
                        shadow_flags = zeros(size(set1.east_lim));
                        for m = 1:set1.ny
                            for l = 1:n_segs
                                if (set1.east_lim(m) <= baseseg_n(l) && ...
                                        set1.west_lim(m) >= baseseg_s(l))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;
                        loc = find(shadow_flags == 0);
                        if (~isempty(loc))
                            if (length(loc) < set1.ny)
                                tmp_n = set1.east_lim;
                                tmp_s = set1.west_lim;
                                set1.east_lim = [];
                                set1.west_lim = [];
                                set1.east_lim = tmp_n(loc);
                                set1.west_lim = tmp_s(loc);
                                set1.ny = length(set1.east_lim);
                                clear tmp_n;
                                clear tmp_s;
                            end;
                            n_segs = n_segs+set1.ny;
                            baseseg_n = [baseseg_n set1.east_lim];
                            baseseg_s = [baseseg_s set1.west_lim];
                        end;
                        clear loc;
                        clear shadow_flags;
                    end;                          
                end;
            end;
        end;  %@@@ End of accounting for boundaries in cell below

        %@@@ Now moving to boundaries in the cell above (if applicable)
                
        if (~no_boundary)
            
            for off = 1:offset_right
                kk = k+off;
                if (kk <= Ny)
                    if (cell(sl(kk,j)).ny ~= 0)
                        set1 = cell(sl(kk,j)); 

                        %@@@ Again follow the same procedure to merge boundaries 
                                
                        shadow_flags = zeros(size(baseseg_n));
                        for m = 1:n_segs
                            for l = 1:set1.ny
                                if (set1.east_lim(l) >= baseseg_n(m) && ...
                                        set1.west_lim(l) <= baseseg_s(m))
                                    shadow_flags(m) = 1;
                                    break;
                                end;
                            end;                                
                        end;

                        loc = find(shadow_flags == 0);
                        if (isempty(loc))
                            no_boundary = 1;
                            n_segs = 0;
                            baseseg_n = [];
                            baseseg_s = [];
                        elseif (length(loc) < n_segs)
                            tmp_n = baseseg_n;
                            tmp_s = baseseg_s;
                            baseseg_n = [];
                            baseseg_s = [];
                            baseseg_n = tmp_n(loc);
                            baseseg_s = tmp_s(loc);
                            n_segs = length(baseseg_n);
                            clear tmp_n tmp_s;
                        end;
                        clear loc shadow_flags;
                                
                        if (~no_boundary)
                            shadow_flags = zeros(size(set1.east_lim));
                            for m = 1:set1.ny
                                for l = 1:n_segs
                                    if (set1.east_lim(m) <= baseseg_n(l) && ...
                                            set1.west_lim(m) >= baseseg_s(l))
                                        shadow_flags(m) = 1;
                                        break;
                                    end;
                                end;                                
                            end;

                            loc = find(shadow_flags == 0);
                            if (~isempty(loc))
                                if (length(loc) < set1.ny)
                                    tmp_n = set1.east_lim;
                                    tmp_s = set1.west_lim;
                                    set1.east_lim = [];
                                    set1.west_lim = [];
                                    set1.east_lim = tmp_n(loc);
                                    set1.west_lim = tmp_s(loc);
                                    set1.ny = length(set1.east_lim);
                                    clear tmp_n tmp_s;
                                end;
                                n_segs = n_segs+set1.ny;
                                baseseg_n = [baseseg_n set1.east_lim];
                                baseseg_s = [baseseg_s set1.west_lim];
                            end;
                            clear loc shadow_flags;

                        end;                           
                    end;
                end;
            end;
        end;  %@@@ End of accounting for cell above
                
        %@@@ Building obstruction grid in y  
                
        if (~no_boundary)

            if (n_segs==1)
            
                sy(k,j) = baseseg_n(1)-baseseg_s(1);
            
            else    

                ind_segs=0;
                indseg_n = [];
                indseg_s = [];

                while (n_segs > 0)
                    overlap_found = 0;
                    if (n_segs > 1)
                        for l = 2:n_segs
                            if (baseseg_n(1) >= baseseg_s(l) && baseseg_s(1) ...
                                    <= baseseg_n(l))
                                baseseg_n(1) = max([baseseg_n(1) baseseg_n(l)]);
                                baseseg_s(1) = min([baseseg_s(1) baseseg_s(l)]);
                                overlap_found = 1;
                                if (l == n_segs)
                                    n_segs = n_segs-1;
                                else
                                    for m = l+1:n_segs
                                        baseseg_n(m-1) = baseseg_n(m);
                                        baseseg_s(m-1) = baseseg_s(m);
                                    end;
                                    n_segs = n_segs-1;
                                end;
                                break;
                            end;
                        end;
                    end;
                        
                    if (n_segs == 1)
                        ind_segs = ind_segs+1;
                        indseg_n(ind_segs) = baseseg_n(1);
                        indseg_s(ind_segs) = baseseg_s(1);
                        n_segs = n_segs-1;
                    else
                        if (overlap_found == 0)
                            ind_segs = ind_segs+1;
                            indseg_n(ind_segs) = baseseg_n(1);
                            indseg_s(ind_segs) = baseseg_s(1);
                            for l = 2:n_segs                                    
                                baseseg_n(l-1) = baseseg_n(l);
                                baseseg_s(l-1) = baseseg_s(l);
                            end;
                            n_segs = n_segs-1;
                        end;
                    end;
                end;
      
                for l = 1:ind_segs
                    sy(k,j) = sy(k,j) + indseg_n(l)-indseg_s(l);
                end;
                clear baseseg_n baseseg_s baseseg_bound ;
            end;
        end;
        
    end; %@@@ End of obstruction grid in y                                 

    %@@@ Setting the obstruction grid to zero if neighboring cells are
    %@@@ dry cells to prevent spurious swell attenuation near the coast 

     if (j < Nx && mask(k,j+1) == 0)
        sx(k,j) = 0;
    end;
    if (j > 1 && mask(k,j-1) == 0)
        sx(k,j) = 0;
    end;
    if (k < Ny && mask(k+1,j) == 0)
        sy(k,j) = 0;
    end;
    if (k > 1 && mask(k-1,j) == 0)
        sy(k,j) = 0;
    end;

end; %@@@ End of wet cell (with enclosed boundaries) loop

return;
