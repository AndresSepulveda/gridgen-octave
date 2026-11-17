 function depth_sub = oct_generate_grid(x,y,ref_dir,bathy_source,limit,cut_off,dry,varargin)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  27-Jan-2017 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| This function creates a 2D bathymetry data set from high resolution    | 
%| "ETOPO1" or "ETOPO2" global bathymetry sets. Global bathymetry data    |
%| sets are assumed to be stored in Netcdf formats                        |
%|                                                                        |
%| depth = oct_generate_grid(x,y,ref_dir,bathy_source,limit,cut_off,dry,[var])|
%|                                                                        |
%| INPUT                                                                  |
%|  x            : A 2D array specifying the longitudes of each cell      | 
%|  y            : A 2D array specifying the lattitudes of each cell      |
%|  ref_dir      : PATH string to where the global reference bathymetry   |
%|                    data sets are stored                                |
%|  bathy_source : String file to indicate which type of bathymetry is    |
%|                 being used. User needs to make sure that the bathymetry|
%|                 data files corresponding to the options are available  |
%|                 in ref_dir                                             |
%|                 Options are --                                         |
%|                     'etopo1' -- ETOPO1 bathymetry (etopo1.nc)          |
%|                     'etopo2' -- ETOPO2 bathymetry (etopo2.nc)          |
%|  limit        : Value ranging between 0 and 1 indicating what fraction |
%|                 of a grid cell needs to be covered by wet cells (from  |
%|                 the base grid) for the cell to be marked wet           | 
%|  cut_off      : Cut_off depth to distinguish between dry and wet cells.|
%|                 All depths below the cut_off depth are marked wet      | 
%|  dry          : Depth value assigned to the dry cells                  |
%|  var          : These are optional string arrrays for variable         |
%|                 definition names for lon (x), lat (y) and depth        |
%|                 respectively. If ommitted default names are used.      | 
%|                   For the etopo2.nc file these are 'x',   'y'   and 'z'|
%|                   For the etopo1.nc file these are 'lon', 'lat' and 'z'|
%|                                                                        |
%| OUTPUT                                                                 |
%|  depth        : A 2D array of dimensions (Nx,Ny) consisting of the grid|
%|                 depths                                                 |
%|                                                                        |
%| NOTES                                                                  |
%|    a. This version uses the in-house matlab NETCDF functions which is  |
%|       available with Matlab 2008a or higher. If you are using an older |
%|       version of MATLAB then you will have to install a NETCDF package |
%|       and change portions of this script that handle netcdf files.     |  
%|       Keep in mind that the NETCDF functions used in this package      |
%|       start their index from 0, This may not be the case in other      |
%|       NETCDF packages                                                  |
%|    b. The default bathymetric sets are etopo1 and etopo2, and while    |
%|       this package will allow you to create grids finer than these     |
%|       base grids, make sure that features are defined in the base      |
%|       bathymetries before using them.                                  |
%|    c. While the code allows for wrapping around in the base grid it    |
%|       assumes that the target grid will be monotonically increasing.   |
%|    d. The latest version now works with curvilinear grids. The function|
%|       no longer assumes constant increments in x (lon) and y (lat)     |
%|       directions. It generates the depth values based on user provided |
%|       lat/lon arrays. Locally the depths are averaged (interpolated)   |
%|       from base grid depending on whether the cell dimnensions are     |
%|       larger (or smaller) than the base grid cells                     |
%|	e. This version is compatible with Octave + the netcdf package 	  |
%|        which has to be installed, for more information see here:	  |
%|		 <https://octave.sourceforge.io/netcdf/>				  |
%|                                                                        |
%| BUG FIXES                                                              |
%|    06/25/2012 : Fixed erroneous generation of NaNs in grids that wrap  |
%|                 around the globe                                       |
%|    10/23/2012 : Removed coord choice (-180 to 180 or 0 360). Now grids |
%|                 are defined only from 0 to 360 (to avoid spurious      |
%|                 boundaries near the edges)                             |
%|    05/02/2013 : Modified the algorithm for curvilinear grids.          |
%|                 (See Notes)							  |
%|	  01/27/2017 : Octave compatible                                      |
%|    06/27/2017 : Modified for very fine grids, Deanna's update          |
% -------------------------------------------------------------------------
%@@@ Octave Friendly
vers=ver;
for i1=1:1:length(vers)
    if strcmpi (vers(i1).Name, 'Octave')
        pkg load netcdf;
        import_netcdf;
    end
end

%@@@ Read additional input variables if they exist (read arbitrary base depth)
narg=nargin;

%@@@ Determine if third input variable present (requesting inland features)
if narg == 10 %@@@ Extra 3 arguments define the lat, lon and depth var names
    var_x = varargin{1};
    var_y = varargin{2};
    var_z = varargin{3};
    bathy_custom = bathy_source;    
elseif narg == 7
    bathy_custom='none';
elseif narg < 7
  disp('Too few input arguments, exiting')
else
  disp('Too many input arguments, exiting')
end

%@@@ Initialize the corners of the grid domain and the depth values

lats = min(min(y));
lons = min(min(x));
late = max(max(y));
lone = max(max(x));

depth_sub = zeros(size(x));

%@@@ Compute cell corners (vectorised: one call for the whole grid instead
%@@@ of one call per cell and a Nx*Ny struct array)

[Ny,Nx] = size(x);

if (Nx >= 3 && Ny >= 3)
    [c1x,c1y,c2x,c2y,c3x,c3y,c4x,c4y,Wc,Hc] = oct_compute_cellcorner_all(x,y);
else
    c1x = zeros(Ny,Nx); c1y = c1x; c2x = c1x; c2y = c1x;
    c3x = c1x; c3y = c1x; c4x = c1x; c4y = c1x; Wc = c1x; Hc = c1x;
    for j = 1:Nx
        for k = 1:Ny
            [c1,c2,c3,c4,wdth,hgt] = oct_compute_cellcorner(x,y,j,k,Nx,Ny);
            c1x(k,j)=c1(1); c1y(k,j)=c1(2); c2x(k,j)=c2(1); c2y(k,j)=c2(2);
            c3x(k,j)=c3(1); c3y(k,j)=c3(2); c4x(k,j)=c4(1); c4y(k,j)=c4(2);
            Wc(k,j)=wdth; Hc(k,j)=hgt;
        end;
    end;
end;

dx = max(Wc(:));
dy = max(Hc(:));

%@@@ Determine the file name for source bathymetry

if (strcmp(bathy_source,'etopo2'))
    fname_base = [ref_dir,'/etopo2.nc'];
    var_x = 'x';
    var_y = 'y';
    var_z = 'z';    
elseif (strcmp(bathy_source,'etopo1'))
    fname_base = [ref_dir,'/etopo1.nc'];
    var_x = 'lon';
    var_y = 'lat';
    var_z = 'z';    
elseif (strcmp(bathy_source,bathy_custom))
    fname_base = [ref_dir,'/',bathy_custom,'.nc'];
else
    fprintf(1,'Unrecognized source bathymetry option \n');
    return;
end;

%@@@ Determine dimensions and ranges of base bathymetry coords

f = netcdf.open(fname_base,'nowrite');

varid_lon = netcdf.inqVarID(f,var_x);
varid_lat = netcdf.inqVarID(f,var_y);
varid_dep = netcdf.inqVarID(f,var_z);

[~,Nx_base]=netcdf.inqDim(f,varid_lon);
[~,Ny_base]=netcdf.inqDim(f,varid_lat);

lat_range= netcdf.getAtt(f,varid_lat,'actual_range');
dy_base = diff(lat_range)/(Ny_base-1);

lats_base=lat_range(1);
late_base=lat_range(2);

lon_range= netcdf.getAtt(f,varid_lon,'actual_range');
dx_base = diff(lon_range)/(Nx_base-1);

lons_base=lon_range(1);
lone_base=lon_range(2); 

if (lats < lats_base || lats > late_base || late < lats_base || late > ...
        late_base)
    fprintf(1,'ERROR : Lattitudes (%d,%d) beyond range (%d,%d) \n',lats,...
        late,lats_base,late_base);
    return;
end;
  
if (lons < lons_base || lons > lone_base || lone < lons_base || lone > ...
        lone_base)
    fprintf(1,'ERROR : Longitudes (%d,%d) beyond range (%d,%d) \n',lons,...
        lone,lons_base,lone_base);
    return;
end;

%@@@ Determine the starting and end points for extracting lattitude data
%@@@ from NETCDF

lat_start = floor(( (lats-2*dy) - lats_base)/dy_base);

if (lat_start < 1)
    lat_start = 1;
end;
 
lat_end = ceil(((late+2*dy) - lats_base)/dy_base) +1;

if (lat_end > Ny_base)
    lat_end = Ny_base;
end;

%@@@ Determine the starting and end points for extracting longitude data 
%@@@ from NETCDF

lon_start = floor(((lons-2*dx) - lons_base)/dx_base);
lon_end = ceil(((lone+2*dx) - lons_base)/dx_base) +1;

if (lon_start < 1)
    lon_start = 1;
end;

if (lon_start > Nx_base)
    lon_start = Nx_base;
end;

if (lon_end < 1)
    lon_end = 1;
end;

if (lon_end >Nx_base)
    lon_end = Nx_base;
end;    

%@@@ Extract data from Netcdf files
%@@@ The next few lines assume that your Matlab
%@@@ version has ability to read NETCDF files
%@@@ !! NOTE: Indexing from 0 for these functions!!


count_lat = (lat_end - lat_start) + 1;
lat_base = netcdf.getVar(f,varid_lat,lat_start-1,count_lat);             
                                                                  
if (lon_end <= lon_start)
    count_lon1 = (Nx_base - lon_start) + 1;                       
    count_lon2 = (lon_end - 2) + 1;                               
    lon1 = netcdf.getVar(f,varid_lon,lon_start-1,count_lon1);    
    lon2 = netcdf.getVar(f,varid_lon,1,count_lon2);               
    lon_base = [lon1;lon2];
    dep1 = ( netcdf.getVar(f,varid_dep,...
          [lon_start-1 lat_start-1],[count_lon1 count_lat]) )';   
    dep2 = ( netcdf.getVar(f,varid_dep,...
          [1 lat_start-1],[count_lon2 count_lat]) )';             
    depth_base = [dep1 dep2];
else
    count_lon = (lon_end - lon_start) + 1;                        
    lon_base = netcdf.getVar(f,varid_lon,lon_start-1,count_lon);  
    depth_base = ( netcdf.getVar(f,varid_dep,...
            [lon_start-1 lat_start-1],[count_lon count_lat]) )';  
end;

fprintf(1,'read in the base bathymetry \n');
netcdf.close(f);

%loc = find(lon_base < 0);
%lon_base(loc) = lon_base(loc)+360;
%if (lon_base(end) == 0)
%    lon_base(end) = lon_base(end)+360;
%end;

%Remove overlapped regions (occurs when longitudes wrap around)
%(skipped when lon_base is already strictly increasing: unique/intersect
% would return the identical array, but copying the whole bathymetry
% matrix is expensive)

if (~all(diff(lon_base) > 0))
    clear lon_base_tmp depth_base_tmp ib;
    lon_base_tmp = unique(lon_base);
    [~,~,ib] = intersect(lon_base_tmp,lon_base);
    depth_base_tmp = depth_base(:,ib);

    clear lon_base depth_base;
    lon_base = lon_base_tmp;
    depth_base = depth_base_tmp;
    clear depth_base_tmp;
end;

%@@@ Obtaining data from base bathymetry. If desired grid is coarser than
%@@@ base grid then 2D averaging of bathymetry else grid is interpolated
%@@@ from base grid. Checks if grid cells wrap around in Longitudes. Does
%@@@ not do so for Lattitudes

ndx_c = round(Wc/dx_base);
ndy_c = round(Hc/dy_base);
interp_cell = (ndx_c <= 1 & ndy_c <= 1);

%@@@ ---------------------------------------------------------------------
%@@@ FAST PATH. Cells that are plain rectangles on a rectilinear grid, that
%@@@ do not wrap around in longitude and that need cell averaging are all
%@@@ treated at once: the number of base points, the number of wet base
%@@@ points and their sum inside every cell are obtained from running sums
%@@@ along the base longitudes (one base-latitude strip per row of cells).
%@@@ That gives the same sets of points that the inpolygon test of the
%@@@ cell-by-cell algorithm selects (points on the cell edge included).
%@@@ All remaining cells go through the original cell-by-cell code below.
%@@@ ---------------------------------------------------------------------

fast_cell = false(Ny,Nx);
tolb = 1e-10;     % cells whose edges agree to this level are the same box
tolt = 1e-6;      % a base node closer than this to a cell edge is a "tie"

if (Nx >= 3 && Ny >= 3 && any(~interp_cell(:)) && ...
        all(diff(lat_base) > 0) && all(diff(lon_base) > 0) && ...
        isequal(x,repmat(x(1,:),Ny,1)) && isequal(y,repmat(y(:,1),1,Nx)))

    box_cell = (c1x==c2x) & (c3x==c4x) & (c1y==c4y) & (c2y==c3y) & ...
               (c1x>c4x) & (c2y>c1y);
    wrap_cell = (min(min(c4x,c3x),min(c1x,c2x)) < 0) | ...
                (max(max(c4x,c3x),max(c1x,c2x)) > 360);

    %@@@ Representative (interior) cell edges for every column / row

    xlo = c4x(2,:);  xhi = c1x(2,:);
    ylo = c4y(:,2);  yhi = c2y(:,2);

    same_box = (abs(c4x-repmat(xlo,Ny,1)) <= tolb) & ...
               (abs(c1x-repmat(xhi,Ny,1)) <= tolb) & ...
               (abs(c4y-repmat(ylo,1,Nx)) <= tolb) & ...
               (abs(c2y-repmat(yhi,1,Nx)) <= tolb);

    %@@@ Ties: when a base-grid node lies (numerically) ON a cell edge, whether
    %@@@ it belongs to the cell is decided by round-off in inpolygon (and by
    %@@@ its tolerance, which depends on the MATLAB/Octave version). Such
    %@@@ columns / rows are left to the original cell-by-cell code so the
    %@@@ result stays exactly what inpolygon gives.

    tie_col = false(1,Nx);
    for j = 1:Nx
        tie_col(j) = near_node(lon_base,xlo(j),tolt) || near_node(lon_base,xhi(j),tolt);
    end;
    tie_row = false(Ny,1);
    for k = 1:Ny
        tie_row(k) = near_node(lat_base,ylo(k),tolt) || near_node(lat_base,yhi(k),tolt);
    end;

    fast_cell = ~interp_cell & box_cell & ~wrap_cell & same_box & ...
                ~repmat(tie_col,Ny,1) & ~repmat(tie_row,1,Nx);
    clear box_cell wrap_cell same_box;
end;

if (any(fast_cell(:)))

    fprintf(1,'Generating grid bathymetry (vectorised cell averaging) ....\n');

    %@@@ Range of base-grid longitudes covered by every column of cells
    %@@@ (cell edges included)

    col_i1 = zeros(1,Nx);  col_i2 = zeros(1,Nx);
    for j = 1:Nx
        a = find(lon_base >= xlo(j)-tolb,1,'first');
        b = find(lon_base <= xhi(j)+tolb,1,'last');
        if (isempty(a) || isempty(b) || b < a)
            col_i1(j) = 1;  col_i2(j) = 0;
        else
            col_i1(j) = a;  col_i2(j) = b;
        end;
    end;
    ncol_pts = col_i2 - col_i1 + 1;

    itmp_prev = 0;
    for k = 1:Ny
        cols = find(fast_cell(k,:));
        if (isempty(cols))
            continue;
        end;
        a = find(lat_base >= ylo(k)-tolb,1,'first');
        b = find(lat_base <= yhi(k)+tolb,1,'last');
        if (isempty(a) || isempty(b) || b < a)
            depth_sub(k,cols) = dry;
            continue;
        end;
        nrow_pts = b - a + 1;

        %@@@ Running sums along longitude of (number of wet points) and
        %@@@ (sum of wet depths) for the strip of base latitudes of this row

        dstrip = double(depth_base(a:b,:));
        wet = (dstrip <= cut_off);
        dstrip(~wet) = 0;
        cc = [0 cumsum(sum(wet,1))];
        cs = [0 cumsum(sum(dstrip,1))];

        i1 = col_i1(cols);  i2 = col_i2(cols);
        nwet = cc(i2+1) - cc(i1);
        swet = cs(i2+1) - cs(i1);
        Nt = nrow_pts*ncol_pts(cols);

        isw = (nwet > 0) & (nwet./Nt > limit);
        val = repmat(dry,1,numel(cols));
        val(isw) = swet(isw)./nwet(isw);
        depth_sub(k,cols) = val;

        itmp = floor(k/Ny*100);
        if (mod(itmp,10) == 0 && itmp_prev ~= itmp)
            fprintf(1,'Completed %d per cent of the cell rows \n',itmp);
        end;
        itmp_prev = itmp;
    end;
    clear dstrip wet cc cs;
end;

%@@@ Remaining cells: original cell-by-cell algorithm (interpolation of the
%@@@ base grid, wrap-around cells, curvilinear grids, ...)

slow_list = find(~fast_cell)';
Nslow = numel(slow_list);

if (Nslow > 0)

    fprintf(1,'Generating grid bathymetry for %d cells ....\n',Nslow);
    itmp_prev = 0;

    for ic = 1:Nslow

        lin = slow_list(ic);
        k = mod(lin-1,Ny)+1;
        j = (lin-k)/Ny+1;

        cpx = [c4x(k,j) c1x(k,j) c2x(k,j) c3x(k,j) c4x(k,j)]';
        cpy = [c4y(k,j) c1y(k,j) c2y(k,j) c3y(k,j) c4y(k,j)]';
        cwidth  = Wc(k,j);
        cheight = Hc(k,j);

        ndx = round(cwidth/dx_base);
        ndy = round(cheight/dy_base);

        if (ndx <= 1 && ndy <= 1)

        %@@@ Interpolating from base grid

            den = dx_base*dy_base;
%            [lon_prev,~] = min(abs(lon_base-x(k,j)));
%For very fine grids
            [lon_prev,tmp] = min(abs(lon_base-x(k,j)));
            if ( round(lon_prev) == 0 )
              lon_prev = tmp-1;
            end;
%
            if ( lon_base(lon_prev) > x(k,j) )
                lon_prev = lon_prev - 1;
            end;
            lon_next = lon_prev + 1;
            if (lon_prev == 0 )
                lon_prev = Nx_base - 1;
            end;
            if (lon_next > Nx_base)
                lon_next = 2;
            end;
%            [lat_prev,~] = min(abs(lat_base-y(k,j)));
%For very fine grids
             [lat_prev,tmp] = min(abs(lat_base-y(k,j)));
             if ( round(lat_prev) == 0 )
               lat_prev = tmp-1;
             end;
%
            if ( lat_base(lat_prev) > y(k,j) )
                lat_prev = lat_prev - 1;
            end;
            lat_next = lat_prev + 1;

            dx1 = double(min(abs(x(k,j)-lon_base(lon_prev)),abs(x(k,j)-...
                lon_base(lon_prev)-360*sign(x(k,j)-lon_base(lon_prev)))));
            dx2 = double(dx_base - dx1);
            dy1 = double(y(k,j) - lat_base(lat_prev));
            dy2 = double(dy_base - dy1);

            %@@@ Four point interpolation

            a11 = double(depth_base(lat_prev,lon_prev));
            a12 = double(depth_base(lat_prev,lon_next));
            a21 = double(depth_base(lat_next,lon_prev));
            a22 = double(depth_base(lat_next,lon_next));
            depth_sub(k,j) = (a11*dy2*dx2 + a12*dy2*dx1 + a21*dy1*dx2 + ...
                a22*dx1*dy1)/den;

        else

            %@@@ Cell averaging
            %@@@ Determine the base bathymetry region that covers the cell

            lon_start = min(cpx);
            lon_end = max(cpx);
            lat_start = min(cpy);
            lat_end = max(cpy);

            overlap = 0;
            if (lon_end > 360)
                overlap = 1;
                lon_end = lon_end - 360;
            end;
            if (lon_start < 0)
                overlap = 2;
                lon_start = lon_start + 360;
            end;

            clear depth_tmp lon_tmp lat_tmp lon1 lon2 depth_tmp1 depth_tmp2;

            [~,lat_start_pos] = min(abs(lat_base-lat_start));
            [~,lat_end_pos] = min(abs(lat_base-lat_end));
            lat_tmp = lat_base(lat_start_pos:lat_end_pos);
            [~,lon_start_pos] = min(abs(lon_base-lon_start));
            [~,lon_end_pos] = min(abs(lon_base-lon_end));

            switch overlap
                case 0
                    depth_tmp = depth_base(lat_start_pos:lat_end_pos,...
                        lon_start_pos:lon_end_pos);
                    lon_tmp = lon_base(lon_start_pos:lon_end_pos);
                case 1
                    depth_tmp1 = depth_base(lat_start_pos:lat_end_pos,...
                        lon_start_pos:end);
                    depth_tmp2 = depth_base(lat_start_pos:lat_end_pos,...
                        2:lon_end_pos);
                    lon1 = lon_base(lon_start_pos:end);
                    lon2 = lon_base(2:lon_end_pos)+360.0;
                    lon_tmp = [lon1;lon2];
                    depth_tmp = [depth_tmp1 depth_tmp2];
                case 2
                    depth_tmp1 = depth_base(lat_start_pos:lat_end_pos,...
                       lon_start_pos:end);
                    depth_tmp2 = depth_base(lat_start_pos:lat_end_pos,...
                        2:lon_end_pos);
                    lon1 = lon_base(lon_start_pos:end) - 360;
                    lon2 = lon_base(2:lon_end_pos);
                    lon_tmp = [lon1;lon2];
                    depth_tmp = [depth_tmp1 depth_tmp2];
            end;

            %@@@ Compute the average depth from points that lie inside the
            %@@@ the cell and are below the cut off


            clear lon_tmp2d lat_tmp2d in_cell depth_tmp_incell;
            [lon_tmp2d,lat_tmp2d] = meshgrid(lon_tmp,lat_tmp);
            in_cell = inpolygon(lon_tmp2d,lat_tmp2d,cpx,cpy);
            depth_tmp_incell = depth_tmp(in_cell > 0);
            Nt = numel(depth_tmp_incell);
            clear loc;
            loc = find(depth_tmp_incell <= cut_off);
            if (~isempty(loc))
                Ntt = length(loc);
                if (Ntt/Nt > limit)
                    depth_sub(k,j) = mean(depth_tmp_incell(loc));
                else
                    depth_sub(k,j) = dry;
                end;
            else
                depth_sub(k,j) = dry;
            end;

        end %@@@ End of check to see if it will interpolate or average

        %@@@ Counter to check proportion of cells completed

        itmp = floor(ic/Nslow*100);
        if (mod(itmp,5) == 0 && itmp_prev ~= itmp)
            fprintf(1,'Completed %d per cent of the remaining cells \n',itmp);
        end;
        itmp_prev = itmp;

    end; %@@@ Loop through the remaining cells

end;

return;

function tf = near_node(v,val,tol)
%@@@ true if a node of the sorted vector v lies within tol of val
n = numel(v);
a = find(v >= val,1,'first');
if (isempty(a))
    tf = abs(v(n)-val) <= tol;
elseif (a == 1)
    tf = abs(v(1)-val) <= tol;
else
    tf = (abs(v(a)-val) <= tol) || (abs(v(a-1)-val) <= tol);
end;
