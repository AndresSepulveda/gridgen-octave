function bound_ingrid = oct_split_boundary(bound,lim)

% -------------------------------------------------------------------------
%|                                                                        |
%|                    +----------------------------+                      |
%|                    | GRIDGEN          NOAA/NCEP |                      |
%|                    |                            |                      |
%|                    | Last Update :  23-Oct-2012 |                      |
%|                    +----------------------------+                      |
%|                     Distributed with WAVEWATCH III                     |
%|                                                                        |
%|                 Copyright 2009 National Weather Service (NWS),         |
%|  National Oceanic and Atmospheric Administration.  All rights reserved.|
%|                                                                        |
%| DESCRIPTION                                                            |
%| This function splits up large boundary segments into smaller ones so   |
%| that they are more managable                                           |
%|                                                                        |
%| bound_ingrid = oct_split_boundary(bound,lim,[bflg])                        |
%|                                                                        |
%| INPUT                                                                  |
%|  bound : Data structure array of boundary polygons that lie inside the |
%|          grid domain                                                   |
%|  lim   : Limiting size to determine if a polygon needs to be split     |
%|                                                                        |
%| OUTPUT                                                                 |
%| bound_ingrid : A new data structure of boundary polygons where the     |
%|                larger polygons have been split up to more managable    |
%|                smaller sizes                                           |
%|                                                                        |
% -------------------------------------------------------------------------

eps = 1e-5;

N = length(bound);
in_coord = 1;
itmp = 0;

%@@@ The pieces are collected in a cell array and concatenated once at the
%@@@ end. Growing a struct array one element at a time (as in
%@@@ "bound_ingrid(in_coord) = bound(i)") copies the whole array, with all
%@@@ the coordinates it holds, on every single boundary.

parts = cell(1,N);

for i = 1:N  
    if (bound(i).width > lim || bound(i).height > lim)
        low = floor(bound(i).west);
        high = ceil(bound(i).east);
        x_axis = [low:lim:high];
        if x_axis(end) < high
            x_axis(end+1) = high;
        end;
        low = floor(bound(i).south);
        high = ceil(bound(i).north);
        y_axis = [low:lim:high];
        if y_axis(end) < high
            y_axis(end+1) = high;
        end;
        
        Nx = length(x_axis);
        Ny = length(y_axis);

        sub = cell(1,(Nx-1)*(Ny-1));
        ns = 0;

        for lx = 1:Nx-1
            for ly = 1:Ny-1
                lat_start = y_axis(ly);
                lon_start = x_axis(lx);
                lat_end = y_axis(ly+1);
                lon_end = x_axis(lx+1);
                
                [bt,Nb] = oct_compute_boundary([lat_start lon_start lat_end ...
                    lon_end],bound(i),bound(i).level); 
                if (Nb > 0)
                    ns = ns+1;
                    sub{ns} = bt;
                    in_coord = in_coord + Nb;
                end;
                clear bt;
            end;
        end;
        if (ns > 0)
            parts{i} = [sub{1:ns}];
        end;
        clear sub;
    else
        parts{i} = bound(i);
        in_coord = in_coord+1;
    end;
    itmp_prev = itmp;
    itmp = floor(i/N*100);
    if (mod(itmp,5)==0 && itmp_prev ~= itmp && N > 100)
        fprintf(1,'Completed %d per cent of %d boundaries and split into %d boundaries \n',...
            itmp,N,in_coord-1);
    end;
end;   

%@@@ Single concatenation of all the pieces, in the original order

nonempty = ~cellfun('isempty',parts);
if (any(nonempty))
    bound_ingrid = [parts{nonempty}];
else
    bound_ingrid = [];
end;

return;
