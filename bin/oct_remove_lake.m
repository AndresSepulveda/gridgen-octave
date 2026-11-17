function [mask_mod,mask_map] = oct_remove_lake(mask,lake_tol,igl)

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
%| This routine groups wet cells into independant water bodies with all   | 
%| the wet cells connected to each other sharing the same unique ID       |
%|                                                                        |
%|  [mask_mod,mask_map] = oct_remove_lake(mask,lake_tol,igl)                  | 
%|                                                                        |
%| INPUT                                                                  |
%|  mask     : Input 2D land/sea mask                                     |
%|  lake_tol : Tolerance value that determines all the wet cells          | 
%|             corresponding to a particular wet body should be flagged   | 
%|             dry or not. If the value is a +ve number then all water    | 
%|             bodies having less than this value of total wet cells will |
%|             be flagged dry. If the value is 0 then the output and input|
%|             masks are unchanged. If the value is -ve then all but the  |
%|             largest water body is flagged dry                          |
%|  igl      : Switch to determine if the grid is global or regional. This| 
%|             is needed to determine if the first and last cells along a | 
%|             particular row (longitude/x) are connected or not.         |
%|             Options are --                                             |
%|                igl = 0 for regional (not connected) grids              |
%|                igl = 1 for global (connected) grids                    |
%|                                                                        |
%| OUTPUT                                                                 |
%|  mask_mod : Modified 2D land/sea mask based on the value of lake_tol   |
%|  mask_map : 2D array that has a value of -1 for all land (dry) cells   |
%|             and unique IDs for wet cells that are part of a water body.|
%|                                                                        |
% -------------------------------------------------------------------------

[Ny,Nx] = size(mask);

%@@@ Fast path: label the 4-connected water bodies with a sparse-matrix
%@@@ connected-components search instead of a cell-by-cell flood fill. The
%@@@ labels are numbered exactly as the flood fill numbered them (in order
%@@@ of the first wet cell of each water body in column-major order), so
%@@@ the outputs are identical. Only valid for 0/1 masks with at least one
%@@@ wet cell; anything else uses the original routine.

wet = (mask == 1);
nw = nnz(wet);

if (~(all(mask(:) == 0 | mask(:) == 1) && nw > 0))
    [mask_mod,mask_map] = oct_remove_lake_orig(mask,lake_tol,igl);
    return;
end;

%@@@ Wet cells numbered in column-major order

idx = zeros(Ny,Nx);
idx(wet) = 1:nw;

%@@@ Links between 4-connected wet neighbours (wrap around in x if igl == 1)

mh = wet(:,1:Nx-1) & wet(:,2:Nx);
mv = wet(1:Ny-1,:) & wet(2:Ny,:);
ia = idx(:,1:Nx-1);   ib = idx(:,2:Nx);
ja = idx(1:Ny-1,:);   jb = idx(2:Ny,:);
ea = [ia(mh); ja(mv)];
eb = [ib(mh); jb(mv)];

if (igl == 1 && Nx > 2)
    mw = wet(:,Nx) & wet(:,1);
    ea = [ea; idx(mw,Nx)];
    eb = [eb; idx(mw,1)];
end;

A = sparse([ea; eb; (1:nw)'],[eb; ea; (1:nw)'],1,nw,nw);

%@@@ Connected components (blocks of the Dulmage-Mendelsohn decomposition)

[p,~,r] = dmperm(A);
blk = zeros(nw,1);
blk(r(1:end-1)) = 1;
blk = cumsum(blk);
comp = zeros(nw,1);
comp(p) = blk;

%@@@ Number the water bodies by the order of their first wet cell

[~,first] = unique(comp,'first');
[~,ord] = sort(first);
rnk = zeros(numel(ord),1);
rnk(ord) = 1:numel(ord);
lab = rnk(comp);

last_mask = numel(ord);
N1 = accumarray(lab,1)';

mask_map = mask-1;
mask_map(wet) = lab;

for i = 1:last_mask
    fprintf(1,'%d Wet cells set to flag id %d \n',N1(i),i);
end;

%@@@ Modify mask based on lake_tol value

mask_mod = mask;

if (lake_tol < 0)                    
    [~,pos] = max(N1);
    for i = 1:last_mask
        if (i ~= pos)
            fprintf(1,'Masking out cells with flag set to %d\n',i);
        end;
    end;
    mask_mod(mask_map ~= pos) = 0;
else                                 
    small = find(N1 < lake_tol);
    for i = small
        fprintf(1,'Masking out cells with flag set to %d\n',i);
    end;
    mask_mod(ismember(mask_map,small)) = 0;
end;

return;
