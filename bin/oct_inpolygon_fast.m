function [in,on] = oct_inpolygon_fast(x,y,xv,yv)

% -------------------------------------------------------------------------
%| oct_inpolygon_fast                                                         |
%|                                                                        |
%| Same result as the built-in  [in,on] = inpolygon(x,y,xv,yv)  but much  |
%| faster when the polygon has many vertices and the points lie on a      |
%| small number of different latitudes (as the sampling points of the     |
%| land-sea mask clean up do).                                            |
%|                                                                        |
%| The built-in routine tests every point against every edge of the      |
%| polygon. An edge can only matter for a point if the latitude of the    |
%| point is inside the latitude range of the edge, so here the points are |
%| grouped by latitude and each group is only tested against the edges    |
%| that span that latitude. The tests that are applied to those edges are |
%| the standard crossing-number test ("in") and the point-on-edge test    |
%| ("on"), exactly the same expressions the built-in uses; points         |
%| inside the polygon or on its boundary are flagged in "in".             |
%|                                                                        |
%| Falls back to the built-in routine for small problems, where it is     |
%| already fast, or when there are too many different latitudes for the   |
%| grouping to pay off.                                                   |
%| (Implicit array expansion: MATLAB R2016b or later, or Octave.)         |
% -------------------------------------------------------------------------

xv = xv(:);  yv = yv(:);
npol = numel(xv);
npts = numel(x);

if (npol < 500 || npts == 0)
    [in,on] = inpolygon(x,y,xv,yv);
    return;
end;

[ylev,~,grp] = unique(y(:));
nlev = numel(ylev);

if (nlev > 0.5*npts)
    [in,on] = inpolygon(x,y,xv,yv);
    return;
end;

%@@@ Edges (xi,yi) -> (xj,yj) as in the built-in routine: j is the vertex
%@@@ before i (cyclically)

xi = xv;                 yi = yv;
xj = xv([npol 1:npol-1]); yj = yv([npol 1:npol-1]);
dxv = xj - xi;           dyv = yj - yi;
emin = min(yi,yj);       emax = max(yi,yj);
exlo = min(xi,xj);       exhi = max(xi,xj);

in = false(size(x));
on = false(size(x));

%@@@ Group the points by latitude level (grp gives the level of each point)

[~,order] = sort(grp);
cnt = accumarray(grp,1);
last = cumsum(cnt);
first = last-cnt+1;

maxel = 4e6;     % largest (points x edges) block handled at once

for L = 1:nlev

    yl = ylev(L);
    cand = find(emin <= yl & emax >= yl);
    pid = order(first(L):last(L));
    if (isempty(cand))
        continue;
    end;

    xi_c = xi(cand)';  yi_c = yi(cand)';  yj_c = yj(cand)';
    dx_c = dxv(cand)'; dy_c = dyv(cand)';
    xlo_c = exlo(cand)'; xhi_c = exhi(cand)';

    %@@@ conditions that only depend on the latitude of the level
    c1 = ((yi_c <= yl & yl < yj_c) | (yj_c <= yl & yl < yi_c));
    c2 = ((yi_c <= yl & yl <= yj_c) | (yj_c <= yl & yl <= yi_c));

    nc = numel(cand);
    blk = max(1,floor(maxel/nc));

    for b0 = 1:blk:numel(pid)
        pb = pid(b0:min(b0+blk-1,numel(pid)));
        xp = x(pb);  xp = xp(:);

        distance = dx_c .* (yl - yi_c) - (xp - xi_c) .* dy_c;

        cross = c1 & (0 < distance.*dy_c);
        in(pb) = (mod(sum(cross,2),2) == 1);

        onedge = c2 & (xlo_c <= xp & xp <= xhi_c) & (0 == distance | ~dx_c);
        on(pb) = any(onedge,2);
    end;
end;

in = in | on;

return;
