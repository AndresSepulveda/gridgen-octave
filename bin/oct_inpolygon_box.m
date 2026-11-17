function [in,on] = oct_inpolygon_box(x,y,xv,yv)

% -------------------------------------------------------------------------
%| oct_inpolygon_box                                                          |
%|                                                                        |
%| Same result as  [in,on] = inpolygon(x,y,xv,yv)  for the particular case |
%| where the polygon is the closed, axis aligned rectangle                 |
%|                                                                        |
%|     xv = [xa xb xb xa xa] ,  yv = [ya ya yb yb ya]                      |
%|                                                                        |
%| that gridgen uses for grid domains and for the sub-boxes of            |
%| oct_split_boundary. For such a polygon the point-in-polygon test is a pair |
%| of interval tests, so the whole set of points is classified with a few |
%| comparisons instead of a loop over the polygon edges.                  |
%|                                                                        |
%| Any other polygon (or a degenerate rectangle) is passed on to the      |
%| built-in routine, so the caller can use this for every polygon.        |
% -------------------------------------------------------------------------

xv = xv(:);  yv = yv(:);

isbox = false;
if (numel(xv) == 5 && numel(yv) == 5)
    xa = xv(1);  xb = xv(2);  ya = yv(1);  yb = yv(3);
    isbox = isequal(xv,[xa;xb;xb;xa;xa]) && isequal(yv,[ya;ya;yb;yb;ya]) ...
            && xb > xa && yb > ya;
end;

if (~isbox)
    [in,on] = inpolygon(x,y,xv,yv);
    return;
end;

%@@@ Strictly inside, on one of the four edges (corners included)

inx = (x > xa) & (x < xb);
iny = (y > ya) & (y < yb);

onx = ((x == xa) | (x == xb)) & (y >= ya) & (y <= yb);
ony = ((y == ya) | (y == yb)) & (x >= xa) & (x <= xb);

on = onx | ony;
in = (inx & iny) | on;

return;
