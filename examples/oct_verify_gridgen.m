function oct_verify_gridgen(orig_dir, opt_dir, ref_dir, box, dxdeg, varargin)
% OCT_VERIFY_GRIDGEN  Compare the ORIGINAL gridgen against this oct_ version.
%
%   oct_verify_gridgen(orig_dir, opt_dir, ref_dir, box, dxdeg)
%   oct_verify_gridgen(orig_dir, opt_dir, ref_dir, box, dxdeg, 'boundary','full')
%
%   orig_dir : folder of the original gridgen   (contains bin/ with the
%              functions under their original names)
%   opt_dir  : folder of this version           (contains bin/ with the
%              functions under their oct_ names)
%   ref_dir  : folder with the reference data (etopo1.nc, coastal_bound_*.mat,
%              optional_coastal_polygons.mat, user_polygons.flag)
%   box      : [lat_min lon_min lat_max lon_max]  (lon 0..360) of a SMALL
%              test region, e.g. [-40 270 -30 285] -- the original code is
%              very slow, so keep the region small (a few thousand cells)
%   dxdeg    : grid resolution in degrees (dx = dy), e.g. 0.25
%
% It runs the same steps as examples/create_glo_15m.m (bathymetry, boundaries,
% polygon splitting, mask clean up, lake removal, obstructions, and the
% WAVEWATCH III text files) with both code bases, prints the time of every
% step and reports whether every result is IDENTICAL.
%
% Because the two code bases no longer share function names, both bin folders
% can sit on the path at the same time: the original is called as
% generate_grid, clean_mask, ... and this version as oct_generate_grid,
% oct_clean_mask, ...  The dispatch is done with str2func on a prefix.
%
% Note: for global grids (lon 0..360) the original code is far too slow to
% run; use a regional box for the comparison, and use
% examples/bench_create_glo_15m.m to time this version on the full case.

p = inputParser;
addParameter(p,'boundary','full');
addParameter(p,'maxpts',Inf);
parse(p,varargin{:});
boundary = p.Results.boundary;

res = cell(1,2);
names    = {'original','oct_ version'};
dirs     = {orig_dir,opt_dir};
prefixes = {'','oct_'};

for ic = 1:2
    fprintf('\n=========== %s ===========\n',names{ic});
    addpath(fullfile(dirs{ic},'bin'));
    res{ic} = run_case(prefixes{ic},ref_dir,box,dxdeg,boundary,...
        strrep(names{ic},' ','_'));
    rmpath(fullfile(dirs{ic},'bin'));
end

fprintf('\n=========== comparison ===========\n');
a = res{1};  b = res{2};
items = {'depth','m2','m3','m4','mask_map','sx1','sy1'};
allok = true;
for i = 1:numel(items)
    ok = isequaln(a.(items{i}),b.(items{i}));
    allok = allok && ok;
    mx = max(abs(a.(items{i})(:)-b.(items{i})(:)));
    fprintf('%-10s identical: %d   max abs difference: %g\n',items{i},ok,mx);
end
files = {'depth_ascii','maskorig_ascii','obstr_lev1'};
for i = 1:numel(files)
    ta = fileread(a.files.(files{i}));  tb = fileread(b.files.(files{i}));
    ok = strcmp(ta,tb);  allok = allok && ok;
    fprintf('file %-16s identical: %d\n',files{i},ok);
end
fprintf('\n%-22s %12s %12s %9s\n','step','original (s)','oct_ (s)','speed-up');
f = fieldnames(a.T);
for i = 1:numel(f)
    fprintf('%-22s %12.2f %12.2f %8.1fx\n',f{i},a.T.(f{i}),b.T.(f{i}),a.T.(f{i})/max(b.T.(f{i}),1e-3));
end
ta = sum(cellfun(@(s) a.T.(s),f));  tb = sum(cellfun(@(s) b.T.(s),f));
fprintf('%-22s %12.2f %12.2f %8.1fx\n','TOTAL',ta,tb,ta/tb);
if allok
    fprintf('\nALL RESULTS IDENTICAL\n');
else
    fprintf('\nWARNING: differences found (see above)\n');
end
end


function R = run_case(pre,ref_dir,box,dxdeg,boundary,tag)
% Same flow as examples/create_glo_15m.m, for a single box. "pre" is the
% prefix of the code base being exercised: '' for the original, 'oct_' for
% this version.

f_optional_bound   = str2func([pre 'optional_bound']);
f_generate_grid    = str2func([pre 'generate_grid']);
f_compute_boundary = str2func([pre 'compute_boundary']);
f_split_boundary   = str2func([pre 'split_boundary']);
f_clean_mask       = str2func([pre 'clean_mask']);
f_remove_lake      = str2func([pre 'remove_lake']);
f_create_obstr     = str2func([pre 'create_obstr']);
f_write_ww3file    = str2func([pre 'write_ww3file']);
f_write_ww3obstr   = str2func([pre 'write_ww3obstr']);

T = struct();
S = load(fullfile(ref_dir,['coastal_bound_',boundary,'.mat']));  bound = S.bound;
[bound_user,Nu] = f_optional_bound(ref_dir,fullfile(ref_dir,'user_polygons.flag'));

y = box(1):dxdeg:box(3);  x = box(2):dxdeg:box(4);
[X,Y] = meshgrid(x,y);  lon = x;  lat = y;

t0 = tic;
depth = f_generate_grid(X,Y,ref_dir,'etopo1',.5,-1,0);        T.bathymetry = toc(t0);

t0 = tic;
[b,N1] = f_compute_boundary(box,bound);
if (Nu > 0), [b_opt,N2] = f_compute_boundary(box,bound_user); else, N2 = 0; end
T.compute_boundary = toc(t0);

m = ones(size(depth));  m(depth == 999) = 0;
t0 = tic;  b_split = f_split_boundary(b,10);                  T.split_boundary = toc(t0);

t0 = tic;
m2 = f_clean_mask(X,Y,m,b_split,0.5,0.5);
if (Nu > 0 && N2 ~= 0), m3 = f_clean_mask(X,Y,m2,b_opt,0.5,0.5); else, m3 = m2; end
T.clean_mask = toc(t0);

t0 = tic;  [m4,mask_map] = f_remove_lake(m3,-1,1);           T.remove_lake = toc(t0);
t0 = tic;  [sx1,sy1] = f_create_obstr(X,Y,b,m4,1,1);         T.create_obstr = toc(t0);

t0 = tic;
out_dir = fullfile(tempdir,['gridgen_verify_',tag]);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
R.files.depth_ascii    = fullfile(out_dir,'t.depth_ascii');
R.files.maskorig_ascii = fullfile(out_dir,'t.maskorig_ascii');
R.files.obstr_lev1     = fullfile(out_dir,'t.obstr_lev1');
f_write_ww3file(R.files.depth_ascii,round(depth*1000));
f_write_ww3file(R.files.maskorig_ascii,m4);
f_write_ww3obstr(R.files.obstr_lev1,round(sx1*100),round(sy1*100));
T.write_files = toc(t0);

R.depth = depth;  R.m2 = m2;  R.m3 = m3;  R.m4 = m4;  R.mask_map = mask_map;
R.sx1 = sx1;  R.sy1 = sy1;  R.T = T;
end
