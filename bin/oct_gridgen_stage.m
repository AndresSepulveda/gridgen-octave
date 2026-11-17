function oct_gridgen_stage(name)
% GRIDGEN_STAGE  Simple stage timer used by the bench_* example scripts.
%   oct_gridgen_stage('<stage name>')  marks the beginning of a stage
%   oct_gridgen_stage('report')        prints the time taken by each stage
persistent names times t0
if isempty(names), names = {}; times = []; t0 = tic; end
if strcmp(name,'report')
    fprintf(1,'\n==== gridgen stage timing (s) ====\n');
    d = diff(times);
    for i = 1:numel(d), fprintf(1,'%-30s %10.1f\n', names{i}, d(i)); end
    fprintf(1,'%-30s %10.1f\n','TOTAL',times(end)-times(1));
    names = {}; times = [];
    return;
end
names{end+1} = name; times(end+1) = toc(t0);
