function out = aqua_export_events(matPath)
%AQUA_EXPORT_EVENTS  Export the events of ANY AQuA2 result file (v37).
%   aqua_export_events('D:\...\<run>_bin31f_1p98s\<run>_bin31f_1p98s_AQuA2.mat')
%
%   Writes next to the .mat, exactly as the pipeline's CFU batch does for its
%   own results: <name>_events.csv (one row per event: frames, centre, footprint
%   size and box, peak dF/F, area), <name>_eventmap.tif and <name>_meanimg.tif.
%   Meant for a result saved BY HAND from the AQuA2 GUI (its default save
%   folder is a subfolder named after the movie): the tuning report then uses
%   it as the REFERENCE for that run — per swept setting, how many of the
%   hand-accepted events each setting recovers and how many others it adds
%   (python aqua_tune.py <mouse> --report; GUI tab 3 "3 Report").
%   The AQuA2 folder must be on the MATLAB path (startup of AQuA2), as for the
%   batch scripts. Returns the events csv path.
if nargin < 1 || isempty(matPath)
    [f, p] = uigetfile('*_AQuA2.mat', 'AQuA2 result to export');
    if isequal(f, 0), out = ''; return; end
    matPath = fullfile(p, f);
end
if ~exist(matPath, 'file'), error('not found: %s', matPath); end
[filepath, f1, ~] = fileparts(matPath);
base = fullfile(filepath, f1);
fprintf('loading %s ...\n', matPath);
S = load(matPath, 'res');
res = S.res;
H = size(res.datOrg1, 1); W = size(res.datOrg1, 2);
imwrite(uint16(round(mean(single(res.datOrg1), 4))), [base, '_meanimg.tif']);
fts = res.fts1;
nEvt = 0;
if isstruct(fts) && isfield(fts, 'basic') && isfield(fts.basic, 'area')
    nEvt = numel(fts.basic.area);
end
ctr = zeros(nEvt, 2); npx = zeros(nEvt, 1); bb = zeros(nEvt, 4);
t0 = zeros(nEvt, 1); t1 = zeros(nEvt, 1); dffMax = zeros(nEvt, 1); area = zeros(nEvt, 1);
evtMap = zeros(H, W, 'uint16');
for k = 1:nEvt
    c = fts.basic.center{k}; ctr(k, :) = c(1:2);
    if isfield(fts.loc, 'xSpa')
        pix = fts.loc.xSpa{k};
    else
        [ih0, iw0, ~, ~] = ind2sub(size(res.datOrg1), res.evt1{k});
        pix = unique(sub2ind([H, W], ih0, iw0));
    end
    npx(k) = numel(pix);
    [ih, iw] = ind2sub([H, W], pix);
    bb(k, :) = [min(ih) max(ih) min(iw) max(iw)];
    evtMap(pix) = evtMap(pix) + 1;
    t0(k) = fts.loc.t0(k); t1(k) = fts.loc.t1(k);
    dffMax(k) = fts.curve.dffMax(k); area(k) = fts.basic.area(k);
end
evT = table((1:nEvt)', t0, t1, ctr(:, 1), ctr(:, 2), npx, ...
    bb(:, 1), bb(:, 2), bb(:, 3), bb(:, 4), dffMax, area, ...
    'VariableNames', {'event', 't0_frame', 't1_frame', 'y_px', 'x_px', 'n_pixels', ...
    'y0', 'y1', 'x0', 'x1', 'dff_max', 'area_um2'});
out = [base, '_events.csv'];
writetable(evT, out);
imwrite(evtMap, [base, '_eventmap.tif']);
fprintf('%d events -> %s\n', nEvt, out);
end
