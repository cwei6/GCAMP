function out = legacy_register_run(sbxPath, outDir, varargin)
%LEGACY_REGISTER_RUN  The lab's legacy MATLAB registration on ONE Scanbox run -> the pipeline's binned tiff.
%
%   legacy_register_run(sbxPath, outDir, 'legacyDir', 'D:\...\lagacy_matlab\Matlab', 'binFrames', 31, ...
%                       'tag', 'bin31f_1p98s', 'behaviorCsv', '...\<run>_bin31f_1p98s_behavior.csv')
%
% reg_compare.py writes the call for you (run_reg_compare_legacy.m); the output is the
% D_legacy variant of the registration comparison: <outDir>\<sbx stem>_<tag>.tif plus the
% sidecars the pipeline expects (.json, _meanproj.tif, _behavior.csv) and
% <stem>_legacy_shifts.csv / _legacy.png (what the two stages found per raw frame).
%
% WHAT IT REPRODUCES — the chain of ProcessScanboxData_concatFirst.m for dataSet 'Afferents'
% with regParam.method = 'rigid' (DefaultProcessingParams: edges [80 80 40 40], binXY 2, no
% filters, chunkSize 100), calling the same legacy functions with the same settings:
%
%   stage 1  CorrectExpt            every raw frame against the mean of its own 15-frame chunk:
%            (DetermineXYShiftsFBS,  edges cropped, imresize 1/2, central 95 %, Gaussian blur 1 px,
%             dftregistrationAlex)   dftregistrationAlex usfac 100 -> sub-pixel row/column shift x 2,
%                                    global median removed, applied with imtranslate (bilinear) to the
%                                    full-resolution frame -> the frames the '.sbxdft' file would hold
%   stage 2  AlignPlanes /           reference = mean of the still epoch of those frames (the longest
%            RegisterSBX 'rigid'     zero-speed epoch >= 30 s with 5 s trimmed off each end, as
%            (binXY, imregtform,     DetermineReference / BinStillEpochs 'speed'); every frame: edges
%             imwarp)                cropped, 2x binned (binXY), imregtform 'rigid' against the binned
%                                    reference (RegularStepGradientDescent: MaximumIterations 1000,
%                                    MinimumStepLength 1e-6, MaximumStepLength 0.001, RelaxationFactor
%                                    0.75; MeanSquares), translation x 2, applied with imwarp (bilinear)
%                                    to the full frame -> the frames the '.sbxreg' file would hold
%   stage 3  GenerateExptProjections mean of every bin of binFrames registered frames, uint16
%            (WriteSbxPlaneTif)      (the legacy script bins to 1 Hz; 31 frames here = the pipeline's bins)
%
% WHAT IS DIFFERENT, on purpose (so the result lines up with the other variants of reg_compare.py):
%   * one run and one channel instead of the concatenated day: the reference comes from this run
%     (its binned _behavior.csv: still_fraction >= 0.9 per 2-s bin stands in for speed == 0);
%   * no .sbxdft / .sbxreg files are written — the stages run chunk by chunk in memory (same maths,
%     same uint16 rounding between the stages as the files would impose);
%   * the first bin of the run is kept (GenerateExptProjections skips it), the bins are the
%     pipeline's bins, and the output is cropped by the shifts (99.5th percentile, every side)
%     instead of the fixed projection edges [80 80 40 40] ('outputCrop' changes that);
%   * every binned frame's brightness level is normalised (median of the frame -> the run's median),
%     which the pipeline's unpack does for the other variants (afferent2p.globalfix step 2 only — no
%     frame flagging here).
%
% PARAMETERS (name/value; defaults = the legacy configuration)
%   legacyDir     folder holding MovieProcessing\, Misc\, Locomotion\ (added to the path); '' = already on it
%   binFrames 31  raw frames per output frame      tag ''         output name tag (default bin<N>f_<sec>s)
%   secondsPerFrame []  (default binFrames / frame rate)           chan 'green' | 'red'
%   edges [80 80 40 40]  [L R T B] px cropped for the estimates    binXY 2        spatial binning for the estimates
%   dft true      stage 1 on/off                   dftChunk 15, dftBlur 1, dftKeep 0.95
%   method 'rigid' | 'translation' | 'none'  stage 2               rigidChunk 100
%   refScan []    raw frames of the reference (1-based), else from behaviorCsv (minStillDur 30, trim 5 s),
%                 else the stillest 40 s by velocity, else the middle 40 s
%   outputCrop 'shift' | [L R T B] | 'none'   cropPct 99.5         normalise true
%   maxScans []   process only the first K raw frames (a quick test; sidecars are cut to match)
%   overwrite false   verbose true
%
% Needs the Image Processing Toolbox (imregtform, imwarp, imtranslate, imresize, imgaussfilt). A
% parallel pool speeds up stage 2 (parfor over the frames of a chunk); without one, expect roughly
% 0.2-0.5 s per raw frame for imregtform (an hour or two for a 16-min run).

t_all = tic;
IP = inputParser;
addRequired(IP, 'sbxPath', @ischar);
addRequired(IP, 'outDir', @ischar);
addParameter(IP, 'legacyDir', '', @ischar);
addParameter(IP, 'binFrames', 31, @isnumeric);
addParameter(IP, 'tag', '', @ischar);
addParameter(IP, 'secondsPerFrame', [], @isnumeric);
addParameter(IP, 'chan', 'green', @ischar);
addParameter(IP, 'edges', [80 80 40 40], @isnumeric);
addParameter(IP, 'binXY', 2, @isnumeric);
addParameter(IP, 'dft', true, @islogical);
addParameter(IP, 'dftChunk', 15, @isnumeric);
addParameter(IP, 'dftBlur', 1, @isnumeric);
addParameter(IP, 'dftKeep', 0.95, @isnumeric);
addParameter(IP, 'method', 'rigid', @ischar);
addParameter(IP, 'rigidChunk', 100, @isnumeric);
addParameter(IP, 'refScan', [], @isnumeric);
addParameter(IP, 'minStillDur', 30, @isnumeric);
addParameter(IP, 'trim', 5, @isnumeric);
addParameter(IP, 'behaviorCsv', '', @ischar);
addParameter(IP, 'outputCrop', 'shift');
addParameter(IP, 'cropPct', 99.5, @isnumeric);
addParameter(IP, 'normalise', true, @islogical);
addParameter(IP, 'maxScans', [], @isnumeric);
addParameter(IP, 'overwrite', false, @islogical);
addParameter(IP, 'verbose', true, @islogical);
parse(IP, sbxPath, outDir, varargin{:});
P = IP.Results;
N = round(P.binFrames);
e = round(P.edges(:))';
verbose = P.verbose;

% ---------------------------------------------------------------- the legacy functions
if ~isempty(P.legacyDir)
    for sub = {'MovieProcessing', 'Misc', 'Locomotion'}
        d = fullfile(P.legacyDir, sub{1});
        if isfolder(d), addpath(d); end
    end
end
need = {'readSBX', 'MakeChunkLims', 'DetermineXYShiftsFBS', 'dftregistrationAlex', 'binXY'};
missing = need(cellfun(@(f) exist(f, 'file') == 0, need));
if ~isempty(missing)
    error('legacy_register_run:missing', 'legacy functions not on the MATLAB path: %s (give ''legacyDir'', the folder with MovieProcessing\\)', strjoin(missing, ', '));
end
for f = {'imregtform', 'imwarp', 'imtranslate', 'imresize', 'imgaussfilt'}
    if exist(f{1}, 'file') == 0 && exist(f{1}, 'builtin') == 0
        error('legacy_register_run:toolbox', '%s not found: the Image Processing Toolbox is needed', f{1});
    end
end

% ---------------------------------------------------------------- the recording
if ~exist(sbxPath, 'file'), error('legacy_register_run:nofile', '%s does not exist', sbxPath); end
info = sbx_info(sbxPath);                                                             % LoadSBXinfo's logic in double precision (see there)
if info.Nplane ~= 1, error('legacy_register_run:planes', 'single-plane recordings only (this one has %d planes)', info.Nplane); end
T = info.totScan;
if ~isempty(P.maxScans), T = min(T, round(P.maxScans)); end
H = info.sz(1); W = info.sz(2);
fs = frame_rate(info);
spf = P.secondsPerFrame;
if isempty(spf), spf = N / fs; end
tag = P.tag;
if isempty(tag), tag = sprintf('bin%df_%ss', N, strrep(sprintf('%.2f', N / fs), '.', 'p')); end
[~, stem0] = fileparts(sbxPath);
stem = sprintf('%s_%s', stem0, tag);
if ~isfolder(outDir), mkdir(outDir); end
outTif = fullfile(outDir, [stem '.tif']);
if exist(outTif, 'file') && ~P.overwrite
    fprintf('legacy_register_run: %s exists - kept (overwrite=true to redo)\n', outTif);
    out = struct('tif', outTif, 'skipped', true);
    return
end
pmt = which_pmt(P.chan, info);
nOut = floor(T / N);
if nOut < 1, error('legacy_register_run:short', 'only %d frames: fewer than one bin of %d', T, N); end
Hc = H - e(3) - e(4); Wc = W - e(1) - e(2);
if Hc < 32 || Wc < 32, error('legacy_register_run:edges', 'edges %s leave no frame (%dx%d)', mat2str(e), H, W); end
log_(verbose, '%s: %d raw frames @ %.3f Hz, %dx%d, channel %s (pmt %d) -> %d bins of %d (%.3f s) in %s', ...
     [stem0 '.sbx'], T, fs, H, W, P.chan, pmt, nOut, N, spf, outDir);

% ---------------------------------------------------------------- behaviour (reference epoch, sidecar)
beh = struct('ok', false);
behCsv = P.behaviorCsv;
if isempty(behCsv)
    cand = dir(fullfile(fileparts(sbxPath), sprintf('%s_%s_behavior.csv', stem0, tag)));
    if isempty(cand), cand = dir(fullfile(fileparts(sbxPath), sprintf('%s_bin*_behavior.csv', stem0))); end
    if ~isempty(cand), behCsv = fullfile(cand(1).folder, cand(1).name); end
end
if ~isempty(behCsv) && exist(behCsv, 'file')
    beh = read_behavior(behCsv);
    if ~beh.ok, log_(verbose, 'behaviour csv %s could not be read (no still_fraction column?)', behCsv); end
end

% ---------------------------------------------------------------- stage 1: DFT jitter within 15-frame chunks (CorrectExpt)
t1 = tic;
RS = zeros(1, T); CS = zeros(1, T);
if P.dft
    [chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, T, T, 'size', P.dftChunk, 'allowPartial', true);
    log_(verbose, 'stage 1 (CorrectExpt): DFT shifts against the chunk mean, %d chunks of %d frames ...', Nchunk, P.dftChunk);
    tnext = tic;
    for c = 1:Nchunk
        raw = readSBX(sbxPath, info, chunkLims(c, 1), chunkLength(c), pmt, []);      % H x W x n, uint16 (65535 - stored)
        raw = raw(e(3) + 1:end - e(4), e(1) + 1:end - e(2), :);
        n = size(raw, 3);
        small = imresize(raw(:, :, 1), 1 / P.binXY);                                 % as CorrectExpt: imresize(raw_chunk, 1/binXY)
        chunk = zeros(size(small, 1), size(small, 2), 1, n, class(raw));
        chunk(:, :, 1, 1) = small;
        for k = 2:n
            chunk(:, :, 1, k) = imresize(raw(:, :, k), 1 / P.binXY);
        end
        ref1 = mean(chunk, 4);
        [RS1, CS1] = DetermineXYShiftsFBS(chunk, P.dftBlur, P.dftKeep, ref1);       % 1 x n, sub-pixel, binned px
        RS(chunkLims(c, 1):chunkLims(c, 2)) = RS1 * P.binXY;
        CS(chunkLims(c, 1):chunkLims(c, 2)) = CS1 * P.binXY;
        if toc(tnext) > 30
            log_(verbose, '   chunk %d / %d (%.0f s)', c, Nchunk, toc(t1)); tnext = tic;
        end
    end
    RS = RS - median(RS); CS = CS - median(CS);                                     % CorrectExpt second stage (2D: no inter-chunk term)
    log_(verbose, 'stage 1 done in %.0f s: |shift| median %.2f px, max %.2f px', toc(t1), median(hypot(RS, CS)), max(hypot(RS, CS)));
else
    log_(verbose, 'stage 1 (DFT) skipped');
end

% ---------------------------------------------------------------- the reference epoch (DetermineReference)
[refScan, refWhy] = reference_scans(P.refScan, beh, N, T, spf, P.minStillDur, P.trim);
log_(verbose, 'reference: %d frames (%d-%d) = %s', numel(refScan), refScan(1), refScan(end), refWhy);

% ---------------------------------------------------------------- stage 2 + 3: rigid registration to the reference, binning
doRigid = ~strcmpi(P.method, 'none');
refImage = mean_of_frames(sbxPath, info, pmt, refScan, RS, CS);                   % uint16, as AlignPlanes' refVol
if doRigid
    rawRef = uint16(refImage(e(3) + 1:end - e(4), e(1) + 1:end - e(2)));           % RegisterSBX: crop the reference
    preRef = binXY(rawRef, P.binXY);                                                % spatial downsampling
    [optimizer, metric] = legacy_optimizer();
    method = lower(P.method);
end
ty = zeros(1, T); tx = zeros(1, T); ang = zeros(1, T);                              % stage 2 translation (full-res px) and rotation (deg)
nFailed = 0;
binned = zeros(H, W, nOut, 'uint16');
accBin = zeros(H, W);
Tbin = nOut * N;
[chunkLims, Nchunk, chunkLength] = MakeChunkLims(1, Tbin, Tbin, 'size', P.rigidChunk, 'allowPartial', true);
log_(verbose, 'stage 2 (RegisterSBX ''%s'') + binning: %d chunks of %d frames ...', P.method, Nchunk, P.rigidChunk);
t2 = tic; tnext = tic;
imRef = imref2d([H, W]);
for c = 1:Nchunk
    s0 = chunkLims(c, 1); n = chunkLength(c);
    raw = readSBX(sbxPath, info, s0, n, pmt, []);
    if n == 1, raw = reshape(raw, H, W, 1); end
    dftFr = raw;
    if P.dft
        for k = 1:n
            dftFr(:, :, k) = imtranslate(raw(:, :, k), [CS(s0 + k - 1), RS(s0 + k - 1)]);   % CorrectExpt: imtranslate(data, [CS, RS]) -> uint16
        end
    end
    reg = dftFr;
    if doRigid
        pre = binXY(dftFr(e(3) + 1:end - e(4), e(1) + 1:end - e(2), :), P.binXY);
        tf = cell(1, n); failed = false(1, n);
        parfor k = 1:n
            try
                tf{k} = imregtform(pre(:, :, k), preRef, method, optimizer, metric);      % RegisterSBX, rigid branch
            catch
                tf{k} = affine2d(eye(3)); failed(k) = true;                                  % (a blank / unusable frame: left where it is)
            end
        end
        nFailed = nFailed + sum(failed);
        for k = 1:n
            Tm = tform_matrix(tf{k});                                                        % 3x3, [x y 1] * T convention (affine2d.T)
            Tm(3, 1:2) = Tm(3, 1:2) * P.binXY;                                               % "Correct reg_transforms for spatial downsampling"
            tform = affine2d(Tm);
            reg(:, :, k) = imwarp(dftFr(:, :, k), tform, 'OutputView', imRef);               % AlignPlanes: imwarp(data, regTform, 'OutputView', imRef)
            s = s0 + k - 1;
            tx(s) = Tm(3, 1); ty(s) = Tm(3, 2); ang(s) = atan2d(Tm(1, 2), Tm(1, 1));
        end
    end
    for k = 1:n
        s = s0 + k - 1;
        accBin = accBin + double(reg(:, :, k));
        if mod(s, N) == 0
            binned(:, :, s / N) = uint16(accBin / N);                                        % WriteSbxPlaneTif: uint16(mean(chunk, 3))
            accBin(:) = 0;
        end
    end
    if toc(tnext) > 30
        log_(verbose, '   chunk %d / %d, frame %d / %d (%.0f s)', c, Nchunk, s0 + n - 1, Tbin, toc(t2)); tnext = tic;
    end
end
log_(verbose, 'stage 2 done in %.0f s: rigid |translation| median %.2f px, max %.2f px, |rotation| max %.3f deg%s', ...
     toc(t2), median(hypot(tx, ty)), max(hypot(tx, ty)), max(abs(ang)), failed_note(nFailed));

% ---------------------------------------------------------------- output crop
totY = RS(1:Tbin) + ty(1:Tbin); totX = CS(1:Tbin) + tx(1:Tbin);
if ischar(P.outputCrop) && strcmpi(P.outputCrop, 'shift')
    crop = ceil(pctile_(abs([totY, totX]), P.cropPct));
    cropLRTB = [crop crop crop crop];
    cropWhy = sprintf('%.1fth percentile of the shifts, every side', P.cropPct);
elseif ischar(P.outputCrop)
    cropLRTB = [0 0 0 0]; cropWhy = 'none';
else
    cropLRTB = round(P.outputCrop(:))'; cropWhy = 'fixed [L R T B] (the legacy projection edges)';
end
cropLRTB = max(cropLRTB, 0);
binned = binned(cropLRTB(3) + 1:end - cropLRTB(4), cropLRTB(1) + 1:end - cropLRTB(2), :);
Ho = size(binned, 1); Wo = size(binned, 2);

% ---------------------------------------------------------------- level normalisation (afferent2p.globalfix, step 2)
lvl = zeros(1, nOut);
for b = 1:nOut
    fr = binned(:, :, b); lvl(b) = median(double(fr(:)));
end
G = median(lvl);
if P.normalise
    for b = 1:nOut
        binned(:, :, b) = uint16(double(binned(:, :, b)) * (G / max(lvl(b), 1e-9)));
    end
end

% ---------------------------------------------------------------- write the movie and the sidecars
write_tif_stack(outTif, binned, spf);
proj = uint16(mean(double(binned), 3));
write_tif_stack(fullfile(outDir, [stem '_meanproj.tif']), proj, []);
shiftsCsv = fullfile(outDir, [stem '_legacy_shifts.csv']);
fid = fopen(shiftsCsv, 'w');
fprintf(fid, 'frame,time_s,dft_dy,dft_dx,rigid_dy,rigid_dx,rigid_deg,total_dy,total_dx,running\n');
isRun = -ones(1, T);
if beh.ok
    sf = repelem(beh.still_fraction(:)', N); sf = sf(1:min(end, T));
    isRun(1:numel(sf)) = double(sf < 0.5);
end
for s = 1:T
    fprintf(fid, '%d,%.4f,%.3f,%.3f,%.3f,%.3f,%.4f,%.3f,%.3f,%d\n', s - 1, (s - 1) / fs, RS(s), CS(s), ty(s), tx(s), ang(s), RS(s) + ty(s), CS(s) + tx(s), isRun(s));
end
fclose(fid);
if beh.ok
    write_behavior(beh, nOut, fullfile(outDir, [stem '_behavior.csv']));
end
amp = hypot(totY, totX);
wy = reshape(totY, N, nOut); wx = reshape(totX, N, nOut);
stats = struct('shift_median_px', median(amp), 'shift_p95_px', pctile_(amp, 95), 'shift_max_px', max(amp), ...
               'dft_median_px', median(hypot(RS, CS)), 'dft_max_px', max(hypot(RS, CS)), ...
               'rigid_median_px', median(hypot(tx, ty)), 'rigid_max_px', max(hypot(tx, ty)), 'rotation_max_deg', max(abs(ang)), ...
               'jump_p95_px', pctile_(hypot(diff(totY), diff(totX)), 95), ...
               'within_bin_range_median_px', median(max(max(wy, [], 1) - min(wy, [], 1), max(wx, [], 1) - min(wx, [], 1))));
J = struct();
J.source_sbx = sbxPath; J.channel = pmt; J.channels_in_file = info.nchan; J.frames_in = T; J.fs_in_hz = fs;
J.bin_frames = N; J.seconds_per_frame = spf; J.fs_out_hz = fs / N; J.frames_out = nOut;
J.crop_LRTB_px = cropLRTB; J.height = Ho; J.width = Wo; J.dtype = 'uint16';
J.values = 'mean of (65535 - raw) over each bin, every raw frame registered first (legacy MATLAB chain)';
J.registration = struct('method', sprintf('legacy MATLAB chain: CorrectExpt DFT (chunks of %d, sub-pixel, imtranslate) then RegisterSBX ''%s'' (imregtform on edges %s, binXY %d; imwarp), legacy_register_run.m', ...
                                          P.dftChunk, P.method, mat2str(e), P.binXY), ...
                        'reference', sprintf('%d raw frames %d-%d: %s', numel(refScan), refScan(1), refScan(end), refWhy), ...
                        'stages', struct('dft', P.dft, 'rigid', doRigid, 'rigid_failed_frames', nFailed), 'stats', stats, 'crop_in_raw_px', cropLRTB(1), ...
                        'crop_rule', cropWhy, 'filled_border', 'zeros (imtranslate / imwarp fill)', 'shifts_csv', [stem '_legacy_shifts.csv'], ...
                        'seconds', round(toc(t_all)));
J.global_correction = struct('method', 'level normalisation only (MATLAB port of afferent2p.globalfix step 2: frame median -> run median); no frame flagging', ...
                             'level', G, 'applied', P.normalise, 'written', datestr(now, 'yyyy-mm-dd HH:MM'));
J.tool = 'legacy_register_run.m';
write_json(fullfile(outDir, [stem '.json']), J);
try
    shift_figure(fullfile(outDir, [stem '_legacy.png']), stem, fs, RS, CS, ty, tx, beh, N, stats);
catch err
    log_(verbose, '   figure failed (%s)', err.message);
end
log_(verbose, 'wrote %s (%d bins, %dx%d after cropping %s px; |shift| median %.1f px, max %.1f; %.0f s)', ...
     outTif, nOut, Ho, Wo, mat2str(cropLRTB), stats.shift_median_px, stats.shift_max_px, toc(t_all));
out = struct('tif', outTif, 'json', fullfile(outDir, [stem '.json']), 'shifts_csv', shiftsCsv, 'frames_out', nOut, ...
             'crop_LRTB_px', cropLRTB, 'stats', stats, 'reference_frames', [refScan(1) refScan(end)], 'skipped', false);
end


% ====================================================================== helpers

function log_(verbose, fmt, varargin)
if verbose, fprintf(['[%s] ' fmt '\n'], datestr(now, 'HH:MM:SS'), varargin{:}); end
end


function s = failed_note(n)
if n > 0, s = sprintf(' (imregtform failed on %d frame(s): left unshifted)', n); else, s = ''; end
end


function info = sbx_info(sbxPath)
% The Scanbox info as the legacy readSBX needs it (nchan, nsamples, Nplane, totScan, optotune_used, ...),
% built like LoadSBXinfo / MakeInfoStruct but with every numeric field cast to double first: an info
% .mat that was re-saved from Python stores sz / recordsPerBuffer / resfreq as uint16, and the legacy
% loader's uint16 arithmetic would saturate (nsamples = 65535) or error on mixed integer classes.
[d, n] = fileparts(sbxPath);
S = load(fullfile(d, [n '.mat']));
if ~isfield(S, 'info'), error('legacy_register_run:info', '%s has no info structure', fullfile(d, [n '.mat'])); end
info = to_double(S.info);
if ~isfield(info, 'scanmode'), info.scanmode = 1; end
if info.scanmode == 0, info.recordsPerBuffer = info.recordsPerBuffer * 2; end
switch info.channels
    case 1, info.nchan = 2;                        % both PMTs
    case 2, info.nchan = 1;                        % PMT 0
    case 3, info.nchan = 1;                        % PMT 1
    case -1, info.nchan = info.chan.nchan;         % scanbox 3: the chan structure says
    otherwise, error('legacy_register_run:info', 'unexpected info.channels = %g', info.channels);
end
if info.nchan == 1, factor = 2; else, factor = 1; end
info.nsamples = info.sz(2) * info.recordsPerBuffer * 2 * info.nchan;   % bytes per record
info.Nrow = info.sz(1); info.Ncol = info.sz(2);
info.path = sbxPath;
if ~isfield(info, 'max_idx')
    f = dir(sbxPath);
    if isfield(info, 'scanbox_version') && info.scanbox_version >= 2
        info.max_idx = f.bytes / info.recordsPerBuffer / info.sz(2) * factor / 4 - 1;
    else
        info.max_idx = f.bytes / info.bytesPerBuffer * factor - 1;
    end
end
info.nframes = floor(info.max_idx + 1);
info.optotune_used = false; info.otlevels = 1;
if isfield(info, 'volscan') && info.volscan > 0, info.optotune_used = true; end
if ~isfield(info, 'volscan') && isfield(info, 'otwave') && ~isempty(info.otwave), info.optotune_used = true; end
if info.optotune_used
    if isfield(info, 'otwave') && ~isempty(info.otwave), info.otlevels = numel(info.otwave);
    elseif isfield(info, 'otparam') && numel(info.otparam) >= 3 && isfinite(info.otparam(3)) && info.otparam(3) > 0, info.otlevels = info.otparam(3);
    else, info.otlevels = 30;
    end
end
if ~isfield(info, 'framerate')
    info.framerate = 15.49; if info.scanmode == 0, info.framerate = 30.98; end
end
info.height = info.sz(2); info.width = info.recordsPerBuffer;
info.Nplane = info.otlevels;
info.Nscan = floor(info.nframes / info.Nplane);
info.totScan = info.Nscan;
end


function x = to_double(x)
% numeric fields of a struct (recursively) as double; everything else untouched
if isstruct(x)
    if numel(x) ~= 1, return; end
    f = fieldnames(x);
    for k = 1:numel(f)
        try
            x.(f{k}) = to_double(x.(f{k}));
        catch
        end
    end
elseif isnumeric(x) || islogical(x)
    x = double(x);
end
end


function fs = frame_rate(info)
% The frame rate as the pipeline computes it (sbxreader): resonant frequency / lines per frame,
% doubled for bidirectional scanning, per plane. LoadSBXinfo's framerate (15.49 default) as a fallback.
fs = [];
try
    if isfield(info, 'resfreq') && isfield(info, 'config') && isfield(info.config, 'lines') && info.config.lines > 0
        factor = 1; if isfield(info, 'scanmode') && info.scanmode == 0, factor = 2; end
        fs = factor * info.resfreq / info.config.lines / info.Nplane;
    end
catch
    fs = [];
end
if isempty(fs) || ~isfinite(fs) || fs <= 0
    if isfield(info, 'framerate'), fs = info.framerate; else, fs = 15.49; end
end
end


function pmt = which_pmt(chan, info)
% DeterminePMT when the info has the PMT gains, else green = 1, red = 2 (or the only channel).
pmt = [];
if exist('DeterminePMT', 'file') && isfield(info, 'config') && isfield(info.config, 'pmt0_gain') && isfield(info.config, 'pmt1_gain')
    try
        pmt = DeterminePMT(chan, info);
    catch
        pmt = [];
    end
end
if isempty(pmt)
    if strcmpi(chan, 'red'), pmt = 2; else, pmt = 1; end
end
if info.nchan == 1, pmt = 1; end
pmt = pmt(1);
end


function v = pctile_(x, p)
% percentile (linear interpolation between order statistics, as MATLAB's prctile) without the toolbox
x = sort(double(x(:)));
x = x(isfinite(x));
n = numel(x);
if n == 0, v = NaN; return; end
if n == 1, v = x; return; end
r = p / 100 * n + 0.5;                      % MATLAB prctile: the k-th sorted value is the 100*(k-0.5)/n percentile
r = min(max(r, 1), n);
lo = floor(r); hi = ceil(r);
v = x(lo) + (r - lo) * (x(hi) - x(lo));
end


function beh = read_behavior(csvPath)
% The pipeline's binned behaviour csv (one row per output frame): header names + numeric columns
% (NaN allowed), as fields of a struct.
beh = struct('ok', false, 'path', csvPath, 'header', {{}}, 'data', []);
try
    txt = fileread(csvPath);
    lines = regexp(txt, '\r?\n', 'split');
    lines = lines(~cellfun(@isempty, strtrim(lines)));
    names = strtrim(strsplit(lines{1}, ','));
    data = NaN(numel(lines) - 1, numel(names));
    for r = 2:numel(lines)
        vals = str2double(strsplit(lines{r}, ','));
        m = min(numel(vals), numel(names));
        data(r - 1, 1:m) = vals(1:m);
    end
    beh.header = names; beh.data = data; beh.ok = true;
    for k = 1:numel(names)
        nm = regexprep(names{k}, '[^A-Za-z0-9_]', '_');
        if ~isempty(nm) && ~any(strcmp(nm, {'ok', 'path', 'header', 'data'}))
            beh.(nm) = data(:, k)';
        end
    end
    if ~isfield(beh, 'still_fraction'), beh.ok = false; end
catch
    beh.ok = false;
end
end


function write_behavior(beh, nOut, outPath)
fid = fopen(outPath, 'w');
fprintf(fid, '%s\n', strjoin(beh.header, ','));
d = beh.data(1:min(nOut, size(beh.data, 1)), :);
fmt = [strjoin(repmat({'%.6g'}, 1, size(d, 2)), ',') '\n'];
for r = 1:size(d, 1), fprintf(fid, fmt, d(r, :)); end
fclose(fid);
end


function [refScan, why] = reference_scans(given, beh, N, T, spf, minStillDur, trim)
% The raw frames of the reference: as given, else the longest still epoch of the behaviour csv
% (>= minStillDur s, shorter if none, trimmed by `trim` s at both ends — DetermineReference), else
% the 40 s with the least running, else the middle 40 s.
if ~isempty(given)
    refScan = round(given(:)'); refScan = refScan(refScan >= 1 & refScan <= T);
    why = 'given'; return
end
nb = floor(T / N);
if beh.ok && numel(beh.still_fraction) >= 1 && nb >= 1
    sf = beh.still_fraction(1:min(nb, numel(beh.still_fraction)));
    still = sf >= 0.9;
    epochs = runs_of(still);                                                         % [first last] bins
    minDur = minStillDur; tr = trim; pick = [];
    while isempty(pick) && minDur > 0 && ~isempty(epochs)
        dur = (epochs(:, 2) - epochs(:, 1) + 1) * spf;
        ok = find(dur >= minDur);
        if ~isempty(ok)
            [~, j] = max(dur(ok)); pick = epochs(ok(j), :);
        else
            minDur = minDur - 1;
            if tr >= minDur / 4, tr = minDur / 10; end
        end
    end
    if ~isempty(pick)
        trimBins = round(tr / spf);
        if pick(2) - pick(1) + 1 - 2 * trimBins < 1, trimBins = floor((pick(2) - pick(1)) / 2); end
        b0 = pick(1) + trimBins; b1 = pick(2) - trimBins;
        refScan = (b0 - 1) * N + 1:b1 * N;
        why = sprintf('longest still epoch of the behaviour csv (%.0f s, %.0f s trimmed each side)', (pick(2) - pick(1) + 1) * spf, trimBins * spf);
        return
    end
    if isfield(beh, 'velocity_mean')
        v = abs(beh.velocity_mean(1:numel(sf))); v(~isfinite(v)) = 0;
        win = max(1, round(40 / spf));
        if numel(v) >= win
            cs = cumsum([0 v]); s = cs(win + 1:end) - cs(1:end - win);
            [~, b0] = min(s);
            refScan = (b0 - 1) * N + 1:min((b0 + win - 1) * N, T);
            why = 'the 40 s with the least running (no still epoch in the behaviour csv)';
            return
        end
    end
end
mid = floor(T / 2); half = 20 * N;                                                   % 20 bins each side
refScan = max(1, mid - half):min(T, mid + half);
why = 'middle of the run (no behaviour csv)';
end


function ep = runs_of(mask)
mask = logical(mask(:)');
d = diff([false mask false]);
ep = [find(d == 1)', find(d == -1)' - 1];
end


function img = mean_of_frames(sbxPath, info, pmt, scans, RS, CS)
% uint16 mean of the stage-1-corrected frames `scans` (AlignPlanes: refVol = uint16(WriteSbxProjection(... .sbxdft ...)))
H = info.sz(1); W = info.sz(2);
acc = zeros(H, W); n = 0;
scans = scans(:)';
k = 1;
while k <= numel(scans)
    s0 = scans(k); j = k;
    while j < numel(scans) && scans(j + 1) == scans(j) + 1 && j - k < 199, j = j + 1; end
    cnt = j - k + 1;
    raw = readSBX(sbxPath, info, s0, cnt, pmt, []);
    if cnt == 1, raw = reshape(raw, H, W, 1); end
    for q = 1:cnt
        s = s0 + q - 1;
        fr = raw(:, :, q);
        if RS(s) ~= 0 || CS(s) ~= 0, fr = imtranslate(fr, [CS(s), RS(s)]); end
        acc = acc + double(fr); n = n + 1;
    end
    k = j + 1;
end
img = uint16(acc / max(n, 1));
end


function Tm = tform_matrix(tf)
% The 3x3 matrix of a geometric transformation in affine2d's convention ([x y 1] * T):
% affine2d (R2022a and earlier imregtform output) has T; rigidtform2d / affinetform2d
% (R2022b+, premultiply convention) have A = T.'.
if isprop(tf, 'T')
    Tm = tf.T;
elseif isprop(tf, 'A')
    Tm = tf.A.';
elseif isstruct(tf) && isfield(tf, 'T')
    Tm = tf.T;
else
    error('legacy_register_run:tform', 'unknown transformation object of class %s', class(tf));
end
Tm = double(Tm);
end


function [optimizer, metric] = legacy_optimizer()
% RegisterSBX, non-TurboReg branch (used for method 'rigid' / 'translation')
optimizer = registration.optimizer.RegularStepGradientDescent;
optimizer.MaximumIterations = 1000;
optimizer.MinimumStepLength = 1.0000e-06;
optimizer.MaximumStepLength = 0.001;
optimizer.RelaxationFactor = 0.75;
metric = registration.metric.MeanSquares;
end


function write_tif_stack(path, stack, spf)
% uint16 multi-page TIFF written directly: classic little-endian TIFF, uncompressed, one strip per
% page, all IFDs first and the page data contiguous after them (what AQuA2's reader, tifffile —
% including its memmap — and ImageJ expect); no Tiff class / imwrite quirks.
H = size(stack, 1); W = size(stack, 2); n = size(stack, 3);
if ~isa(stack, 'uint16'), stack = uint16(stack); end
if isempty(spf)
    desc = sprintf('legacy_register_run.m frames=%d', n);
else
    desc = sprintf('legacy_register_run.m frames=%d seconds_per_frame=%.6f', n, spf);
end
descBytes = [uint8(desc) 0];
nd = numel(descBytes);
bytesPerPage = H * W * 2;
ifdBytes = 2 + 12 * 12 + 4;                     % 12 entries per page, 13 on the first (ImageDescription)
ifdOff = 8 + cumsum([0, repmat(ifdBytes, 1, n)]) + [0, repmat(12, 1, n)];   % the first IFD is 12 bytes longer
descOff = ifdOff(n + 1);
dataOff0 = descOff + nd;
total = dataOff0 + n * bytesPerPage;
if total >= 2^32 - 1
    error('legacy_register_run:tiff', '%s would be %.1f GB: too large for a classic TIFF (bin more frames or crop)', path, total / 1e9);
end
if exist(path, 'file'), delete(path); end
fid = fopen(path, 'w', 'ieee-le');
if fid < 0, error('legacy_register_run:tiff', 'cannot write %s', path); end
fwrite(fid, uint8('II'), 'uint8'); fwrite(fid, 42, 'uint16'); fwrite(fid, 8, 'uint32');
for k = 1:n
    dataOff = dataOff0 + (k - 1) * bytesPerPage;
    nextIfd = ifdOff(k + 1); if k == n, nextIfd = 0; end
    fwrite(fid, 12 + (k == 1), 'uint16');
    tiff_entry(fid, 254, 4, 1, 0);                 % NewSubfileType
    tiff_entry(fid, 256, 4, 1, W);                 % ImageWidth
    tiff_entry(fid, 257, 4, 1, H);                 % ImageLength
    tiff_entry(fid, 258, 3, 1, 16);                % BitsPerSample
    tiff_entry(fid, 259, 3, 1, 1);                 % Compression: none
    tiff_entry(fid, 262, 3, 1, 1);                 % Photometric: BlackIsZero
    if k == 1, tiff_entry(fid, 270, 2, nd, descOff); end   % ImageDescription
    tiff_entry(fid, 273, 4, 1, dataOff);           % StripOffsets
    tiff_entry(fid, 277, 3, 1, 1);                 % SamplesPerPixel
    tiff_entry(fid, 278, 4, 1, H);                 % RowsPerStrip
    tiff_entry(fid, 279, 4, 1, bytesPerPage);      % StripByteCounts
    tiff_entry(fid, 284, 3, 1, 1);                 % PlanarConfiguration: chunky
    tiff_entry(fid, 339, 3, 1, 1);                 % SampleFormat: unsigned integer
    fwrite(fid, nextIfd, 'uint32');
end
fwrite(fid, descBytes, 'uint8');
for k = 1:n
    fwrite(fid, stack(:, :, k).', 'uint16');       % rows contiguous
end
fclose(fid);
end


function tiff_entry(fid, tag, typ, count, value)
% one 12-byte IFD entry; typ 2 = ASCII (value = offset), 3 = SHORT, 4 = LONG
fwrite(fid, tag, 'uint16'); fwrite(fid, typ, 'uint16'); fwrite(fid, count, 'uint32');
if typ == 3 && count == 1
    fwrite(fid, value, 'uint16'); fwrite(fid, 0, 'uint16');
else
    fwrite(fid, value, 'uint32');
end
end


function write_json(path, J)
try
    txt = jsonencode(J, 'PrettyPrint', true);
catch
    txt = jsonencode(J);
end
fid = fopen(path, 'w'); fwrite(fid, txt, 'char'); fclose(fid);
end


function shift_figure(pngPath, stem, fs, RS, CS, ty, tx, beh, N, stats)
T = numel(RS); t = (0:T - 1) / fs / 60;
fig = figure('Visible', 'off', 'Position', [100 100 1100 520], 'Color', 'w');
ax1 = subplot(3, 1, 1); plot(t, RS, '-', 'LineWidth', 0.5); hold on; plot(t, CS, '-', 'LineWidth', 0.5);
ylabel('stage 1: DFT shift (px)'); legend({'dy (rows)', 'dx (columns)'}, 'Location', 'northeast', 'Box', 'off');
title(sprintf('%s - legacy registration: |shift| median %.1f px, max %.1f (DFT within 15-frame chunks + rigid to the still reference)', stem, stats.shift_median_px, stats.shift_max_px), 'Interpreter', 'none', 'FontSize', 9);
ax2 = subplot(3, 1, 2); plot(t, ty, '-', 'LineWidth', 0.5); hold on; plot(t, tx, '-', 'LineWidth', 0.5);
ylabel('stage 2: rigid translation (px)');
ax3 = subplot(3, 1, 3);
if beh.ok && isfield(beh, 'velocity_mean')
    tb = ((0:numel(beh.velocity_mean) - 1) + 0.5) * N / fs / 60;
    v = abs(beh.velocity_mean); v(~isfinite(v)) = 0;
    area(tb, v, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
ylabel('running speed'); xlabel('minutes');
linkaxes([ax1 ax2 ax3], 'x'); xlim([0, max(t)]);
print(fig, '-dpng', '-r100', pngPath);
close(fig);
end
