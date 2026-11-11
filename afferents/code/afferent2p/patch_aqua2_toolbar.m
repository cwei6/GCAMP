function patch_aqua2_toolbar(aquaDir)
%PATCH_AQUA2_TOOLBAR  Make aqua_gui start on MATLAB releases where the axes
% toolbar object is created lazily.
%
%   patch_aqua2_toolbar('D:\AQuA2')
%
% AQuA2 hides the per-axes hover toolbar with lines like
%       pMov1.Toolbar.Visible = 'off';
% (src\+ui\+com\addCon_dat.m). On some MATLAB releases/display setups the
% AxesToolbar object does not exist yet at that point, so the assignment fails:
%       Property assignment is not allowed when the object is empty.
% Setting the property to [] is the documented way to remove the toolbar and
% works whether or not the object exists, so every such line becomes
%       pMov1.Toolbar = [];
% A backup <file>.m.bak is written next to each changed file. Re-run this after
% every AQuA2 update (a fresh download brings the original lines back).

if nargin < 1 || isempty(aquaDir)
    aquaDir = fileparts(which('aqua_gui'));
    if isempty(aquaDir)
        error('patch_aqua2_toolbar:notFound', ...
            'Pass the AQuA2 folder, e.g. patch_aqua2_toolbar(''D:\\AQuA2'')');
    end
end
srcDir = fullfile(aquaDir, 'src');
if ~isfolder(srcDir)
    error('patch_aqua2_toolbar:noSrc', 'No src folder under %s — is this the AQuA2 folder?', aquaDir);
end
files = dir(fullfile(srcDir, '**', '*.m'));
pat = '(\w+)\.Toolbar\.Visible\s*=\s*''off''\s*;';
rep = '$1.Toolbar = [];  % patched by patch_aqua2_toolbar: toolbar object can be empty';
nFiles = 0; nLines = 0;
for k = 1:numel(files)
    p = fullfile(files(k).folder, files(k).name);
    txt = fileread(p);
    hits = regexp(txt, pat, 'match');
    if isempty(hits)
        continue
    end
    copyfile(p, [p, '.bak']);
    fid = fopen(p, 'w');
    fwrite(fid, regexprep(txt, pat, rep));
    fclose(fid);
    nFiles = nFiles + 1; nLines = nLines + numel(hits);
    fprintf('patched %d line(s): %s\n', numel(hits), p);
end
if nFiles == 0
    fprintf('nothing to patch under %s (already patched, or a different AQuA2 version)\n', srcDir);
else
    fprintf('done: %d line(s) in %d file(s). Run "clear classes; rehash" or restart MATLAB, then aqua_gui.\n', nLines, nFiles);
end
end
