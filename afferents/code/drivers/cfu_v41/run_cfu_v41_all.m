%% run_cfu_v41_all.m — the CFU step of the three RearChr animals with AQuA2's published CFU rule,
% overlap 0.3 / min events 1 (afferent2p v41, decided 2 Nov 2026). Every result is redone (its saved
% cfuOpts differ from the new values). Output in run_cfu_v41_all_log.txt next to this file.
logFile = 'D:\2photon\Chao\GCAMP_affrerents\run_cfu_v41_all_log.txt';
try
    diary(logFile);
catch
    diary(strrep(logFile, '.txt', ['_' datestr(now, 'yyyymmdd_HHMM') '.txt']));
end
setappdata(0, 'cfu_v41_t0', now);          % survives the clearvars of the batch scripts
fprintf('CFU v41 batches started %s\n', datestr(now));
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-1\aqua_batch\run_cfu_batch_v41.m');
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2608-2\aqua_batch\run_cfu_batch_v41.m');
run('D:\2photon\Chao\GCAMP_affrerents\NavRearChrTIGRE-2609-1\aqua_batch\run_cfu_batch_v41.m');
fprintf('CFU v41 batches done for the three animals: %s (%.0f min).\n', datestr(now), (now - getappdata(0, 'cfu_v41_t0')) * 24 * 60);
diary off
