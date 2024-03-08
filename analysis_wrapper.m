datapath = 'example_data';
outputpath = './';
numframes = 37;
titlepost = '';
%%
% postfix = '';
% fprintf('Making video...\n')
% run('analysis_videography.m');
% %%
% postfix = '';
% fprintf('Analyzing PBL...\n')
% run('analysis_pbl');
%%
postfix = '';
titlepost = '';
fprintf('Analyzing Energy...\n')
run('analysis_energetics.m');
%%
% postfix = '_0.1';
% titlepost = 'at z = 10% PBL';
% zfactor = .1;
% fprintf('Analyzing Spectra...\n')
% run('analysis_spectra.m');
% 
% postfix = '_0.5';
% titlepost = 'at z = 50% PBL';
% zfactor = .5;
% fprintf('Analyzing Spectra...\n')
% run('analysis_spectra.m');
% 
% postfix = '_0.8';
% titlepost = 'at z = 80% PBL';
% zfactor = .8;
% fprintf('Analyzing Spectra...\n')
% run('analysis_spectra.m');
% 
