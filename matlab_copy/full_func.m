function full_func(iDir,iD)
datapathsuf = 'out/';
datapathpre = ['0708_bomex/run2/',iD{1}];
outputsuf   = 'plots/';
datapathlist = dir(datapathpre);

dataprefix = 'cm1out_';

numframes = 193;
titlepost = '';

numframes = 193;
if ~(datapathlist(iDir).isdir)
        return
end
datapath = fullfile(datapathpre,datapathlist(iDir).name,datapathsuf);
outputpath = fullfile(datapathpre,datapathlist(iDir).name,outputsuf);
if ~isfolder(datapath) | strcmp(datapath(1),'.')
        return
end

if false
try
z_spectra = 1000;
zfactor = 0;
budget_spectra = false;
compute_3d = false;
postfix = 'init';
titlepost = 't=0 at z = 1km';
k_spectra = 1;
fprintf('Analyzing Spectra...\n')
run('analysis_spectra.m');
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end
end
if true
try
        z_spectra = 1000;
        zfactor = 0;
        for iN=[73,numframes]
        postfix = num2str(iN);
        titlepost = ['t=',num2str(min(iN,73)*5/60 + max([(iN-73)*.5/60,0]),2),'hr at z = 1km'];
        k_spectra = iN:iN;
        budget_spectra = false;
        compute_3d = false;
        if iN == 73
                compute_3d = true;
        end
        if iN == numframes
          compute_3d = true;
        end
        fprintf('Analyzing Spectra...%i\n',iN)
        run('analysis_spectra.m');
        end
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end
end
if true
try
postfix = '';
fprintf('Making video...\n')
run('analysis_videography.m');
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end
end

if false
try
postfix = '';
titlepost = '';
fprintf('Analyzing Energy...\n')
run('analysis_energetics.m');
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end
end

try
postfix = '';
titlepost = '';
fprintf('Analyzing Clouds...\n');
run('analysis_cloudvid.m');
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end


try
numframes = 73;
postfix = '';
titlepost = '';
fprintf('Analyzing Clouds...\n');
run('analysis_clouds.m');
catch me
        disp(getReport(me,'extended','hyperlinks','off'))
end
end
