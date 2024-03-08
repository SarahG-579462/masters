names =  {'default','dry'}

for iDd =1:length(names)
iD = names(iDd);
datapathsuf = 'out/';
datapathpre = ['0708_bomex/run2/',iD{1}];
outputsuf   = 'plots/';
datapathlist = dir(datapathpre);

dataprefix = 'cm1out_';

numframes = 193;
titlepost = '';
parpool(10);
parfor iDir=1:length(datapathlist)
  full_func(iDir,iD);
end
end
