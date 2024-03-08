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

c = mfilename
b = regexp(c,'\d+','match');
iDir = str2num(b{1});
full_func(2+iDir,iD);
end
