% -------------------------------------------------------------------------
% Please place this m-file in an empty folder and execute it.
% -------------------------------------------------------------------------
%
% Input : Data will be downloaded from OSF.
%
% Output:
%       - pigExVivo_Maps.mat (contains preFix and postFix structures)
%       - pigExVivo_Stat.mat (contains)
%
% This script follows the variable naming and structure convention 
% that is compatible with following (online executable) Jupyter Notebook:
%
% http://rebrand.ly/pigExVivo
% 
%
% Written by: Agah Karakuzu 
% =========================================================================

function pigExVivo

try
disp('Please wait. Downloading data...');    
filename = 'pigExVivo_Map.zip';
url = 'https://osf.io/6nt7e/download/';
path = websave(filename,url);
disp('Data has been downloaded.');    
catch
warning('Data could not be downloaded.');
warning('Please visit the OSF and download manually');
warning('Remaning section can be executed after asigning path.');
error('Cannot download data.');
end

% path = .../ if data cannot be downloaded automatically assign this var.
idx = max(strfind(path,filesep));
path = path(1:idx-1);
cd(path);
unzip(filename);

inf = imread([path filesep 'dataSchematic.png']);
imshow(inf);

preFix = struct();
postFix = struct();

% Pre-fixation scans: SCAN01 SCAN02 SCAN03 
% Post-fixation scans: SCAN04 SCAN05 SCAN06

for i=1:3
preFix(i).scanName  = ['SCAN0' num2str(i)];
postFix(i).scanName = ['SCAN0' num2str(i+3)] ;
end

for i =1:3 % for each pre/post-fixation scan 
    
    for j=1:3 % for each batch 
    
    % Read preFixation data
    
    preFix(i).MTR.(['Batch' num2str(j)]) = getData('MTR',preFix(i).scanName,j,'pre',path);
    preFix(i).IR.(['Batch' num2str(j)]) = getData('IR',preFix(i).scanName,j,'pre',path);
    preFix(i).MOLLI.(['Batch' num2str(j)]) = getData('MOLLI',preFix(i).scanName,j,'pre',path);
    preFix(i).SHMOLLI.(['Batch' num2str(j)]) = getData('SHMOLLI',preFix(i).scanName,j,'pre',path);
    preFix(i).SASHA.(['Batch' num2str(j)]) = getData('SASHA',preFix(i).scanName,j,'pre',path);
    preFix(i).T2.(['Batch' num2str(j)]) = getData('T2',preFix(i).scanName,j,'pre',path);

    
    % Read postFixation data
    
    postFix(i).MTR.(['Batch' num2str(j)]) = getData('MTR',postFix(i).scanName,j,'post',path);
    postFix(i).IR.(['Batch' num2str(j)]) = getData('IR',postFix(i).scanName,j,'post',path);
    postFix(i).MOLLI.(['Batch' num2str(j)]) = getData('MOLLI',postFix(i).scanName,j,'post',path);
    postFix(i).SHMOLLI.(['Batch' num2str(j)]) = getData('SHMOLLI',postFix(i).scanName,j,'post',path);
    postFix(i).SASHA.(['Batch' num2str(j)]) = getData('SASHA',postFix(i).scanName,j,'post',path);
    postFix(i).T2.(['Batch' num2str(j)]) = getData('T2',postFix(i).scanName,j,'post',path);
    
    
    end
    
    % Read masks (m x n x 6)
    
    scnID = preFix(i).scanName;
    preFix(i).Mask = getMatData([path filesep 'preFixation'  filesep scnID filesep 'MASKS' filesep scnID '_MASKS.mat']);
    
    scnID = postFix(i).scanName;
    postFix(i).Mask = getMatData([path filesep 'postFixation'  filesep scnID filesep 'MASKS' filesep scnID '_MASKS.mat']);
    
    
    
end

 save pigExVivo_Maps.mat preFix postFix


     
end   

function out = getData(modName,scanName,j,fixState,path)


mainDir = path;
if strcmp(fixState,'pre')
fixDir = [mainDir filesep 'preFixation'];
elseif strcmp(fixState,'post')
fixDir = [mainDir filesep 'postFixation'];
end

switch modName
    case 'MTR'
folName = 'MTR';
    case 'IR'
folName = 'T1_IR';
    case 'MOLLI'
folName = 'T1_MOLLI';
    case 'SASHA'
folName = 'T1_SASHA';
    case 'SHMOLLI'
folName = 'T1_SHMOLLI';
    case 'T2'
folName = 'T2_SE_MC';
end

scanFol = [fixDir filesep scanName];
curName = [modName '_' scanName 'BATCH0' num2str(j) '.mat'];
out = getMatData([scanFol filesep folName filesep curName]);


end
     
function im = getMatData(filename)

str = load(filename);
fields = fieldnames(str);

   
im = double(str.(fields{1}));    



end
     
     