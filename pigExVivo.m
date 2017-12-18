% -------------------------------------------------------------------------
% Please place this m-file in an empty folder and execute it.
% -------------------------------------------------------------------------
%
% Input : Data will be downloaded from OSF.
%
% Output:
%       - phantomData.mat
%       - pigMyocardium.mat
%       - pigExVivo_Maps.mat
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
% ========================================================================
% Download data from Open Science Framework
% Extract zip files

try
    disp('Please wait. Downloading data...');
    disp('Downloading heart dataset...')
    filename_heart = 'pigExVivo_Map.zip';
    url_heart = 'https://osf.io/6nt7e/download/';
    websave(filename_heart,url_heart);
    disp('Downloading phantom dataset...')
    filename_phantom = 'Phantom.zip';
    url_phantom = 'https://osf.io/f2ycr/download/';
    path = websave(filename_phantom,url_phantom);
    disp('Data has been downloaded successfully...');
catch
    warning('Data could not be downloaded.');
    warning('Please visit the OSF and download manually');
    warning('Remaning section can be executed after asigning path.');
    error('Cannot download data.');
end

% path = .../ if data cannot be downloaded automatically assign this var.
idx = max(strfind(path,filesep));
path = path(1:idx-1);
cd(path); % This is the main folder created by user.
mkdir('Outputs');
mkdir('Phantom');
disp('Extracting files...');
unzip(filename_heart);
cd('Phantom');
unzip([path filesep filename_phantom]);
cd(path);
% Remove zip files.
delete('Phantom.zip');
delete('pigExVivo_Map.zip');
% Display data tree
inf = imread([path filesep 'dataSchematic.png']);
imshow(inf);


%% ========================================================================
% Read heart data into structure - Type I
% This section will generate two structs:
%
%     - preFix(3X1 struct #: Scan1, Scan2, Scan3).methodName.Batch#(1,2,3)
%     - postFix(3X1 struct #: Scan4, Scan5, Scan6).methodName.Batch#(1,2,3)
% -------------------------------------------------------------------------

disp('Reading pigExVivo Maps Data: Type I...');


preFix = struct();
postFix = struct();

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
cd([path filesep 'Outputs']);
save pigExVivo_Maps.mat preFix postFix
%% ========================================================================
% Read phantom data into structure
%
% This section will generate one struct:
%
%     - phantomData(mX1 struct #: Number of selected ROIs).methodName_vec
%     ==> These entries contain masked values for selected method and ROI.
% -------------------------------------------------------------------------

disp('Reading Phantom Maps Data');

cd([path filesep 'Phantom']);

selectedROI = [8 10 11 13 14 16 17]; % You can modify here (1 to 18 ROI)
tmpName = 'MT-ROI_';
cd('MASKS');
% Read all masks into a matrix.
tempMask = getMatData('MT-ROI_1.mat');
sz = size(tempMask);
phantomMask = zeros(sz(1),sz(2),length(selectedROI));

for i =1:length(selectedROI)
    phantomMask(:,:,i) = getMatData([tmpName num2str(selectedROI(i)) '.mat']);
    
end

cd([path filesep 'Phantom']);

phantomData = struct();

phanMask = logical(phantomMask);

imIR = getMatData('IR_Phantom.mat');
imMTR = getMatData('MTR_Phantom.mat');
imMOLLI = getMatData('MOLLI_Phantom.mat');
imSHMOLLI = getMatData('SHMOLLI_Phantom.mat');
imSASHA = getMatData('SASHA_Phantom.mat');
imT2 = getMatData('T2_Phantom.mat');

for i =1:length(selectedROI)
    
    phantomData(i).IRvec = imIR(phanMask(:,:,i));
    phantomData(i).MTRvec = imMTR(phanMask(:,:,i));
    phantomData(i).MOLLIvec = imMOLLI(phanMask(:,:,i));
    phantomData(i).SHMOLLIvec = imSHMOLLI(phanMask(:,:,i));
    phantomData(i).SASHAvec = imSASHA(phanMask(:,:,i));
    phantomData(i).T2vec = imT2(phanMask(:,:,i));
    
end


cd([path filesep 'Outputs']);

save phantomData.mat phantomData imIR imMTR imMOLLI imSHMOLLI imSASHA imT2 phanMask

%% ========================================================================
% Read heart data into structure - Type II & Statistical Operations
%
% This section will generate one struct:
%
%     - pigMyocardium(Heart#: 1 to 6).Week#(#: 1 to 6).methodName_vec
%     ==> Contains masked values for selected week and heart.
%
%
% Jupyter Notebook inputs this structures.
% -------------------------------------------------------------------------

disp('Reading pigExVivo Maps Data: Type II...');

cd(path);

pigMyocardium  = struct();
assignin('base','pigMyocardium',pigMyocardium);

for k=1:2
    
    if k==1
        curDat = preFix;
    elseif k==2
        curDat = postFix;
    end
    
    
    
    for j = 1:3
        
        maskHeart(curDat,j,k,'MTR');
        maskHeart(curDat,j,k,'IR');
        maskHeart(curDat,j,k,'T2');
        maskHeart(curDat,j,k,'MOLLI');
        maskHeart(curDat,j,k,'SHMOLLI');
        maskHeart(curDat,j,k,'SASHA');
        
    end
    
    
    
end

pigMyocardium = evalin('base', 'pigMyocardium');

cd([path filesep 'Outputs']);
save pigMyocardium.mat pigMyocardium

disp(['All outputs have been saved to ' path filesep 'Outputs']);
disp('---------------------------------------------- DONE.');

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

if strcmp(fields{1},'nlvolume')
    
    nlvolume = str.nlvolume;
    curDat = double(nlvolume.theData);
    
    if length(size(curDat))>2
        
        im = squeeze(curDat(:,:,1,1));
        
    else
        
        im = curDat;
    end
    
else
    
    im = double(str.(fields{1}));
    
end

end

function maskHeart(curDat,j,k,mapName)

if k==1
    wk = j;
elseif k==2
    wk = j+3;
end

pigMyocardium = evalin('base', 'pigMyocardium');
lut = [1,1,2,2,3,3];
for k=1:6
    msk = curDat(j).Mask(:,:,k);
    im = curDat(j).(mapName).(['Batch' num2str(lut(k))]);
    pigMyocardium(k).(['Week' num2str(wk)]).([mapName 'vec']) = im(logical(msk));
end


assignin('base','pigMyocardium',pigMyocardium);


end




