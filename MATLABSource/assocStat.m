function assocStat

load pigMyocardium.mat

[preLab, postLab] = getLabels;

pooLab = [preLab; postLab];

% Initialize structs for analysis and visualization (Plotly)

IR = struct();
MOLLI = struct();
SHMOLLI = struct();
SASHA = struct();
MTR = struct();
T2 = struct();

IR.prefix = getPrefixHearts(pigMyocardium,'IRvec');
MOLLI.prefix.Bias = getPrefixHearts(pigMyocardium,'MOLLIvec')- IR.prefix;
SHMOLLI.prefix.Bias = getPrefixHearts(pigMyocardium,'SHMOLLIvec') - IR.prefix;
SASHA.prefix.Bias = getPrefixHearts(pigMyocardium,'SASHAvec') - IR.prefix;
MTR.prefix =  getPrefixHearts(pigMyocardium,'MTRvec');
T2.prefix =  getPrefixHearts(pigMyocardium,'T2vec');

IR.postfix = getPostfixHearts(pigMyocardium,'IRvec');
MOLLI.postfix.Bias = getPostfixHearts(pigMyocardium,'MOLLIvec') - IR.postfix;
SHMOLLI.postfix.Bias = getPostfixHearts(pigMyocardium,'SHMOLLIvec') - IR.postfix;
SASHA.postfix.Bias = getPostfixHearts(pigMyocardium,'SASHAvec') - IR.postfix;
MTR.postfix =  getPostfixHearts(pigMyocardium,'MTRvec');
T2.postfix =  getPostfixHearts(pigMyocardium,'T2vec');

MOLLI.pool.Bias = [MOLLI.prefix.Bias MOLLI.postfix.Bias];
SHMOLLI.pool.Bias = [SHMOLLI.prefix.Bias SHMOLLI.postfix.Bias];
SASHA.pool.Bias = [SASHA.prefix.Bias SASHA.postfix.Bias];
T2.pool = [T2.prefix T2.postfix];
MTR.pool = [MTR.prefix MTR.postfix];


fixState = [{'prefix'},{'postfix'}];
label = [preLab postLab];

% T2
disp('===============================================');

% T2 pre & post -----------------------------------------------------

for i = 1:length(fixState)
    Lab = label(:,2*i-1:2*i);
    
    % MOLLI
    [MOLLI.(fixState{i}).T2.outlierBin,MOLLI.(fixState{i}).T2.elX, MOLLI.(fixState{i}).T2.elY, MOLLI.(fixState{i}).T2.CI, MOLLI.(fixState{i}).T2.pval] = ...
        runSkipCor(MOLLI.(fixState{i}).Bias, T2.(fixState{i}), 'MOLLI', 'T2', fixState{i}, Lab);
    
    % SHMOLLI
    [SHMOLLI.(fixState{i}).T2.outlierBin,SHMOLLI.(fixState{i}).T2.elX, SHMOLLI.(fixState{i}).T2.elY,SHMOLLI.(fixState{i}).T2.CI, SHMOLLI.(fixState{i}).T2.pval] = ...
        runSkipCor(SHMOLLI.(fixState{i}).Bias, T2.(fixState{i}), 'SHMOLLI', 'T2', fixState{i}, Lab);
    
    % SASHA
    [SASHA.(fixState{i}).T2.outlierBin,SASHA.(fixState{i}).T2.elX, SASHA.(fixState{i}).T2.elY, SASHA.(fixState{i}).T2.CI, SASHA.(fixState{i}).T2.pval] = ...
        runSkipCor(SASHA.(fixState{i}).Bias, T2.(fixState{i}), 'SASHA', 'T2', fixState{i}, Lab);
end

% MTR pre & post -----------------------------------------------------

for i = 1:length(fixState)
    Lab = label(:,2*i-1:2*i);
    
    % MOLLI
    [MOLLI.(fixState{i}).MTR.outlierBin,MOLLI.(fixState{i}).MTR.elX, MOLLI.(fixState{i}).MTR.elY, MOLLI.(fixState{i}).MTR.CI, MOLLI.(fixState{i}).MTR.pval] = ...
        runSkipCor(MOLLI.(fixState{i}).Bias, MTR.(fixState{i}), 'MOLLI', 'MTR', fixState{i}, Lab);
    
    % SHMOLLI
    [SHMOLLI.(fixState{i}).MTR.outlierBin,SHMOLLI.(fixState{i}).MTR.elX, SHMOLLI.(fixState{i}).MTR.elY, SHMOLLI.(fixState{i}).MTR.CI, SHMOLLI.(fixState{i}).MTR.pval] = ...
        runSkipCor(SHMOLLI.(fixState{i}).Bias, MTR.(fixState{i}), 'SHMOLLI', 'MTR', fixState{i}, Lab);
    
    % SASHA
    [SASHA.(fixState{i}).MTR.outlierBin,SASHA.(fixState{i}).MTR.elX, SASHA.(fixState{i}).MTR.elY, SASHA.(fixState{i}).MTR.CI, SASHA.(fixState{i}).MTR.pval] = ...
        runSkipCor(SASHA.(fixState{i}).Bias, MTR.(fixState{i}), 'SASHA', 'MTR', fixState{i}, Lab);
end


% Pooled T2 -----------------------------------------------------

fixState = 'pool';

% MOLLI
[MOLLI.(fixState).T2.outlierBin,MOLLI.(fixState).T2.elX, MOLLI.(fixState).T2.elY, MOLLI.(fixState).T2.CI, MOLLI.(fixState).T2.pval] = ...
    runSkipCor(MOLLI.(fixState).Bias, T2.(fixState), 'MOLLI', 'T2', fixState, pooLab);

% SHMOLLI
[SHMOLLI.(fixState).T2.outlierBin,SHMOLLI.(fixState).T2.elX, SHMOLLI.(fixState).T2.elY, SHMOLLI.(fixState).T2.CI, SHMOLLI.(fixState).T2.pval] = ...
    runSkipCor(SHMOLLI.(fixState).Bias, T2.(fixState), 'SHMOLLI', 'T2', fixState, pooLab);

% SASHA
[SASHA.(fixState).T2.outlierBin,SASHA.(fixState).T2.elX, SASHA.(fixState).T2.elY, SASHA.(fixState).T2.CI, SASHA.(fixState).T2.pval] = ...
    runSkipCor(SASHA.(fixState).Bias, T2.(fixState), 'SASHA', 'T2', fixState, pooLab);

% Pooled MTR -----------------------------------------------------

% MOLLI
[MOLLI.(fixState).MTR.outlierBin,MOLLI.(fixState).MTR.elX, MOLLI.(fixState).MTR.elY, MOLLI.(fixState).MTR.CI, MOLLI.(fixState).MTR.pval] = ...
    runSkipCor(MOLLI.(fixState).Bias, MTR.(fixState), 'MOLLI', 'MTR', fixState, pooLab);

% SHMOLLI
[SHMOLLI.(fixState).MTR.outlierBin,SHMOLLI.(fixState).MTR.elX, SHMOLLI.(fixState).MTR.elY, SHMOLLI.(fixState).MTR.CI, SHMOLLI.(fixState).MTR.pval] = ...
    runSkipCor(SHMOLLI.(fixState).Bias, MTR.(fixState), 'SHMOLLI', 'MTR', fixState, pooLab);

% SASHA
[SASHA.(fixState).MTR.outlierBin,SASHA.(fixState).MTR.elX, SASHA.(fixState).MTR.elY, SASHA.(fixState).MTR.CI, SASHA.(fixState).MTR.pval] = ...
    runSkipCor(SASHA.(fixState).Bias, MTR.(fixState), 'SASHA', 'MTR', fixState, pooLab);

disp('===============================================');

save exvivoCorrelation.mat IR T2 MTR MOLLI SHMOLLI SASHA preLab postLab pooLab

end

function [preHr] = getPrefixHearts(pigMyocardium, method)
preHr = [];
for hr_iter = 1:6
    for wk_iter = 1:3
        
        eval(['preHr = [preHr ' 'mean(pigMyocardium(' num2str(hr_iter) ').Week' num2str(wk_iter) '.' method ')];'])
        
    end
end


end


function [postHr] = getPostfixHearts(pigMyocardium, method)
postHr = [];
for hr_iter = 1:6
    for wk_iter = 4:6
        
        eval(['postHr = [postHr ' 'mean(pigMyocardium(' num2str(hr_iter) ').Week' num2str(wk_iter) '.' method ')];'])
        
    end
end

end

function [preLab, postLab] = getLabels

preLab = cell(18,2);
postLab = preLab;

k = 1;
for hr_iter = 1:6
    for wk_iter = 1:3
        
        preLab(k,1) = {['Heart' num2str(hr_iter)]};
        preLab(k,2) = {['Week' num2str(wk_iter)]};
        
        k = k+1;
    end
end

k = 1;
for hr_iter = 1:6
    for wk_iter = 4:6
        
        postLab(k,1) =   {['Heart' num2str(hr_iter)]};
        postLab(k,2) = {['Week' num2str(wk_iter)]};
        
        k = k+1;
        
    end
end



end

function [outlierBin,elX, elY, ciOut, pval] = runSkipCor(input2, input1, nm1, nm2, prepost, label, figFlag)

% Input order switched for input2 and input1 at function definition
% to keep T2 and MTR on the horizontal axis.

disp('--------------------------------------');
disp(['Skipped correlation ' prepost ' ' nm1  ' vs ' nm2]);


[~,~,~,outid] = skipped_correlation(input1',input2');
fgr = gcf;
saveas(fgr, [prepost '_' nm1 'vs' nm2 '.png']);
close(gcf);


outlierBin = zeros(length(input1),1);

if  not(isempty(outid{1}))
    
    
    disp([label(outid{1}',2) label(outid{1}',1)]);
    outlierBin(outid{1}) = 1;
    
else
    disp('No outliers can be detected for this distribution.');
end


% These are assigned into base by skipped_correlation.m

elX = evalin('base', 'ellipseX');
elY = evalin('base', 'ellipseY');
ciOut = evalin('base','ciOut'); 
pval = evalin('base','pval'); 

disp('--------------------------------------');
end


