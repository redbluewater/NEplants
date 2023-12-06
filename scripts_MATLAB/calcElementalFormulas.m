%calculate elemental formulas for the Bowen/NEplants2 dataset
%Remeber this will require older MATLAB (2017b will work fine)
%Krista Longnecker, 13 November 2023
clear all
close all

%set the path to the scripts, note these scripts are available on GitHub 
%here: https://github.com/redbluewater/findformula
path(path,'C:\Users\klongnecker\Documents\Current projects\MSdataAnalysis\listOrgCpds');

load NEplants_neg_aligned.2023.12.06.mat

%Need to prune out the isomers as the script assumes all the mz values are
%unique, which may or may not be true. Turns out that in this dataset,
%there are no isomers, which is handy. Leave the error check here.
uPeaks = unique(Peaks);
if ~isequal(length(uPeaks),length(Peaks))
    error('need to remove isomers first before proceeding')
end
clear uPeaks

%Convert the peaks to neutral mass
H = 1.007825032;
elec = 5.4858e-4;
IUPACpeak = Peaks + H - elec; % for negative ion mode
clear H elec

load LongneckerKujawinski_fullCompoundList.2016.11.21.mat
[formulas elementOrder] = findformula_useList_KL17(IUPACpeak, zeros(size(IUPACpeak)), 1, 20, 500,fullCompoundList,'HAcap',1);

%go through and find the C13-containing formulas
FormulasC13 = quick13C_KL_1(formulas,IUPACpeak,0);

%calc DBE, AI, AImod, NOSC
[DBE AI AImod NOSC] = calcDBEandAI_3(FormulasC13);

%Matt Costa is an R user, so save the result as a CSV file so he can open
%the file; 
csvForm = array2table(FormulasC13);
csvForm.Properties.VariableNames = elementOrder;
csvForm.DBE = DBE;
csvForm.AI = AI;
csvForm.AImod = AImod;
csvForm.NOSC = NOSC;

writetable(csvForm,'NEplants2_formulas.2023.12.06.csv');


