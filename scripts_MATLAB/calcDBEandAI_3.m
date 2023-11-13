function [DBE AI AImod NOSC] = calcDBEandAI_3(formulas);
%function [DBE AI AImod NOSC] = calcDBEandAI_3(formulas);
%do the additional calculations of DBE and AI
%input
%'formulas'
% should be an 9-column matrix in the following order: 
% {'C','H','O','N','C13','S','P','Na','Error';}
%output
% 'DBE' is the double-bond equivalent calculated as:
% DBE = 1 + 1/2(2C - H + N + P). Original reference is 
%McClafferty and Turecek, 1993)
% 'AI' from Koch and Dittmar 2006 which is:
%%corrected 12/31/2015
% AI = (1 + C - O - S - 0.5(N+P+H))/(C - O - S - N - P)...
% the full description of AI is the following: AI = DBEai / Cai where
%DBEai = 1 + C - O - S - 0.5(N+P+H)
%Cai = C - O - S - N - P
%and then the added note that AI = 0 if either DBEai or Cai is <= 0
%AImod 'possible number of carbonyl substitutions is halved' (quoting
%Stubbins et al. L&O 2010); corrected per 2016 erratum
% AImod = (1 + C - 0.5*O - S - 0.5(N+P+H)/(C - 0.5*O - S - N - P)...
%KL 9/2/10
%adding NOSC (KL 12/7/2011)
%NOSC = 4 - [(4c + h - 3n - 2o - 2s) / c]
% updated 12/31/2015 to reflect 2016 erratum published on the AI
% calculation

%first calculate DBE
DBE = 1 + 0.5*(2*(formulas(:,1) + formulas(:,5)) - formulas(:,2) + formulas(:,4) + formulas(:,7));
%but if there is no elemental formula, this will spit back a one...
k = find(formulas(:,1)==0); % if no carbon, no formula -> use that to set those DBE to zero
DBE(k) = 0;
clear k 

%and also want AI:
num = 1 + (formulas(:,1) + formulas(:,5)) - formulas(:,3) - formulas(:,6) - ...
    0.5*(formulas(:,2) + formulas(:,7) + formulas(:,4));

den = (formulas(:,1) + formulas(:,5)) - formulas(:,3) - formulas(:,6) - formulas(:,4) - formulas(:,7);
k = find(num > 0 & den > 0);
AI = zeros(size(formulas,1),1);
AI(k) = num(k)./den(k);
clear num den

%and the modified AI:
num = 1 + (formulas(:,1) + formulas(:,5)) - 0.5*formulas(:,3) - formulas(:,6) - ...
    0.5*(formulas(:,2) + formulas(:,7) + formulas(:,4));

den = (formulas(:,1) + formulas(:,5)) - 0.5*formulas(:,3) - formulas(:,6) - formulas(:,4) - formulas(:,7);
k = find(num > 0 & den > 0);
AImod = zeros(size(formulas,1),1);
AImod(k) = num(k)./den(k);
clear num den

%and the NOSC
%NOSC = 4 - [(4c + h - 3n - 2o - 2s) / c]
num = (4*formulas(:,1) + formulas(:,2) - 3*formulas(:,4) - 2*formulas(:,3) - 2*formulas(:,6));
den = formulas(:,1);
NOSC = find(num > 0 & den > 0);
NOSC = zeros(size(formulas,1),1);
NOSC(k) = 4 - num(k)./den(k);
clear num den



