% June 19, 2016
% script to format data to run ITC hyperbolic model

clear all

% read in behavioral data from csv file
cd('/data/joy/BBL/projects/pehlivanovaPncItc/pehlivanovaPncItcScripts/itcScripts')

importingItcData

% compiling list of relevant delay discounting variables
% combine individual question variables; invert the matrix and change
% options form 1=chose immediate option / 2=chose delayed option to 0/1 to
% match script format

choicesITC=[kddisc_q_01 kddisc_q_02 kddisc_q_03 kddisc_q_04 kddisc_q_05 ...
    kddisc_q_06 kddisc_q_07 kddisc_q_08 kddisc_q_09 kddisc_q_10 ...
    kddisc_q_11 kddisc_q_12 kddisc_q_13 kddisc_q_14 kddisc_q_15 ...
    kddisc_q_16 kddisc_q_17 kddisc_q_18 kddisc_q_19 kddisc_q_20 ...
    kddisc_q_21 kddisc_q_22 kddisc_q_23 kddisc_q_24 kddisc_q_25 ...
    kddisc_q_26 kddisc_q_27 kddisc_q_28 kddisc_q_29 kddisc_q_30 ...
    kddisc_q_31 kddisc_q_32 kddisc_q_33 kddisc_q_34]'-1;

rtsITC=[kddisc_trr_01 kddisc_trr_02 kddisc_trr_03 kddisc_trr_04 kddisc_trr_05 ...
    kddisc_trr_06 kddisc_trr_07 kddisc_trr_08 kddisc_trr_09 kddisc_trr_10 ...
    kddisc_trr_11 kddisc_trr_12 kddisc_trr_13 kddisc_trr_14 kddisc_trr_15 ...
    kddisc_trr_16 kddisc_trr_17 kddisc_trr_18 kddisc_trr_19 kddisc_trr_20 ...
    kddisc_trr_21 kddisc_trr_22 kddisc_trr_23 kddisc_trr_24 kddisc_trr_25 ...
    kddisc_trr_26 kddisc_trr_27 kddisc_trr_28 kddisc_trr_29 kddisc_trr_30 ...
    kddisc_trr_31 kddisc_trr_32 kddisc_trr_33 kddisc_trr_34]';

%% save data to use later
% save('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/rawItcAug2016.mat','bblid','choicesITC','rtsITC')

% there are no missing values here because of how the original file was
% constructed

%% running hyperbolic delay discounting model to obtain k-values
% load data about amounts and delays for each question
load('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/itemOrderITC.mat')

clearvars -except bblid choicesITC rtsITC itemOrderITC

k_new=nan(size(choicesITC,2),1); % k values
rsq_new=nan(size(choicesITC,2),1); % R2 from hyperbolic model
corr1_new=nan(size(choicesITC,2),1); % correlation between RT and SV
medRT_new=nan(size(choicesITC,2),1); % median reaction time
pctPred_new=nan(size(choicesITC,2),1);
noise_new=nan(size(choicesITC,2),1);
later=nan(size(choicesITC,2),1); % how many delayed choices
rsqMRN=nan(size(choicesITC,2),1); % R2 from logistic model
tjur=nan(size(choicesITC,2),1); % tjur's coefficient of discrimination

tic
for i=1:size(choicesITC,2)
    
i
    if ~isnan(choicesITC(1,i))
    
    % old script, summer 2014    
    %[hyperbolic_old] = ITCanalysis(choicesITC(:,i),itemOrderITC(:,2),itemOrderITC(:,3),...
    %       itemOrderITC(:,4),itemOrderITC(:,5),rtsITC(:,i),bblid(i))
        
    % new script from Arthur modified by Becca, summer 2016
    [hyperbolic_new] = ITCScreenAnalysis_mp(choicesITC(:,i),itemOrderITC(:,2),...
            itemOrderITC(:,4),itemOrderITC(:,5),rtsITC(:,i))
               
    k_new(i,1)=hyperbolic_new.k;
    rsq_new(i,1)=hyperbolic_new.r2;
    corr1_new(i,1)=hyperbolic_new.RTandSubjValueCorr;
    medRT_new(i,1)=hyperbolic_new.medianRT;
    pctPred_new(i,1)=hyperbolic_new.percentPredicted;
    noise_new(i,1)=hyperbolic_new.noise;
    rsqMRN(i,1)=hyperbolic_new.R2mnrfit;
    later(i,1)=sum(choicesITC(:,i));
    tjur(i,1)=hyperbolic_new.tjur;
    
    end
    
end
toc

% writing out ITC results in a table to export into a csv
row_names =  {'bblid'; 'k'; 'logk'; 'kHypR2'; 'kMrnR2'; 'kTjur'; 'kCorr'; 'kMedRT'; 'kNoise'; ...
    'kPctPred'; 'kDelayedCh'}

itc_all=table(bblid, k_new,log(k_new), rsq_new, rsqMRN, tjur, corr1_new, medRT_new, ...
    noise_new, pctPred_new, later, 'VariableNames', row_names)

% saving to a csv file for R 
writetable(itc_all,'/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/n453_Kvalues_08162016.csv','Delimiter',',')

% saving in matlab format, just in case
save('/data/joy/BBL/projects/pehlivanovaPncItc/subjectData/itcRiskData/itcAug152016.mat', ...
    'k_new','rsq_new', 'rsqMRN','tjur', 'noise_new','medRT_new','pctPred_new','corr1_new')





   













