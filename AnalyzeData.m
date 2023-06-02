% Data Analysis for Reeves & Otero-Millan 2023

%% Load invididual data for individual sessions and add and ID column to each file
subjects = {'0101', '0122', '0137', '0138', '0139', '0140', '0141', '0169', '0174', '0180', '0181', '0182', '0183', '0184', '0185', '0186', '0187', '0188', '0189', '0190'};
subjectsFull = {'0101__A', '0122__A', '0137__A', '0138__A', '0139__A', '0140__A', '0141__A', '0169__A', '0174__A', '0180__A', '0181__A', '0182__A', '0183__A', '0184__A', '0185__A', '0186__A', '0187__A', '0188__A', '0189__A', '0190__A'};
SubjList = ["0101","0122","0137","0138","0139","0140","0141","0169","0174","0180","0181","0182","0183","0184","0185","0186","0187","0188","0189","0190"];
numSubjList = [101, 122, 137, 138, 139, 140, 141, 169, 174, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191]'; % should match above line
% projectFolder = 'D:\OneDrive\UC Berkeley\OMlab - OM-lab-share\Projects\MicrosaccadesImageTilt\Data\MicrosaccImTilt_FullStudy';
%projectFolder = 'C:\Users\stephanie_reeves\UC Berkeley\OMlab - OM-lab-share\Projects\MicrosaccadesImageTilt\Data\MicrosaccImTilt_FullStudy';
projectFolder =  '/Users/stephaniereeves/UC Berkeley/OMlab - OM-lab-share/Projects/MicrosaccadesImageTilt/Data/MicrosaccImTilt_FullStudy'; %home

saccade = table();
trialTable = table();
numTrials = 240; % The number of trials for each subject

for i=1:length(subjects)
    qpfile = fullfile(projectFolder, ['FreeViewingFixation__' subjectsFull{i}], 'AnalysisResults_QuickPhases.mat');
    ttfile = fullfile(projectFolder, ['FreeViewingFixation__' subjectsFull{i}], 'trialDataTable.mat');
    stfile = fullfile(projectFolder, ['FreeViewingFixation__' subjectsFull{i}], 'samplesDataTable.mat');
    if (exist(qpfile,'file'))
        d = load(qpfile);
        qps = d.AnalysisResults_QuickPhases;
        qps.Subject = categorical(cellstr(repelem(subjects{i},height(qps),1)));
    else
        disp(['File not found ' qpfile]);
    end
    
    if (exist(ttfile,'file'))
        e = load(ttfile);
        tdt = e.trialDataTable;
        tdt.Subject = categorical(cellstr(repelem(subjects{i},height(tdt),1)));
    else
        disp(['File not found ' ttfile]);
    end

    if (exist(stfile,'file'))
        f = load(stfile);
        sdt = f.samplesDataTable;
        sdt.Subject = categorical(cellstr(repelem(subjects{i},height(sdt),1)));
        zeros(numTrials,1);
        
        % Get the avg position of the two eyes
        sdt.meanT = mean([sdt.RightT sdt.LeftT],2,'omitnan');
        medianT_AccountingForBaselineMeanT = []; 
        
        % Get the baseline X/Y/T (during the 2 s of fixation) to subtract
        % it off
        for j = 1:numTrials 
            start_index = find(sdt.TrialNumber==j & sdt.Subject == SubjList(i),1,'first');
            end_time = sdt.Time(start_index) + 2;
            end_index = find(sdt.Time == end_time & sdt.Subject == SubjList(i));
            medianT = median(sdt.meanT(start_index:end_index),'omitnan'); % first two seconds of the trial where they are fixating
            medianT_AccountingForBaselineMeanT = [medianT_AccountingForBaselineMeanT medianT];
        end
        tdt.medianT_AccountingForBaselineMeanT = tdt.median_T - medianT_AccountingForBaselineMeanT';
       

    else
        disp(['File not found ' stfile]);
    end
     
    saccade = vertcat(saccade, qps); % This is the table now that has all saccades for all subjects. You'll use this for the rest of everything.
    trialTable = vertcat(trialTable, tdt);
    
    % Track how far along we are and print it out
    disp(['We are ', num2str(i/length(SubjList)*100), ' percent done!']);
end
clear d e i qpfile qps ttfile tdt sdt;

% Make the first 2 seconds of each 12 second trial to be the "pre-trial"
% (this was just when the fixation dot appeared and no scene was present
saccade.Task(saccade.TimeFromTrialBegining < 2) = "Pre-Trial";

%% Get polar histograms (using KDE) and do circular cross-correlations 
% If DEBUG = 1, then plots will print out 
DEBUG = 0;

% Initializing things
Tasks = ["Fixation","FreeView"];
LRTilts = [-30,30];
xcorrKDE_Fix_0and30 = [];
xcorrKDE_Fix_0andneg30 = [];
xcorrKDE_Free_0and30 = [];
xcorrKDE_Free_0andneg30 = [];

% You'll get 4 vectors (2 tasks x 2 tilt comparisons)
for atask = 1:length(Tasks)
    
    for imtilt = 1:length(LRTilts)
        if (DEBUG)
            figure;
            sgtitle(sprintf('%s',LRTilts(imtilt), ' and 0 - ',Tasks(atask)));
        end
        
        for subj = 1:length(SubjList)
            
            % Get the saccade direction vector that you are interested in
            switch (true)
                case LRTilts(imtilt) == -30 & Tasks(atask) == "FreeView"
                    h1 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == LRTilts(imtilt) & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                case LRTilts(imtilt) == 30 & Tasks(atask) == "FreeView"
                    h1 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == LRTilts(imtilt) & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                case LRTilts(imtilt) == -30 & Tasks(atask) == "Fixation"
                    h1 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == LRTilts(imtilt) & saccade.Amplitude < 3 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == 0 & saccade.Amplitude < 3 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                case LRTilts(imtilt) == 30 & Tasks(atask) == "Fixation"
                    h1 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == LRTilts(imtilt) & saccade.Amplitude < 3 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Task == Tasks(atask) & saccade.ImTilt == 0 & saccade.Amplitude < 3 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            end
            
            
            % Calculate KDE for both distributions
            [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 radians is the width of the kernel
            [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1);
            
            % Do the circular cross correlation (circular bc we do the
            % repeating)
            hist1Rep = [vfEstimate vfEstimate vfEstimate];
            hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
            [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % only search within 45 degs
            [~,i] = max(x);
            maxlag = lags(i)/10;
            
            % Accumulate these values
            switch true
                case LRTilts(imtilt) == -30 && Tasks(atask) == "Fixation"
                    xcorrKDE_Fix_0andneg30 = [xcorrKDE_Fix_0andneg30 maxlag];
                case LRTilts(imtilt) == 30 && Tasks(atask) == "Fixation"
                    xcorrKDE_Fix_0and30 = [xcorrKDE_Fix_0and30 maxlag];
                case LRTilts(imtilt) == -30 && Tasks(atask) == "FreeView"
                    xcorrKDE_Free_0andneg30 = [xcorrKDE_Free_0andneg30 maxlag];
                case LRTilts(imtilt) == 30 && Tasks(atask) == "FreeView"
                    xcorrKDE_Free_0and30 = [xcorrKDE_Free_0and30 maxlag];
            end
            
            % Get visualization for each subject
            if (DEBUG)
                subplot(3,length(SubjList),subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; 
                polarplot([0:.1:360]/180*pi, vfEstimate2);
                title(sprintf('Subj %s',SubjList(subj)))
                subplot(3,length(SubjList),subj+length(SubjList)); plot(vfEstimate); hold on; 
                plot(vfEstimate2);
                subplot(3,length(SubjList),subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                    text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
            end
            
        end
    end
end

% Save the xcorr values in a table
tab = table();
tab.Subject = transpose(SubjList);
tab.xcorrKDE_Fix_0andneg30(:,1) = xcorrKDE_Fix_0andneg30;
tab.xcorrKDE_Fix_0and30(:,1) = xcorrKDE_Fix_0and30;
tab.xcorrKDE_Free_0andneg30(:,1) = xcorrKDE_Free_0andneg30;
tab.xcorrKDE_Free_0and30(:,1) = xcorrKDE_Free_0and30;


%% Bootstrapping
LRTilts = [-30,30];

% Do this for both tasks
for atask = 1:length(Tasks)
    
    % Do this for both -30 and 30 image tilts
    for atilt = 1:length(LRTilts)
        
        % Do the bootstrapping numSampling times
        numSampling = 1000;
        
        % Initialize conglomoration of displacements
        displacements = NaN([length(SubjList),numSampling]);
        
        % Get sacc dirs needed
        for subj = 1:length(SubjList)
            h1 = saccade.Direction(find(saccade.ImTilt == LRTilts(atilt) & saccade.Task == Tasks(atask) & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            h2 = saccade.Direction(find(saccade.ImTilt == 0 & saccade.Task == Tasks(atask) & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            
            % Get the sampling indicies
            idx_h1 = randi(length(h1),length(h1),numSampling);
            idx_h2 = randi(length(h2),length(h2),numSampling);
            
            % Sample with replacement numSampling times
            for col = 1:numSampling
                newh1 = h1(idx_h1(:,col));
                newh2 = h2(idx_h2(:,col));
                
                % Calculate KDE for both distributions
                [vfEstimate] = circ_ksdensity(newh1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 is arbitrary (it's the width of the kernel), but seems to work well
                [vfEstimate2] = circ_ksdensity(newh2,0:.1/180*pi:2*pi,[-pi pi],.1);
                
                % Do the circular cross correlation (circular bc we do the
                % repeating)
                hist1Rep = [vfEstimate vfEstimate vfEstimate];
                hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
                [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % 45 degs
                [~,i] = max(x);
                maxlag = lags(i)/10;
                
                % Conglomorate
                displacements(subj,col) = maxlag;
            end
        end
        
        
        % Make four tables
        switch (true)
            case LRTilts(atilt) == 30 && Tasks(atask) == "Fixation"
                bootstrapped_fix_right = table();
                bootstrapped_fix_right.Subject = transpose(SubjList);
                bootstrapped_fix_right.SubjNum(:,1) = [1:length(SubjList)];
                bootstrapped_fix_right.Avg(:,1) = 0;
                bootstrapped_fix_right.CI_low(:,1) = 0;
                bootstrapped_fix_right.CI_high(:,1) = 0;
                for subj = 1:length(SubjList)
                    bootstrapped_fix_right.Avg(subj,1) = mean(displacements(subj,:));
                    bootstrapped_fix_right.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                    bootstrapped_fix_right.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
                end
                
            case LRTilts(atilt) == -30 && Tasks(atask) == "Fixation"
                bootstrapped_fix_left = table();
                bootstrapped_fix_left.Subject = transpose(SubjList);
                bootstrapped_fix_left.SubjNum(:,1) = [1:length(SubjList)];
                bootstrapped_fix_left.Avg(:,1) = 0;
                bootstrapped_fix_left.CI_low(:,1) = 0;
                bootstrapped_fix_left.CI_high(:,1) = 0;
                for subj = 1:length(SubjList)
                    bootstrapped_fix_left.Avg(subj,1) = mean(displacements(subj,:))
                    bootstrapped_fix_left.CI_low(subj,1) = prctile(displacements(subj,:),2.5)
                    bootstrapped_fix_left.CI_high(subj,1) = prctile(displacements(subj,:),97.5)
                end
            case LRTilts(atilt) == 30 && Tasks(atask) == "FreeView"
                bootstrapped_free_right = table();
                bootstrapped_free_right.Subject = transpose(SubjList);
                bootstrapped_free_right.SubjNum(:,1) = [1:length(SubjList)];
                bootstrapped_free_right.Avg(:,1) = 0;
                bootstrapped_free_right.CI_low(:,1) = 0;
                bootstrapped_free_right.CI_high(:,1) = 0;
                for subj = 1:length(SubjList)
                    bootstrapped_free_right.Avg(subj,1) = mean(displacements(subj,:));
                    bootstrapped_free_right.CI_low(subj,1) = prctile(displacements(subj,:),2.5);
                    bootstrapped_free_right.CI_high(subj,1) = prctile(displacements(subj,:),97.5);
                end
                
            case LRTilts(atilt) == -30 && Tasks(atask) == "FreeView"
                bootstrapped_free_left = table();
                bootstrapped_free_left.Subject = transpose(SubjList);
                bootstrapped_free_left.SubjNum(:,1) = [1:length(SubjList)];
                bootstrapped_free_left.Avg(:,1) = 0;
                bootstrapped_free_left.CI_low(:,1) = 0;
                bootstrapped_free_left.CI_high(:,1) = 0;
                for subj = 1:length(SubjList)
                    bootstrapped_free_left.Avg(subj,1) = mean(displacements(subj,:))
                    bootstrapped_free_left.CI_low(subj,1) = prctile(displacements(subj,:),2.5)
                    bootstrapped_free_left.CI_high(subj,1) = prctile(displacements(subj,:),97.5)
                end
        end
    end
end

% Sort subjects by combined right/left effect for each task
% For fixation
this = (abs(bootstrapped_fix_left.Avg) + abs(bootstrapped_fix_right.Avg))/2;
[sure,idx] = sortrows(this,'descend')
bootstrapped_fix_right = bootstrapped_fix_right(idx,:);
bootstrapped_fix_left = bootstrapped_fix_left(idx,:);

% For free viewing
this = (abs(bootstrapped_free_left.Avg) + abs(bootstrapped_free_right.Avg))/2;
[sure,idx] = sortrows(this,'descend')
bootstrapped_free_right = bootstrapped_free_right(idx,:);
bootstrapped_free_left = bootstrapped_free_left(idx,:);

%% Looking at Free Viewing and saccades of different sizes -- using quartiles 
% If DEBUG = 1, then plots will print out 
DEBUG = 0;

% Initializing things
LRTilts = ["-30","30"];
Quartiles = ["First","Second","Third","Fourth"];
xcorrKDE_Free_0and30_firstqt = [];
xcorrKDE_Free_0and30_secondqt = [];
xcorrKDE_Free_0and30_thirdqt = [];
xcorrKDE_Free_0and30_fourthqt = [];
xcorrKDE_Free_0andneg30_firstqt = [];
xcorrKDE_Free_0andneg30_secondqt = [];
xcorrKDE_Free_0andneg30_thirdqt = [];
xcorrKDE_Free_0andneg30_fourthqt = [];

% Define amplitude cutoffs using quartile function 
Q_25 = quantile(saccade.Amplitude(find(saccade.Task == "FreeView")),0.25);
Q_50 = quantile(saccade.Amplitude(find(saccade.Task == "FreeView")),0.50);
Q_75 = quantile(saccade.Amplitude(find(saccade.Task == "FreeView")),0.75);


for aquartile = 1:length(Quartiles)
    
    for imtilt = 1:length(LRTilts)
        if (DEBUG)
            figure;
            sgtitle(sprintf('%s',LRTilts(imtilt), ' and 0 - ',Quartiles(aquartile)));
        end
        
        for subj = 1:length(SubjList)
            
            % Get the saccade direction vectors that you are interested in
            switch (true)
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "First"
                    h1 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Second"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Third"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Fourth"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "First"
                    h1 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Second"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Third"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Fourth"
                    h1 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
                    h2 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.Subject == SubjList(subj) & saccade.TrialNumber > 0));
            end
            
            % Calculate KDE for both distributions
            [vfEstimate] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1);  % The .1 radians is the width of the kernel
            [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1);
            
            % Do the circular cross correlation (circular bc we do the
            % repeating)
            hist1Rep = [vfEstimate vfEstimate vfEstimate];
            hist2Rep = [zeros(size(vfEstimate2)) vfEstimate2 zeros(size(vfEstimate2))];
            [x,lags] = xcorr(hist2Rep, hist1Rep, 450); % only search within 45 degs
            [howbig,i] = max(x);
            maxlag = lags(i)/10;
            
            % Accumulate these values
            switch true
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "First"
                    xcorrKDE_Free_0andneg30_firstqt = [xcorrKDE_Free_0andneg30_firstqt maxlag];
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Second"
                    xcorrKDE_Free_0andneg30_secondqt = [xcorrKDE_Free_0andneg30_secondqt maxlag];
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Third"
                    xcorrKDE_Free_0andneg30_thirdqt = [xcorrKDE_Free_0andneg30_thirdqt maxlag];
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Fourth"
                    xcorrKDE_Free_0andneg30_fourthqt = [xcorrKDE_Free_0andneg30_fourthqt maxlag];
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "First"
                    xcorrKDE_Free_0and30_firstqt = [xcorrKDE_Free_0and30_firstqt maxlag];
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Second"
                    xcorrKDE_Free_0and30_secondqt = [xcorrKDE_Free_0and30_secondqt maxlag];
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Third"
                    xcorrKDE_Free_0and30_thirdqt = [xcorrKDE_Free_0and30_thirdqt maxlag];
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Fourth"
                    xcorrKDE_Free_0and30_fourthqt = [xcorrKDE_Free_0and30_fourthqt maxlag];
            end
            
            
            % Get visualization for each subject
            if (DEBUG)
                subplot(3,length(SubjList),subj); polarplot([0:.1:360]/180*pi, vfEstimate); hold on; 
                polarplot([0:.1:360]/180*pi, vfEstimate2);
                title(sprintf('Subj %s',SubjList(subj)))
                subplot(3,length(SubjList),subj+length(SubjList)); plot(vfEstimate); hold on; 
                plot(vfEstimate2);
                subplot(3,length(SubjList),subj+length(SubjList)*2); plot(lags, x); line([maxlag*10 maxlag*10],get(gca,'ylim')); ...
                    %text(maxlag, max(x), sprintf('Max @ %.2f deg',maxlag), 'VerticalAlignment','bottom'); ...
                    set(gca,'ylim',[50 200]), yticklabels({})
            end
            
        end
    end
end

% Save the xcorr values in a table
tab.xcorrKDE_Free_0andneg30_firstqt(:,1) = xcorrKDE_Free_0andneg30_firstqt;
tab.xcorrKDE_Free_0andneg30_secondqt(:,1) = xcorrKDE_Free_0andneg30_secondqt;
tab.xcorrKDE_Free_0andneg30_thirdqt(:,1) = xcorrKDE_Free_0andneg30_thirdqt;
tab.xcorrKDE_Free_0andneg30_fourthqt(:,1) = xcorrKDE_Free_0andneg30_fourthqt;
tab.xcorrKDE_Free_0and30_firstqt(:,1) = xcorrKDE_Free_0and30_firstqt;
tab.xcorrKDE_Free_0and30_secondqt(:,1) = xcorrKDE_Free_0and30_secondqt;
tab.xcorrKDE_Free_0and30_thirdqt(:,1) = xcorrKDE_Free_0and30_thirdqt;
tab.xcorrKDE_Free_0and30_fourthqt(:,1) = xcorrKDE_Free_0and30_fourthqt;

%% Torsion
% For torsion stats
meanTIm30_Free = []; meanTIm0_Free = []; meanTImneg30_Free = [];
meanTIm30_Fix = []; meanTIm0_Fix = []; meanTImneg30_Fix = [];
for atask = 1:length(Tasks)
    for subj = 1:length(SubjList)
        meanTorsionIm30 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==30 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        meanTorsionIm0 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==0 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        meanTorsionImneg30 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==-30 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        
        switch (true)
            case Tasks(atask) == "FreeView"
                meanTIm30_Free = [meanTIm30_Free meanTorsionIm30];
                meanTIm0_Free = [meanTIm0_Free meanTorsionIm0];
                meanTImneg30_Free = [meanTImneg30_Free meanTorsionImneg30];
            case Tasks(atask) == "Fixation"
                meanTIm30_Fix = [meanTIm30_Fix meanTorsionIm30];
                meanTIm0_Fix = [meanTIm0_Fix meanTorsionIm0];
                meanTImneg30_Fix = [meanTImneg30_Fix meanTorsionImneg30];
        end
    end
end
tab.meanTIm30_Free = meanTIm30_Free';
tab.meanTIm0_Free = meanTIm0_Free';
tab.meanTImneg30_Free = meanTImneg30_Free';
tab.meanTIm30_Fix = meanTIm30_Fix';
tab.meanTIm0_Fix = meanTIm0_Fix';
tab.meanTImneg30_Fix = meanTImneg30_Fix';

avgTorsion = mean([tab.meanTIm30_Fix -tab.meanTImneg30_Fix],2)
[h,p,ci,stats] = ttest(avgTorsion)
mean(avgTorsion)



