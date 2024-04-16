% Reeves & Otero-Millan 2023

% This script creates a subject data table with relevant data and saves it as
% "tab" and also saves a couple of bootstrapped analyses and saves them as
% "bootstrappedData"

%% Initialize
SubjList = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"];
projectFolder =  '/Users/stephaniereeves/UC Berkeley/OMlab - OM-lab-share/Projects/MicrosaccadesImageTilt/Manuscript/Reeves_Otero-Millan_JNS_2023'; %edit as needed

% Load the saccade and trialTable
load(fullfile(projectFolder,'saccade.mat'))
load(fullfile(projectFolder,'trialTable.mat'))


%% Get polar histograms (using KDE) and do circular cross-correlations 
% If DEBUG = 1, then plots will print out 
DEBUG = 0;

% Initializing things
Tasks = ["Fixation","FreeView"];
LRTilts = [-30,30];
tab = table();
tab.Subject = transpose(SubjList);

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
                    tab.xcorrKDE_Fix_0andneg30(subj) = maxlag;
                case LRTilts(imtilt) == 30 && Tasks(atask) == "Fixation"
                    tab.xcorrKDE_Fix_0and30(subj) = maxlag;
                case LRTilts(imtilt) == -30 && Tasks(atask) == "FreeView"
                    tab.xcorrKDE_Free_0andneg30(subj) = maxlag;
                case LRTilts(imtilt) == 30 && Tasks(atask) == "FreeView"
                    tab.xcorrKDE_Free_0and30(subj) = maxlag;
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
                    tab.xcorrKDE_Free_0andneg30_firstqt(subj) = maxlag;
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Second"
                    tab.xcorrKDE_Free_0andneg30_secondqt(subj) = maxlag;
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Third"
                    tab.xcorrKDE_Free_0andneg30_thirdqt(subj) = maxlag;
                case LRTilts(imtilt) == "-30" && Quartiles(aquartile) == "Fourth"
                    tab.xcorrKDE_Free_0andneg30_fourthqt(subj) = maxlag;
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "First"
                    tab.xcorrKDE_Free_0and30_firstqt(subj) = maxlag;
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Second"
                    tab.xcorrKDE_Free_0and30_secondqt(subj) = maxlag;
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Third"
                    tab.xcorrKDE_Free_0and30_thirdqt(subj) = maxlag;
                case LRTilts(imtilt) == "30" && Quartiles(aquartile) == "Fourth"
                    tab.xcorrKDE_Free_0and30_fourthqt(subj) = maxlag;
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

%% Torsion
% For torsion stats

for atask = 1:length(Tasks)
    for subj = 1:length(SubjList)
        meanTorsionIm30 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==30 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        meanTorsionIm0 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==0 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        meanTorsionImneg30 = mean(trialTable.medianT_AccountingForBaselineMeanT(trialTable.ImTilt==-30 & trialTable.Task==Tasks(atask) & trialTable.Subject == SubjList(subj)),'omitnan');
        
        switch (true)
            case Tasks(atask) == "FreeView"
                tab.meanTIm30_Free(subj) = [ meanTorsionIm30];
                tab.meanTIm0_Free(subj) = [ meanTorsionIm0];
                tab.meanTImneg30_Free(subj) = [ meanTorsionImneg30];
            case Tasks(atask) == "Fixation"
                tab.meanTIm30_Fix(subj) = [ meanTorsionIm30];
                tab.meanTIm0_Fix(subj) = [ meanTorsionIm0];
                tab.meanTImneg30_Fix(subj) = [ meanTorsionImneg30];
        end
    end
end

avgTorsion = mean([tab.meanTIm30_Fix -tab.meanTImneg30_Fix],2)
[h,p,ci,stats] = ttest(avgTorsion)
mean(avgTorsion)

%% Save
save(fullfile(projectFolder,"tab.mat"),"tab")
save(fullfile(projectFolder,"bootstrappedData.mat"),"bootstrapped_fix_left","bootstrapped_fix_right","bootstrapped_free_left","bootstrapped_free_right")


