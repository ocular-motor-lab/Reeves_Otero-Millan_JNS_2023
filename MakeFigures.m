% Figures and stats for Reeves & Otero-Millan 2023

% Load and initialize
SubjList = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20"];
projectFolder =  '/Users/stephaniereeves/UC Berkeley/OMlab - OM-lab-share/Projects/MicrosaccadesImageTilt/Manuscript/Reeves_Otero-Millan_JNS_2023'; %edit as needed

load(fullfile(projectFolder, 'saccade.mat'))
load(fullfile(projectFolder, 'trialTable.mat'))
load(fullfile(projectFolder, 'samplesDataTable.mat'))
load(fullfile(projectFolder, 'tab.mat'))
load(fullfile(projectFolder, 'bootstrappedData.mat'))


%% MANUSCRIPT FIGURE -- Main sequence for paper 

st = find(samplesDataTable.TrialNumber == 42,1,'first'); % trial 83 for 0186 is good 
en = find(samplesDataTable.TrialNumber == 42,1,'last');
st_fix = find(samplesDataTable.TrialNumber == 40,1,'first'); % trial 84 and 21 for 0186 is ok
en_fix = find(samplesDataTable.TrialNumber == 40,1,'last'); % trial 40 and 10 for 0137 

figure('Position', [415 705 1145 633])
subplot(2,3,1)
scatter(saccade.Amplitude(find(saccade.Task == "FreeView")), saccade.PeakSpeed(find(saccade.Task == "FreeView")),1,'MarkerEdgeColor',[.3 .3 .3]);
xlim([0 20])
xlabel('Amplitude (deg)');
ylabel('Peak Velocity (deg/s)');

subplot(2,3,2)
histAmp = histc(saccade.Amplitude(find(saccade.Task == "FreeView" & saccade.Amplitude < 20)), [0:1:20]);
b = bar([0:1:20],histAmp,'histc');
b.FaceColor = [.6 .6 .6];
xlim([0 20])
xlabel('Amplitude (deg)');
ylabel('N saccades');

subplot(2,3,3)
plot(samplesDataTable.Time(st:en)-samplesDataTable.Time(st), samplesDataTable.RightX(st:en)-mean(samplesDataTable.RightX(st:en),'omitnan'),'Color',[0.850 0.325 0.098]) 
hold
plot(samplesDataTable.Time(st:en)-samplesDataTable.Time(st), samplesDataTable.RightY(st:en)-mean(samplesDataTable.RightY(st:en),'omitnan'),'Color',[0.466 0.674 0.188])
ylim([-10 10]);
xlim([0 12]);
xlabel('Time (s)')
ylabel('Position (deg)')
legend('Horizontal','Vertical')

subplot(2,3,4)
scatter(saccade.Amplitude(find(saccade.Task == "Fixation" & saccade.Amplitude < 1.2)), saccade.PeakSpeed(find(saccade.Task == "Fixation" & saccade.Amplitude < 1.2)),1,'MarkerEdgeColor',[.3 .3 .3]);
ylim([0 150])
xlim([0 1.2])
xticks(0:0.2:1.2)
xticklabels({0, .2, .4, .6, .8, 1, 1.2})
xlabel('Amplitude (deg)');
ylabel('Peak Velocity (deg/s)');

subplot(2,3,5)
histAmp = histc(saccade.Amplitude(find(saccade.Task == "Fixation" & saccade.Amplitude < 1.2)), [0:.05:1.2]);
b = bar([0:.05:1.2],histAmp,'histc');
b.FaceColor = [.6 .6 .6];
xlim([0 1.2])
ylim([0 2000])
xticks(0:0.2:1.2)
xticklabels({0, .2, .4, .6, .8, 1, 1.2})
xlabel('Amplitude (deg)');
ylabel('N saccades');

subplot(2,3,6)
plot(samplesDataTable.Time(st_fix:en_fix)-samplesDataTable.Time(st_fix), samplesDataTable.RightX(st_fix:en_fix)-mean(samplesDataTable.RightX(st_fix:en_fix),'omitnan'),'Color',[0.850 0.325 0.098]) %nanmedfilt is a function that Jorge wrote for Arume and works by taking the sliding median of a certain number of samples, in this case 2000 samples
hold
plot(samplesDataTable.Time(st_fix:en_fix)-samplesDataTable.Time(st_fix), samplesDataTable.RightY(st_fix:en_fix)-mean(samplesDataTable.RightY(st_fix:en_fix),'omitnan'),'Color',[0.466 0.674 0.188]) %nanmedfilt is a function that Jorge wrote for Arume and works by taking the sliding median of a certain number of samples, in this case 2000 samples
ylim([-10 10]);
xlim([0 12]);
xlabel('Time (s)')
ylabel('Position (deg)')
legend('Horizontal','Vertical')



%% MANSUCRIPT FIGURE: Free viewing saccade directions 
binedges = [0:.1:360]/180*pi;
test1 = [];
test2 = [];
test3 = [];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate1] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate1];
    
    h2 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
end

% Prepare the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], length(SubjList)-1); %13 bc 14-1 is the dof 
se1 = std(tab.xcorrKDE_Free_0andneg30,'omitnan') / sqrt(length(tab.xcorrKDE_Free_0andneg30));
ci1 = mean(tab.xcorrKDE_Free_0andneg30,'omitnan') + crit*se1;
se2 = std(tab.xcorrKDE_Free_0and30,'omitnan') / sqrt(length(tab.xcorrKDE_Free_0and30));
ci2 = mean(tab.xcorrKDE_Free_0and30,'omitnan') + crit*se2;

figure('Position',1.0e+03 * [0.0100  0.0610  1.1243  0.9867]);
sp1 = subplot(2,10,1:5)
sp1.Position = sp1.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold
polarplot(binedges,mean(test2),'LineWidth',2, 'Color', get_color('bondiBlue7')) 
leg = legend('0 deg Tilt','-30 deg Tilt')
leg.Location = 'south';

sp2 = subplot(2,10,6:10)
sp2.Position = sp2.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold
polarplot(binedges,mean(test1),':','LineWidth',2, 'Color', get_color('bondiBlue7')) 
leg = legend('0 deg Tilt','30 deg Tilt')
leg.Location = 'south';

sp3 = subplot(2,10,11:15)
bar([1.2:length(SubjList)+.2],bootstrapped_free_right.Avg,'FaceColor',[1 1 1],'EdgeColor',get_color('bondiBlue7'), 'LineWidth',2,'LineStyle',':'); hold on;
bar([1:length(SubjList)],bootstrapped_free_left.Avg,'FaceColor',get_color('bondiBlue7'),'EdgeColor',get_color('bondiBlue7'), 'LineWidth',.5); hold on;
for subj = 1:length(SubjList)
    line([subj+.2 subj+.2],[bootstrapped_free_right.CI_low(subj) bootstrapped_free_right.CI_high(subj) ],'Color',get_color('darkGrey'),'LineWidth', .75)
    line([subj subj],[bootstrapped_free_left.CI_low(subj) bootstrapped_free_left.CI_high(subj) ],'Color',get_color('darkGrey'),'LineWidth', .75)
end
xticks(1:length(SubjList))
xticklabels(bootstrapped_free_right.SubjNum'); %set(gca,'XTickLabelRotation',80);
xlabel('Individual Subjects')
ylabel('Direction Distribution Displacement (deg)')
ylim([-37 37])
set(gca,'FontSize',12);

sp4 = subplot(2,10,16:17)
sp4.Position = sp4.Position - [0.02 0 0 0]
p1y = mean(tab.xcorrKDE_Free_0andneg30);
p2y = mean(tab.xcorrKDE_Free_0and30);
p1 = bar(1, p1y,'FaceColor',get_color('bondiBlue7'),'EdgeColor',get_color('bondiBlue7'), 'LineWidth',2) % edge color used to be .3 .3 .3
hold on
p2 = bar(2, p2y,'FaceColor',[1 1 1],'EdgeColor',get_color('bondiBlue7'), 'LineWidth',2,'LineStyle',':') % edge color used to be .3 .3 .3
for i = 1:length(tab.xcorrKDE_Free_0andneg30)
    scatter(1, tab.xcorrKDE_Free_0andneg30(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(tab.xcorrKDE_Free_0and30)
    scatter(2, tab.xcorrKDE_Free_0and30(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .75)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .75)
xticks(1:2)
xticklabels({'-30','30'})
yticklabels({})
xlabel('Image Tilt')
%ylabel('{\Delta} (deg)','FontSize',13)
ylim([-37 37])
set(gca,'FontSize',12);

sp5 = subplot(2,10,18:20)
%sp4.Position = sp4.Position + [0 0 0.015 0];
scaled_freeview = (mean([tab.xcorrKDE_Free_0and30 -tab.xcorrKDE_Free_0andneg30],2))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 
bp = boxplot(scaled_freeview, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',1.5)
hold
for i = 1:length(scaled_freeview)
    scatter(1, scaled_freeview(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
ylim([-0.1 1.1])
xticklabels({' '})
ylabel('Reference Frame Index','FontSize',13)
txt = {'\leftarrow Egocentric'};
text(.5, 0, txt, 'FontSize',12)
txt2 = {'\leftarrow Image'};
text(.5, 1, txt2, 'FontSize',12)
set(gca,'FontSize',12);

% Stats for manuscript
mean(tab.xcorrKDE_Free_0andneg30)
mean(tab.xcorrKDE_Free_0and30)
mean(scaled_freeview)
[h,p,ci,stats] = ttest(tab.xcorrKDE_Free_0andneg30)
[h,p,ci,stats] = ttest(tab.xcorrKDE_Free_0and30)
[h,p,ci,stats] = ttest(scaled_freeview)
[h,p,ci,stats] = ttest(1-scaled_freeview)

%% MANUSCRIPT FIGURE: Fixation saccade directions 
binedges = [0:.1:360]/180*pi;
test1 = [];
test2 = [];
test3 = [];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "Fixation" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate1] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate1];
    
    h2 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "Fixation" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Subject == SubjList(i) & saccade.Task == "Fixation" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];
end

% Prepare the confidence intervals
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], length(SubjList)-1); %13 bc 14-1 is the dof
se1 = std(tab.xcorrKDE_Fix_0andneg30) / sqrt(length(tab.xcorrKDE_Fix_0andneg30));
ci1 = mean(tab.xcorrKDE_Fix_0andneg30) + crit*se1;
se2 = std(tab.xcorrKDE_Fix_0and30) / sqrt(length(tab.xcorrKDE_Fix_0and30));
ci2 = mean(tab.xcorrKDE_Fix_0and30) + crit*se2;

figure('Position',1.0e+03 * [0.0100  0.0610  1.1243  0.9867]);
sp1 = subplot(2,10,1:5)
sp1.Position = sp1.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k') 
hold
polarplot(binedges,mean(test2),'LineWidth',2,'Color', get_color('salmon3')) 
leg = legend('0 deg Tilt','-30 deg Tilt')
leg.Location = 'south';

sp2 = subplot(2,10,6:10)
sp2.Position = sp2.Position - [0 0.06 0 0];
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', 'k') 
hold
polarplot(binedges,mean(test1),'LineWidth',2, 'LineStyle',':','Color', get_color('salmon3')) 
leg = legend('0 deg Tilt','30 deg Tilt')
leg.Location = 'south';

sp3 = subplot(2,10,11:15)
bar([1.2:length(SubjList)+.2],bootstrapped_fix_right.Avg,'FaceColor',[1 1 1],'EdgeColor',get_color('salmon3'), 'LineWidth',2,'LineStyle',':'); hold on;
bar([1:length(SubjList)],bootstrapped_fix_left.Avg,'FaceColor',get_color('salmon3'),'EdgeColor',get_color('salmon3'), 'LineWidth',.5); hold on;
for subj = 1:length(SubjList)
    line([subj+.2 subj+.2],[bootstrapped_fix_right.CI_low(subj) bootstrapped_fix_right.CI_high(subj) ],'Color','k','LineWidth', .75)
    line([subj subj],[bootstrapped_fix_left.CI_low(subj) bootstrapped_fix_left.CI_high(subj) ],'Color','k','LineWidth', .75)
end
xticks(1:length(SubjList))
xticklabels(bootstrapped_fix_right.SubjNum'); %set(gca,'XTickLabelRotation',80);
xlabel('Individual Subjects')
ylabel('Direction Distribution Displacement (deg)')
ylim([-35 35])
set(gca,'FontSize',12);

sp4 = subplot(2,10,16:17)
sp4.Position = sp4.Position - [0.02 0 0 0]
p1y = mean(tab.xcorrKDE_Fix_0andneg30);
p2y = mean(tab.xcorrKDE_Fix_0and30);
p1 = bar(1, p1y,'FaceColor',get_color('salmon3'),'EdgeColor',get_color('salmon3'), 'LineWidth',2) % edge color used to be .3 .3 .3
hold on
p2 = bar(2, p2y,'FaceColor',[1 1 1],'EdgeColor',get_color('salmon3'), 'LineWidth',2,'LineStyle',':') % edge color used to be .3 .3 .3
for i = 1:length(tab.xcorrKDE_Fix_0andneg30)
    scatter(1, tab.xcorrKDE_Fix_0andneg30(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
for i = 1:length(tab.xcorrKDE_Fix_0and30)
    scatter(2, tab.xcorrKDE_Fix_0and30(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1)
end
line([1 1],[ci1(1) ci1(2)],'Color','k','LineWidth', .75)
line([2 2],[ci2(1) ci2(2)],'Color','k','LineWidth', .75)
xticks(1:2)
xticklabels({'-30','30'})
yticklabels({})
xlabel('Image Tilt')
%ylabel('{\Delta} (deg)','FontSize',13)
ylim([-35 35])
set(gca,'FontSize',12);


sp5 = subplot(2,10,18:20)
%sp5.Position = sp4.Position + [0 0 0.015 0];
scaled_fixation = (mean([tab.xcorrKDE_Fix_0and30 -tab.xcorrKDE_Fix_0andneg30],2))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 
bp = boxplot(scaled_fixation, 'Colors','k','Symbol','+k') %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',1.5)
hold
for i = 1:length(scaled_fixation)
    scatter(1, scaled_fixation(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
ylim([-0.2 1.1])
xticklabels({' '})
ylabel('Reference Frame Index','FontSize',13)
txt = {'\leftarrow Ego'};
text(.5, 0, txt, 'FontSize',12)
txt2 = {'\leftarrow Image'};
text(.5, 1, txt2, 'FontSize',12)
set(gca,'FontSize',12);

% Stats for manuscript
mean(tab.xcorrKDE_Fix_0andneg30)
mean(tab.xcorrKDE_Fix_0and30)
mean(scaled_fixation)
[h,p,ci,stats] = ttest(tab.xcorrKDE_Fix_0andneg30)
[h,p,ci,stats] = ttest(tab.xcorrKDE_Fix_0and30)
[h,p,ci,stats] = ttest(scaled_fixation)
[h,p,ci,stats] = ttest(1-scaled_fixation)



%% MANUSCRIPT FIGURE Amplitude effect during free viewing
binedges = [0:.1:360]/180*pi;
test1 = []; test2 = []; test3 = []; test4 = []; test5 = []; test6 = [];
test7 = []; test8 = []; test9 = []; test10 = []; test11 = []; test12 = [];

scaled_fixation = (mean([tab.xcorrKDE_Fix_0and30 -tab.xcorrKDE_Fix_0andneg30],2))./30; 
scaled_0to25 = (mean([tab.xcorrKDE_Free_0and30_firstqt -tab.xcorrKDE_Free_0andneg30_firstqt],2,'omitnan'))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 
scaled_25to50 = (mean([tab.xcorrKDE_Free_0and30_secondqt -tab.xcorrKDE_Free_0andneg30_secondqt],2,'omitnan'))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 
scaled_50to75 = (mean([tab.xcorrKDE_Free_0and30_thirdqt -tab.xcorrKDE_Free_0andneg30_thirdqt],2,'omitnan'))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 
scaled_75to100 = (mean([tab.xcorrKDE_Free_0and30_fourthqt -tab.xcorrKDE_Free_0andneg30_fourthqt],2,'omitnan'))./30; % divide by 30 bc 30/30 would be 1 which is image-centered 

% colors = [([150, 222, 209]./255); ([138, 154, 91]./255); ([147, 197, 114]./255); ([180, 196, 36]./255);  ([201, 204, 63]./255)  ];
%colors = [ get_color('salmon3'); get_color('bondiBlue4'); get_color('bondiBlue7'); get_color('bondiBlue10'); get_color('bondiBlue13') ];
colors = [ get_color('bondiBlue4'); get_color('bondiBlue7'); get_color('bondiBlue10'); get_color('bondiBlue13') ];

% Get the vals for each subject and stitch them together
for i = 1:length(SubjList)
    h1 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate1] = circ_ksdensity(h1,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test1 = [test1; vfEstimate1];
    
    h2 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate2] = circ_ksdensity(h2,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test2 = [test2; vfEstimate2];
    
    h3 = saccade.Direction(find(saccade.Amplitude < Q_25 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate3] = circ_ksdensity(h3,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test3 = [test3; vfEstimate3];

    h4 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate4] = circ_ksdensity(h4,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test4 = [test4; vfEstimate4];
    
    h5 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate5] = circ_ksdensity(h5,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test5 = [test5; vfEstimate5];
    
    h6 = saccade.Direction(find(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate6] = circ_ksdensity(h6,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test6 = [test6; vfEstimate6];
    
    h7 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate7] = circ_ksdensity(h7,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test7 = [test7; vfEstimate7];
    
    h8 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate8] = circ_ksdensity(h8,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test8 = [test8; vfEstimate8];
    
    h9 = saccade.Direction(find(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate9] = circ_ksdensity(h9,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test9 = [test9; vfEstimate9];
    
    h10 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 30 & saccade.TrialNumber > 0));
    [vfEstimate10] = circ_ksdensity(h10,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test10 = [test10; vfEstimate10];
    
    h11 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == 0 & saccade.TrialNumber > 0));
    [vfEstimate11] = circ_ksdensity(h11,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test11 = [test11; vfEstimate11];
    
    h12 = saccade.Direction(find(saccade.Amplitude >= Q_75 & saccade.Subject == SubjList(i) & saccade.Task == "FreeView" & saccade.ImTilt == -30 & saccade.TrialNumber > 0));
    [vfEstimate12] = circ_ksdensity(h12,0:.1/180*pi:2*pi,[-pi pi],.1) 
    test12 = [test12; vfEstimate12];
end

figure('Position', 1.0e+03 * [0.1950    0.1677    1.3367    0.6213])
sp1 = subplot(2,4,1)
polarplot(binedges,mean(test2),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold on
polarplot(binedges,mean(test1),':','LineWidth',2, 'Color', get_color('bondiBlue4'))
polarplot(binedges,mean(test3),'LineWidth',2, 'Color', get_color('bondiBlue4'))
rlim([0 0.5]);
leg = legend('0 deg Tilt','30 deg Tilt','-30 deg Tilt')
leg.Location = 'south';
title("< 1.1 degree",'FontSize',13)

sp2 = subplot(2,4,2)
polarplot(binedges,mean(test5),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold on
polarplot(binedges,mean(test4),':','LineWidth',2, 'Color', get_color('bondiBlue7'))
polarplot(binedges,mean(test6),'LineWidth',2, 'Color', get_color('bondiBlue7'))
rlim([0 0.5]);
leg = legend('0 deg Tilt','30 deg Tilt','-30 deg Tilt')
leg.Location = 'south';
title("1.1 to 2.3 degrees",'FontSize',13)

sp3 = subplot(2,4,5)
polarplot(binedges,mean(test8),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold on
polarplot(binedges,mean(test7),':','LineWidth',2, 'Color', get_color('bondiBlue10'))
polarplot(binedges,mean(test9),'LineWidth',2, 'Color', get_color('bondiBlue10'))
rlim([0 0.5]);
leg = legend('0 deg Tilt','30 deg Tilt','-30 deg Tilt')
leg.Location = 'south';
title("2.3 to 4.7 degrees",'FontSize',13)

sp4 = subplot(2,4,6)
polarplot(binedges,mean(test11),'LineWidth',2, 'Color', get_color('darkGrey')) 
hold on
polarplot(binedges,mean(test10),':','LineWidth',2, 'Color', get_color('bondiBlue13'))
polarplot(binedges,mean(test12),'LineWidth',2, 'Color', get_color('bondiBlue13'))
rlim([0 0.5]);
leg = legend('0 deg Tilt','30 deg Tilt','-30 deg Tilt')
leg.Location = 'south';
title("> 4.7 degrees",'FontSize',13)

sp4 = subplot(2,4,[3 4])
bp = boxplot([scaled_0to25,scaled_25to50,scaled_50to75,scaled_75to100], 'Colors',colors,'Symbol','+k','Whisker',1) %the +k means use the symbol + for the outlier and use the color k for the outlier
set(bp,'LineWidth',2)
hold on
for i = 1:length(scaled_0to25)
    scatter(1, scaled_0to25(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
for i = 1:length(scaled_25to50)
    scatter(2, scaled_25to50(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
for i = 1:length(scaled_50to75)
    scatter(3, scaled_50to75(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
for i = 1:length(scaled_75to100)
    scatter(4, scaled_75to100(i),20,'MarkerEdgeColor',[.6 .6 .6], 'LineWidth', 1, 'jitter','on','jitterAmount',.01)
end
yline(0,'--')
yticks(-1:.2:1.2)
yticklabels({-1:.2:1.2})
ylabel('Reference Frame Index','FontSize',13)
set(gca,'FontSize',13,'XTick',[],'XTickLabel',[])


sp5 = subplot(2,4,[7 8])
sp5.Position = [0.57 0.125 0.335 0.14]
histAmp1 = histc(saccade.Amplitude(find(saccade.Task == "FreeView" & saccade.Amplitude < Q_25)), [0:.1:7]);
histAmp2 = histc(saccade.Amplitude(find(saccade.Task == "FreeView" & saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50)), [0:.1:7]);
histAmp3 = histc(saccade.Amplitude(find(saccade.Task == "FreeView" & saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75)), [0:.1:7]);
histAmp4 = histc(saccade.Amplitude(find(saccade.Task == "FreeView" & saccade.Amplitude >= Q_75)), [0:.1:7]);
b = bar([0:.1:7],[histAmp1 histAmp2 histAmp3 histAmp4],'stacked');
b(1).FaceColor = get_color('bondiBlue4'); b(1).LineStyle = 'none';
b(2).FaceColor = get_color('bondiBlue7'); b(2).LineStyle = 'none';
b(3).FaceColor = get_color('bondiBlue10'); b(3).LineStyle = 'none';
b(4).FaceColor = get_color('bondiBlue13'); b(4).LineStyle = 'none';
xlim([0 7])
ylim([0 1750])
xticks([1.1,2.3,4.7])
xlabel('Amplitude (deg)')
set(gca,'FontSize',13)

sp4.Position = [0.5689 0.3100 0.3357 0.6150];


mean(scaled_0to25)
mean(scaled_25to50)
mean(scaled_50to75)
mean(scaled_75to100)

% Get the confidence intervals for the paper
pLo = 0.05/2; %0.05 bc 95% conf intervals
pUp = 1 - 0.05/2;
crit = tinv([pLo pUp], 13); %13 bc 14-1 is the dof 

se1 = std(scaled_0to25) / sqrt(length(scaled_0to25));
lowerCI = mean(scaled_0to25) + crit(1)*se1
upperCI = mean(scaled_0to25) + crit(2)*se1

se2 = std(scaled_25to50) / sqrt(length(scaled_25to50));
lowerCI = mean(scaled_25to50) + crit(1)*se1
upperCI = mean(scaled_25to50) + crit(2)*se1

se3 = std(scaled_50to75) / sqrt(length(scaled_50to75));
lowerCI = mean(scaled_50to75) + crit(1)*se1
upperCI = mean(scaled_50to75) + crit(2)*se1

se4 = std(scaled_75to100) / sqrt(length(scaled_75to100));
lowerCI = mean(scaled_75to100) + crit(1)*se1
upperCI = mean(scaled_75to100) + crit(2)*se1


% Stats for manuscript
% Linear Mixed Model
t = table();
t.Subject = [ SubjList'; SubjList'; SubjList'; SubjList'; ];
t.RFIs = [scaled_0to25; scaled_25to50; scaled_50to75; scaled_75to100];
t.Quartile = [ repmat(1,[length(SubjList) 1]); repmat(2,[length(SubjList) 1]); repmat(3,[length(SubjList) 1]); repmat(4,[length(SubjList) 1])];
for subj = 1:length(SubjList)
    t.MeanAmplitude(t.Subject == SubjList(subj) & t.Quartile == 1) = mean(saccade.Amplitude(saccade.Amplitude < Q_25 & saccade.Subject == SubjList(subj)),'omitnan')
    t.MeanAmplitude(t.Subject == SubjList(subj) & t.Quartile == 2) = mean(saccade.Amplitude(saccade.Amplitude >= Q_25 & saccade.Amplitude < Q_50 & saccade.Subject == SubjList(subj)),'omitnan')
    t.MeanAmplitude(t.Subject == SubjList(subj) & t.Quartile == 3) = mean(saccade.Amplitude(saccade.Amplitude >= Q_50 & saccade.Amplitude < Q_75 & saccade.Subject == SubjList(subj)),'omitnan')
    t.MeanAmplitude(t.Subject == SubjList(subj) & t.Quartile == 4) = mean(saccade.Amplitude(saccade.Amplitude >= Q_75 & saccade.Subject == SubjList(subj)),'omitnan')
end
lme = fitlme(t,'RFIs ~ MeanAmplitude + (MeanAmplitude|Subject)') % this is the correct one



