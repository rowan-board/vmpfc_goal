%% Define paths 

homepath = 'path...';
addpath(genpath(homepath)); 

neural_dat = [homepath '/analyse_timeCourses/'];

participants = {'02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', '23'}; 

ROIs = {'VS', 'dmpfc', 'dmpfc_ctrl', 'lFP'};
ROIs = {'amy_rwd_patient', 'amy_rwd_ctrl', 'amy_m_patient', 'amy_m_ctrl'};

%% Load data
% Load behavioural data
subjList = str2double(participants); 

% Load neural data
timefolder = fullfile(neural_dat);

% Load timings and regressor table
load('mri_timings.mat');
load('regressor_table.mat');


%% Upsample time course by event
upsampledData = upsample_tc(timefolder,mri_timings,ROIs,subjList); 
choiceOnset = upsampledData.choiceOnset;
%% Run the regression
regressors= {'trials','reward','chosen_ev','reward_diff','chosen_diff','motivation'};

% Run model
nRegressors = numel(regressors); 
treat_regressors=0; %  0=unnormalised values, 1=normalise, 2 = demean (but
% don't divide by variance), 3=normalise but ignoring all zeros 
analysis_regions = 1:numel(ROIs); 


nModels = 1;
models = cell(nModels,1);
intercept = 0;

[models{1}.betas,models{1}.SE, models{1}.tStat] = run_tc_GLM_byRun(choiceOnset, subjList, data_table,  regressors,intercept, treat_regressors, analysis_regions);
models{1}.name = 'all trials';

%% plot regression of patient vs controls
nModels = numel(models); 
ntimep = size(choiceOnset.Region(1).epoched_timecourses{1}.all_runs,2);  % number of timepoints in each epoch
modelColours = {'r','b','g','k'};

regressorNames = {'Trials','Objective Reward','Objective EV','Devaluation of Reward','Devaluation of EV','Motivation'};

for r = 1:numel(regressorNames)
    regressorNames{r} = strrep(regressorNames{r},'_',' ');
end

regionNames = [];
for reg = 1:numel(analysis_regions)
    regionNames = [regionNames, {strrep(choiceOnset.Region(reg).name,'_',' ')}];
end

changeColour = 1;
for reg = 1:numel(analysis_regions)
    
    figure; 
    sgtitle(regionNames(reg)); 
    
    count = 0;
    for b = 1:nRegressors
        
        count = count+1;
        subplot(4,4,count);
        title(regressorNames{b});  

        for m = 1:nModels

           beta_weights = models{1}.betas; 
           standarderror = std(beta_weights(1:21,reg,:,b))./sqrt(length(subjList)); 
           y = mean(beta_weights(1:21,reg,:,b)); 
           x = [choiceOnset.window/ntimep:choiceOnset.window/ntimep:choiceOnset.window]; 
           xlabel('time (s)');
           ylabel('beta coefficient');
           
           shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{1}});

           hold on 
           standarderror = beta_weights(22,reg,:,b)./sqrt(length(subjList)); 
           y = beta_weights(22,reg,:,b);
           shadedErrorBar(x,y,standarderror,'lineProps',{modelColours{2}});
           legend({'ctrl','pat'},'location','northeastoutside','Orientation','vertical','AutoUpdate','off')
           hold off

           xticks([0:2:choiceOnset.window])
           xlim([0, 12]); 
           hold on; 
        end

        % set title as axis labels
        title(regressorNames{b});
        yline(0,'k'); 
        xlabel('time (s)');
        ylabel('mean signal change (a.u)');

    end
end

%% extract beta weights and plot averages for patient vs control
regressorNames_nT = {'Objective Reward','Objective EV','dmPFC Dev. EV', 'dmPFC_ctrl', 'lFP Dev. EV'};
ROIs_nT = {'Ventral Striatum', 'dmPFC Dev. EV', 'dmPFC_ctrl', 'lFP Dev. EV'};
ROIs_nT = {'amygdala rwd patient', 'amygdala rwd ctrl', 'amygdala m patient', 'amydala m ctrl'};
subj_betas = zeros(length(subjList),nRegressors,length(regionNames));

beta_weights = models{1}.betas; 

for i=1:length(subjList)
    for j=1:nRegressors
        for k=1:length(regionNames)
            subj_betas(i,j,k)=mean(beta_weights(i,k,:,j));
        end
    end
end

% bar plot of patient vs controls

for k=1:length(regionNames)
    ctrl=mean(subj_betas(1:21,2:6,k));
    ctrl_e=std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)));
    pat=[0 0 0 0 0];
    pat_e=[0 0 0 0 0];
    for i=2:nRegressors
        pat(i-1)=mean(subj_betas(22,i,k));
        pat_e(i-1)=mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2;
    end

    figure;
    b=bar([1,2,3,4,5],[ctrl;pat]);
    c=10;
    hold on
    swarmchart([0.6,1.6,2.6,3.6,4.6],[subj_betas(1:21,2:6,k)],'MarkerEdgeColor',[0, 0.4470, 0.74100],'MarkerFaceColor',[0, 0.4470, 0.7410],SizeData=100);
    swarmchart([0.6,1.6,2.6,3.6,4.6],pat,'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],SizeData=100);
    er=errorbar([0.86,1.86,2.86,3.86,4.86],mean(subj_betas(1:21,2:6,k)),(std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)))));
    for i=2:nRegressors
        er_c=errorbar(i-0.86,mean(subj_betas(22,i,k)),(mean(models{1}.betas(22,k,:,i))/mean(models{1}.tStat(22,k,:,i))/2), Color='black');
    end
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    hold off

    title([ROIs_nT{k}],'FontSize',14); 
    yline(0,'k'); 
    xlabel('regressor','FontSize',20);
    ylabel('average beta coefficient','FontSize',20);
    xticklabels(regressorNames_nT)
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14; 

    legend(['Control';'Patient'],Location="northeastoutside")

    figure;
    hold on
    ctrl_dots = swarmchart([1,2,3,4,5], [subj_betas(1:21,2:6,k)] ,'MarkerEdgeColor',[0, 0.4470, 0.74100],'MarkerFaceColor',[0, 0.4470, 0.7410],SizeData=100);
    er=errorbar([1.1,2.1,3.1,4.1,5.1],mean(subj_betas(1:21,2:6,k)),(std(subj_betas(1:21,2:6,k))/sqrt(length(subj_betas(1:21,2:6,k)))));
    ctrl_point=plot([1.1,2.1,3.1,4.1,5.1], mean(subj_betas(1:21,2:6,k)), marker='o', MarkerEdgeColor=[0, 0.4470, 0.74100], MarkerFaceColor=[0, 0.4470, 0.7410], MarkerSize=10);

    er.LineStyle = 'none';
    ctrl_point.LineStyle = 'none';
    er.LineWidth = 2;

    pat_dots = swarmchart([1.25, 2.25, 3.25, 4.25, 5.25],[pat(1:5)],'MarkerEdgeColor',[0.8500, 0.3250, 0.0980],'MarkerFaceColor',[0.8500, 0.3250, 0.0980],SizeData=100);
    for i=1:5
        er_c=errorbar(i+0.25,mean(subj_betas(22,i+1,k)),(mean(models{1}.betas(22,k,:,i+1))/mean(models{1}.tStat(22,k,:,i+1))), Color=[0.8500, 0.3250, 0.0980]);
        %er_c=errorbar(i-2.75,mean(subj_betas(22,i,k)),(mean(models{1}.SE(22,k,:,i))), Color=[0.8500, 0.3250, 0.0980]);
        er_c.LineWidth = 2;
    end
                           
    hold off

    title([ROIs_nT{k}],'FontSize',14); 
    yline(0,'k'); 
    xlabel('regressor','FontSize',20);
    ylabel('average beta coefficient','FontSize',20);
    xticks([1 2 3 4 5]);
    xticklabels(regressorNames_nT(1:5));
    ax = gca;
    ax.XAxis.FontSize = 14; 
    ax.YAxis.FontSize = 14; 

    legend([ctrl_dots(1), pat_dots], {'Control', 'Patient'}, Location="northeastoutside")
end



% patient vs control csv roi
% now lets save the averages for each individual in csv files

for i=1:nRegressors
    for j=1:length(regionNames)
        tmp=subj_betas(:,i,j);
        writematrix(tmp, strjoin([neural_dat, string(ROIs(j)), '/', string(regressors(i)), '_', string(ROIs(j)), '.csv'],''));
    end
end

% lets also save the individual variances for permutation testing 

for i=1:nRegressors
    for j=1:length(regionNames)
        %tmp=mean(models{1}.betas(:,j,:,i),3)./mean(models{1}.tStat(:,j,:,i),3); % the '3' means you calc mean over the 3rd dimension of the matrix
        tmp=mean(models{1}.SE(:,j,:,i),3);
        writematrix(tmp, strjoin([neural_dat, string(ROIs(j)), '/', string(regressors(i)), '_', string(ROIs(j)), '_se.csv'],''));
    end
end


