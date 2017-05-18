% Identify channels where there is a significant MMN effect, using
% non-parametric cluster stats in Fieldtrip
%
% REG and RAND will be pooled (yielding a Dev and noDev condition) and the
% main effect of deviance assessed as (dev_REG + dev_RANd)- (REG + RAND).
% The contrast is performed in time x channels space. Clusters showing a
% deviance effect will then form the ROI time-channel window for the
% contrast (dev_REG-REG) - (dev_RAND - RAND)

close all; clearvars
sub_list        = [1:13 15:21]; % subjects to include
load neighboursBiosemi128.mat;
condlist        = {'REG10', 'RAND10','dev_REG','dev_RAND'}; % conditions in dataset
statpairs       = [1 2; 3 4]; % pairs of conditions to pool for stats
dir_ft = 'FTv3/';
path_in = [dir_ft 'average_dev/'];
file_in = 'F30ABCIOP2Hz' ;
fig_path = 'Plots/FTv3/Dev/'; mkdir(fig_path)
saveplots = 1;
%% load data
si=0;
[subplotRC,nsub]=numSubplots(numel(sub_list));
figure(77);clf % for topo of deviant response

try load([dir_ft 'dev_ALLSUB_' file_in '.mat']);
catch
    allData = [];
    for s = sub_list
        si = si+1;
        
        %% load each condition and put in structure for all subjects
        load([path_in file_in '_s' num2str(s) '.mat']);
        REG{si} = aveData.REG10;
        RAND{si} = aveData.RAND10;
        dev_REG{si} = aveData.dev_REG;
        dev_RAND{si} = aveData.dev_RAND;
        
        % find the differrence waveforms & save into a fieldtrip structure
        diff_REG{si} = REG{si};
        diff_REG{si}.avg = dev_REG{si}.avg - REG{si}.avg;
        
        % find difference waveform for REG and RAND - 'deviant response'
        cfg = [];
        cfg.parameter = 'avg';
        cfg.operation = 'subtract';
        diff_REG{si} = ft_math(cfg,dev_REG{si},REG{si});
        diff_RAND{si} = ft_math(cfg,dev_RAND{si},RAND{si});
        
        %% pool across reg and rand by making average using ft_timelock
        cfg = [];
        RR{si} = ft_appenddata(cfg,REG{si},RAND{si});
        dev_RR{si} = ft_appenddata(cfg,dev_REG{si}, dev_RAND{si});
        cfg = [];    cfg.keeptrials = 'yes';
        RR{si} = ft_timelockanalysis(cfg, RR{si});
        dev_RR{si} = ft_timelockanalysis(cfg, dev_RR{si});
        
        % find the difference waveform
        cfg = [];
        cfg.parameter = 'avg';
        cfg.operation = 'subtract';
        diff_RR{si} = ft_math(cfg,dev_RR{si},RR{si});
        
        % pack the AVG field into a single variable for using in non-fieldtrip
        % plotting functions later
        allData(si,1,:,:) = REG{si}.avg;
        allData(si,2,:,:) = RAND{si}.avg;
        allData(si,3,:,:) = dev_REG{si}.avg;
        allData(si,4,:,:) = dev_RAND{si}.avg;
        
    end
    save([dir_ft 'dev_ALLSUB_' file_in '.mat'] ,'allData','REG','RAND','dev_REG',...
        'dev_RAND','dev_RR','RR','diff_RR','diff_REG','diff_RAND');
end

%% plot topo of the deviant response
si=0;    fh=figure(77); clf
mtit(['All subjects deviant response 0.08 - 0.12 s ' file_in] )
for s = sub_list
    si=si+1;
    % plot each subject error response topo
    subplot(subplotRC(1),subplotRC(2),si);
    cfg=[];
    cfg.latency = [0.09 0.12];
    cfg.layout = 'Biosemi128.lay';
    cfg.parameter = 'avg';
    cfg.style = 'straight'; cfg.comment = 'no';
    ft_topoplotER(cfg,dev_RR{si})
end
if saveplots;
    set(fh,'color',[1 1 1]);
    filename=[fig_path 'TOPO Each Subject 100ms ' file_in];
    export_fig([filename '.pdf'],'-painters','-q101');
    export_fig([filename '.png'],'-dpng');
end
%% bootstrap average over subjects -for plotting time dmain in a selected channel
varNames       = {'REG', 'RAND','dev_REG','dev_RAND'}; % conditions in dataset
B= 1000;
chan = 87; % channel 87 = C23 = FCz
for c = 1:length(varNames)
    x = squeeze(allData(:,c,chan,:));
    %  Force x to be time*repetitions
    if size(x,2) > size(x,1), x=x'; end
    [mn,sd,bsall]=fBootstrapMean(x,B);
    bsmean(c,:)=mn';
    bsstd(c,:)=sd'; % Bootstrap standard deviation IS the estimate of the sample standard ERROR!
    stderror(c,:) = std(x,0,2)/sqrt(numel(sub_list));
end

%% plot grand average time domain
C=[47 85 151;
    192 0 0;
    0 176 240;
    255 147 0;
    0 176 240;
    255 147 0]/255;

time_s = REG{1}.time;
fh=figure(99);clf
ax1 = gca;

for c = 1:length(condlist)
    plot(time_s,bsmean(c,:),'color',C(c,:))
    hold on
end
%     lh = legend(names(toplot(pi,:)),'fontweight','b','Location','best','Interpreter','none');
% Plot SE of the bootstrap resampling
for c = 1:length(condlist)
    
    b = bsstd(c,:)'; % STDEV
    a = bsmean(c,:)';
    Y = [b+a;  flipud(-b+a)]';
    abscissa = time_s';
    X = [abscissa; flipud(abscissa)];
    h = fill(X,Y,C(c,:),'edgecolor','none','facealpha',0.2); hold on;
    plot(abscissa,a*0,'-k'); % plot zero line
end
set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
title(['\bf ERP at FCz. File: ' file_in]); % \bf = bold font
set(gca,'Box','off')
ylabel('\bf Voltage (uV)');
set(gca,'XLim',[-0.2 0.5]);
if saveplots;
    set(fh,'color',[1 1 1]);
    filename=[fig_path 'FCz_' file_in];
    export_fig([filename '.pdf'],'-painters','-q101');
    export_fig([filename '.png'],'-dpng');
end

%% plot time domain each subject
fh =figure(43); clf
si=1;
for s = sub_list
    for c = 1:length(varNames)
        x = squeeze(allData(si,c,chan,:));
        % plot each subject error response topo
        subplot(subplotRC(1),subplotRC(2),si);
        plot(time_s,x,'color',C(c,:))  ; hold on
    end
    set(gca,'Box','off')
    if si == 1
        ylabel('\bf Voltage (uV)');
    end
    set(gca,'XLim',[-0.2 0.5])
    si=si+1;
end
mtit(['\bf ERP at FCz. File: ' file_in]); % \bf = bold font

if saveplots;
    set(fh,'color',[1 1 1]);
    filename=[fig_path 'FCz_eachSub_' file_in];
    export_fig([filename '.pdf'],'-painters','-q101');
    export_fig([filename '.png'],'-dpng');
end

%% Whole scalp whole epoch cluster stats to find DEVIANCE EFFECT
cfg                     = [];
cfg.correctm            = 'cluster'; % method to correct for MCP
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT'; % dependent samples because within-subjects design
cfg.clusteralpha        = 0.025; % use for thresholding to generate clusters
cfg.minnbchan           = 5;
cfg.clusterstatistic    = 'maxsum'; % statistic that is used to evaluate clusters
cfg.alpha               = 0.025;
cfg.tail                = 0; % 2-tailed: MMN could have neg or pos diffs
cfg.clusterttail        = 0; % 2-tailed for assessing cluster level stats
cfg.numrandomization    = 1000;

% experiment specifics
design=zeros(2,2*nsub);
for c = 1:nsub
    design(1,c) = c;
    design(1,c+nsub) = c;
end
design(2,1:nsub)=1;
design(2,(nsub+1):end)=2;
cfg.design              = design; % logical to identify condition membership of trials
cfg.uvar                = 1;        % "subject" is unit of observation
cfg.ivar                = 2;        % "condition" is the independent variable
cfg.parameter = 'avg';
cfg.neighbours = neighbours;
cfg.channel             = {'EEG'};
cfg.latency             = [0 0.5]; % because we have prior knowledge that the effect will be in this window
cfg.tolerance = 0.0001;
[stat_MMN_pooled] = ft_timelockstatistics(cfg, dev_RR{:}, RR{:});

%% Calculate grand averages for plotting later
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.keepindividual = 'no';
cfg.tolerance = 0.0001;
GA_RR         = ft_timelockgrandaverage(cfg,RR{:});
GA_dev_RR        = ft_timelockgrandaverage(cfg,dev_RR{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_diff_pooled = ft_math(cfg,GA_dev_RR,GA_RR);

%% plot in 2D time by channels format
stat = stat_MMN_pooled;
h=figure(42);clf
imagesc(stat.time,1:length(stat.label),stat.mask.*(stat.posclusterslabelmat-stat.negclusterslabelmat))
title('Significant clusters')
ylabel('Channel')
xlabel('Time (s)')
colormap jet
colorbar
set(gcf,'color',[1 1 1]);
if saveplots
    figure(h)
    set(gcf,'Position',[353   359   560   420]);
    filename=([fig_path 'clusterMap_' ...
        file_in]);
    export_fig([filename '.png'],'-dpng');
end
%% Use cluster mask from cluster stats on DEVIANCE effect (pooled across REG and RAND)
% to give a time range and channel selection for further analysis
stat = stat_MMN_pooled;

%% *************** choose cluster to select as ROI***********************
MMNclusti = -1; % needs to be negative if you want to select a negative cluster!
%% *************** choose cluster to select as ROI***********************

if MMNclusti<0 % ie if negative one selected
    MMNclust_mask = stat.negclusterslabelmat==-MMNclusti;
    [m,rc]=min(stat.stat.*(MMNclust_mask));
    [m2,c]=min(m);
else
    MMNclust_mask = stat.posclusterslabelmat==MMNclusti;
    [m,rc]=max(stat.stat.*(MMNclust_mask));
    [m2,c]=max(m);
end
r=rc(c);
MMNtime = stat.time(c); % gives stat peak time for MMN cluster in seconds
MMNtimeRange = stat.time([find(sum(MMNclust_mask,1),1,'first') find(sum(MMNclust_mask,1),1,'last')]);

%% channel selection
% %METHOD 1: determine which are significant channels at MMN peak time
MMNchanns = find( MMNclust_mask(:,c));
% %METHOD 2: determine all channels within cluster at any point within the cluster -
% % more lax than above
% MMNchanns = find(sum(MMNclust_mask,2)>0);
% % %% Method 3: top 10
% MMNchanns = channord(1:10,2);

MMNchanns_labels = stat.label(MMNchanns);
% MMNchanns_labels = 'C23';

%% plot the MAIN EFFECT of DEVIANCE cluster on a topo
cfg = [];
cfg.channel = 'EEG';
cfg.highlight = 'on';
cfg.style = 'straight';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = 'biosemi128.lay';
cfg.colorbar = 'yes';
cfg.xlim=MMNtimeRange; % time limits in s
cfg.zlim = [-0.7 1]; % from inspection
h = figure(55);clf
cfg.highlightchannel = MMNchanns_labels;
ft_topoplotER(cfg, GA_diff_pooled);
title(['Main effect: Deviance. Shown in timerange ' ...
    sprintf('%0.0f',MMNtimeRange(1)*1000) ' : '  sprintf('%0.0f',MMNtimeRange(2)*1000) ' ms '])
if saveplots
    figure(h)
    set(gcf,'Position',[353   359   560   420]);
    filename=([fig_path 'Topo_Deviance_' ...
        file_in '_cluster_' num2str(MMNclusti)  ]);
    export_fig([filename '.png'],'-dpng');
end

%% Stats to test for REG v RAND within ROI from main effect DEVIANCE (time, channels)
cfg = [];
cfg.channel     = MMNchanns_labels;
cfg.latency     = MMNtimeRange;
cfg.avgovertime = 'yes';
cfg.avgoverchan = 'no';
cfg.parameter   = 'avg';
cfg.method              = 'montecarlo';
cfg.statistic =  'depsamplesT';
cfg.alpha       = 0.025;
cfg.correctm    = 'cluster';
cfg.numrandomization    = 1000;
cfg.tolerance = 0.0001;
load('neighboursBiosemi128.mat');
cfg.neighbours = neighbours;
cfg.design(1,1:2*nsub)  = [ones(1,nsub) 2*ones(1,nsub)];
cfg.design(2,1:2*nsub)  = [1:nsub 1:nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
statRvR = ft_timelockstatistics(cfg,diff_REG{:},diff_RAND{:});

%% plot on topo
% calculate grand average MMNs for each condition
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.keepindividual = 'no';
cfg.tolerance = 0.0001;
GA_diff_REG         = ft_timelockgrandaverage(cfg,diff_REG{:});
GA_diff_RAND        = ft_timelockgrandaverage(cfg,diff_RAND{:});
% calculate the grand average difference for plotting on the topo
cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_diff_RvR = ft_math(cfg,GA_diff_REG, GA_diff_RAND);

% get relevant (significant) values
if isfield(statRvR,'posclusters')&&~isempty(statRvR.posclusters)
    pos_clust_pvals = [statRvR.posclusters(:).prob];
    pos_sig_clust = find(pos_clust_pvals < statRvR.cfg.alpha);
    pos = ismember(statRvR.negclusterslabelmat, pos_sig_clust);
else
    pos=zeros(size(MMNchanns));
end

if isfield(statRvR,'negclusters')&&~isempty(statRvR.negclusters)
    neg_clust_pvals = [statRvR.negclusters(:).prob];
    neg_sig_clust = find(neg_clust_pvals < statRvR.cfg.alpha);
    neg = ismember(statRvR.negclusterslabelmat, neg_sig_clust);
else
    neg=zeros(size(MMNchanns));
end

%%
cfg = [];
cfg.channel = MMNchanns_labels;
cfg.highlight = 'on';
cfg.style = 'straight';
cfg.comment = 'xlim';
cfg.commentpos = 'title';
cfg.layout = 'biosemi128.lay';
cfg.latency = MMNtimeRange;
cfg.xlim = cfg.latency;
cfg.zlim = [-0.7 1] ; %set to saim as colorbar limits for main effect
% plot
h = figure(88);
cfg.highlightchannel = MMNchanns_labels(find(neg|pos));
ft_topoplotER(cfg, GA_diff_RvR);
set(gcf,'color',[1 1 1]);
colorbar;
title(['Interaction effect: Regularity * Deviance in ROI shown in timerange ' ...
    sprintf('%0.0f',MMNtimeRange(1)*1000) ' : '  sprintf('%0.0f',MMNtimeRange(2)*1000) ' ms '])

if saveplots
    figure(h)
    set(gcf,'Position',[ 679   264   632   497]);
    filename=([fig_path 'Topo_Interaction_in_ROI_' ...
        file_in '_cluster_' num2str(MMNclusti)]);
    export_fig([filename '.png'],'-dpng');
end

