%% Fieldtrip  sustained RMS REG vs RAND
% 1. Load in data
% 2. calculate RMS over channels for each condition
% 3. calculate mean and standard deviation over subjects using bootstrap
% resampling
% 4. PLOT: individual subjects RMS time-course
% 5. PLOT: group RMS time-course
% 6. STATS: cluster-based permutation for differences between condiitions
%% NB we pool REG and REGdev; RAND and RANDdev for more statistical power
% _____Rosy Southwell 2017-04_____

%% setup
sub_list = [1:13 15:21];
dir_raw = 'Raw/';
dir_ft = 'FTv3/';
condlist = {'REG10', 'REG10dev','RAND10', 'RAND10dev'};   % condition labels
toplot          = [1 3]; % conditions to plot on same graph
stats = 1;
neeg = 128;
load neighboursBiosemi128.mat;
path_in = [dir_ft 'average_sus/'];
file_in = 'F30ABCIOP';

export = 1; % export plots as pdf/png?
plot_path = 'Plots/FTv3/RMS_SusFT/'; mkdir(plot_path);

[subplotRC,unu] = numSubplots(numel(sub_list)+1);

%% Load in
scount = 0; alldatamat = [];alldatacell = [];
for s = sub_list
    scount = scount+1;
    load([path_in file_in '_s' num2str(s) '.mat']);
    for c = [1 3] % for REG and RAND
        data = eval(['aveData.' condlist{c}]);
        alldatacell{scount,c} = data;
        alldatamat(scount,(c+1)/2,:,:) = data.avg;
    end
    for c = [2 4] % for REG and RAND deviant - add and divide by 2 to get average
        data = eval(['aveData.' condlist{c}]);
        alldatacell{scount,c} = data;
        alldatamat(scount,c/2,:,:) = data.avg + squeeze(alldatamat(scount,c/2,:,:))/2 ;
    end
end
rmsall = squeeze(sqrt(mean(alldatamat.^2,3))); % calculate RMS

%% Calculate bootstrap mean and SD over subjects for all conditions
B= 1000;
for c=1:2
    x = squeeze(rmsall(:,c,:));
    %  Force x to be time*repetitions
    if size(x,2) > size(x,1), x=x'; end
    [mn,sd]=fBootstrapRMS(x,B);
    bsmean(c,:)=mn';
    bsstd(c,:)=sd';
end
% colours
C=[47 85 151; % reg10 pooled
    192 0 0]/255; % rand10 pooled

%% %% Plot individual RMS for each subject
time = data.time;
f1 = figure(1); clf;
set(gcf,'Position',[1 1 1400 400])
for s = 1:numel(sub_list) % loop over all subjects
    subplot(subplotRC(1),subplotRC(2),s)
    for k = 1:2 % loop over all conditions
        plot(time,squeeze(rmsall(s,k,:)),'LineWidth',0.2,'color',C(k,:)); hold on;
    end
    hold off;
    xlim([-0.5 time(end)]);
    ylim([0 max(max(max(rmsall(:,k,:))))*1.2]);
    grid off;
    if s == 1 || numel(sub_list)==1;
        ylabel('\bf RMS (uV)');
    end
    if s == 1 || numel(sub_list)==1;
        xlabel('\bf Time (s)');
    end
end
subplot(subplotRC(1),subplotRC(2),s+1)
plot(1,1,'color',C(1,:)); hold on;
plot(1,1,'color',C(2,:));
axis off
lh = legend(condlist(toplot),'fontweight','b','Location','best');
suptitle(['\bf RMS ' num2str(neeg) ' channs. File: ' file_in]); % \bf = bold font
set(f1,'Position',[10 10 900 600]);

if export
    set(f1,'color',[1 1 1]);
    export_fig(f1,'-painters','-q101',[plot_path 'RMS_eachSub_' file_in '_'  ...
        'inclDev' '.pdf'])
    export_fig(f1,'-dpng',[plot_path 'RMS_eachSub_' file_in '_'  ...
        'inclDev' '.png'])
end
%% --------Plot group RMS + bootstrap (+ stats)----------
if numel(sub_list)>1
    f2 = figure(2); clf
    ax1 = gca;
    set(f2,'Position',[1 1 680 310])
    
    for k = 1:2
        a = bsmean(k,:)';
        plot(time,a,'LineWidth',1.2,'color',C(k,:));
        hold on;
    end
    legh = legend(condlist(toplot),'fontweight','b');
    set(legh,'Location','northwest');
    legend boxoff
    % Plot STD of the bootstrap resampling (= SEM)
    for k = 1:2
        b = 2*bsstd(k,:)'./sqrt(numel(sub_list)); % 2* SE of boostrap resampling
        a = bsmean(k,:)';
        style = 'meanbased';
        if strcmp(style,'zerobased'); % plot bs sd on zero-line
            Y = [b;-flipud(b)]';
        elseif strcmp(style,'meanbased');
            Y = [b+a;flipud(-b+a)]';
        end
        abscissa = time(:);
        X = [abscissa;flipud(abscissa)];
        h = fill(X,Y,C(k,:),'edgecolor','none','facealpha',0.2); hold on;
        %plot(abscissa,a*0,'k'); % plot zero line
    end
    set(gca,'children',flipud(get(gca,'children'))); % Send shaded areas in background
    xlim([-0.5 time(end)]);
    xlabel('\bf Time (s)');
    ylabel('\bf RMS (uV)');
    grid off;
    suptitle(['\bf Group RMS ' num2str(numel(sub_list)) ' subs, '  num2str(neeg) ' channs. File: ' file_in]);
    
    %% Calculate stats using randomization and a clustering approach
    if stats;
        aa = ylim; aa = (aa(2)-aa(1))/30;
        Ystat = 3*aa;
        cond1 = squeeze(rmsall(:,1,:))';
        cond2 = squeeze(rmsall(:,2,:))';
        cfg = [];
        cfg.statistic        = 'ft_statfun_depsamplesT';
        cfg.latency = [0 3.5];
        cfg.numrandomization = 1000;
        cfg.correctm         = 'cluster';
        cfg.method           = 'montecarlo';
        cfg.tail             = 1;
        cfg.alpha            = 0.05;
        cfg.design           = [1:numel(sub_list) 1:numel(sub_list) % subject number
            ones(1,numel(sub_list)) 2*ones(1,numel(sub_list))];  % condition number
        cfg.uvar = 1;        % "subject" is unit of observation
        cfg.ivar = 2;        % "condition" is the independent variable
        cfg.dimord = 'time_subj';
        cfg.dim=[1,numel(time)];
        cfg.connectivity =0;
        
        stat = ft_statistics_montecarlo(cfg, [cond1 cond2],cfg.design);
        
        % Find indices of significant clusters
        pos=[]; neg=[];
        if isfield(stat,'posclusters')
            if ~isempty(stat.posclusters)
                pos_cluster_pvals = [stat.posclusters(:).prob];
                pos_signif_clust = find(pos_cluster_pvals < cfg.alpha);
                poss = ismember(stat.posclusterslabelmat, pos_signif_clust);
                pos = [find(diff([0; poss])==1) find(diff([0; poss])==-1)];
                axes(ax1);
                
                for j = 1:size(pos,1)
                    line([time(pos(j,1)) time(pos(j,2))],[Ystat Ystat],...
                        'LineWidth',4,'Color',C(1,:));hold on
                end
            end
        end
        if isfield(stat, 'negclusters')
            if ~isempty(stat.negclusters)
                neg_cluster_pvals = [stat.negclusters(:).prob];
                neg_signif_clust = find(neg_cluster_pvals <cfg.alpha);
                negs = ismember(stat.negclusterslabelmat, neg_signif_clust);
                neg = [find(diff([0; negs])==1) find(diff([0; negs])==-1)];
                % plot
                axes(ax1);
                for j = 1:size(neg,1)
                    line([time(neg(j,1)) time(neg(j,2))],[Ystat Ystat], ...
                        'LineWidth',4,'Color',C(2,:));
                end
            end
        end
        
        %% store when means differ - use to find when means first diverge stably
        meandiff=mean(cond1,2)-mean(cond2,2);
        Ystat = Ystat+aa;
        hold off;
    end
    
    if export
        set(f2,'color',[1 1 1]);
        export_fig(f2,'-painters','-q101',[plot_path 'group_RMS_' file_in ...
            '_REG10_RAND10_incldev' '.pdf'])
        export_fig(f2,'-dpng',[plot_path 'group_RMS_' file_in ...
            '_REG10_RAND10_incldev' '.png'])
    end
    
    %% Plot the topo of the sustained response
    %1  Calculate grand averages
    for c = 1:4
        cfg = [];
        cfg.channel   = 'all';
        cfg.latency   = 'all';
        cfg.keepindividual = 'no';
        cfg.tolerance = 0.0001;
        GA_sustained{c}         = ft_timelockgrandaverage(cfg,alldatacell{:,c});
    end
    % pool conditions
    cfg = [];
    cfg.operation = 'add';
    cfg.parameter = 'avg';
    % average over REG and REG dev
    GA_REG_inclDev = ft_math(cfg,GA_sustained{1},GA_sustained{3}); % sum
    GA_REG_inclDev.avg = GA_REG_inclDev.avg / 2; % divide by n
    % average over RAND and RAND dev
    GA_RAND_inclDev = ft_math(cfg,GA_sustained{2},GA_sustained{4});
    GA_RAND_inclDev.avg = GA_RAND_inclDev.avg / 2;
    % average over REG and RAND - i.e. all conditions
    GA_seqEvoked = ft_math(cfg, GA_RAND_inclDev, GA_REG_inclDev);
    GA_seqEvoked.avg = GA_seqEvoked.avg / 2;
    
    % difference REG vs RAND
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    GA_RvR_inclDev = ft_math(cfg,GA_REG_inclDev, GA_RAND_inclDev);
    
    %% plot the REG response on topo
    cfg = [];
    cfg.channel = 'EEG';
    cfg.highlight = 'off';
    cfg.style = 'straight';
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    cfg.layout = 'biosemi128.lay';
    cfg.colorbar = 'yes';
    cfg.xlim = [time(pos(1,1)) 3]; % time limits in s; first index is when REG and RAND diverge
    cfg.zlim = [-2 2]; % hard coded to be the same for all conditions for better comparison, but you might need tto adjust this for your data!
    %% plot REG
    h = figure(55);clf
    ft_topoplotER(cfg, GA_REG_inclDev);
    title(['Response to REG including deviant trials. Shown in timerange ' ...
        sprintf('%0.0f',cfg.xlim(1)*1000) ' : '  sprintf('%0.0f',cfg.xlim(2)*1000) ' ms '])
    if export
        figure(h)
        set(gcf,'Position',[353   359   560   420]);
        filename=([plot_path 'Topo_REG_inclDev_' ...
            file_in]);
        export_fig([filename '.png'],'-dpng');
    end
    
    %% plot RAND
    h = figure(56);clf
    ft_topoplotER(cfg, GA_RAND_inclDev);
    title(['Response to RAND including deviant trials. Shown in timerange ' ...
        sprintf('%0.0f',cfg.xlim(1)*1000) ' : '  sprintf('%0.0f',cfg.xlim(2)*1000) ' ms '])
    if export
        figure(h)
        set(gcf,'Position',[353   359   560   420]);
        filename=([plot_path 'Topo_RAND_inclDev_' ...
            file_in]);
        export_fig([filename '.png'],'-dpng');
    end
    
    %% Plot REG - RAND
     h = figure(57);clf
    ft_topoplotER(cfg, GA_RvR_inclDev);
    title(['difference REG vs RAND (includes deviant trials). Shown in timerange ' ...
        sprintf('%0.0f',cfg.xlim(1)*1000) ' : '  sprintf('%0.0f',cfg.xlim(2)*1000) ' ms '])
    if export
        figure(h)
        set(gcf,'Position',[353   359   560   420]);
        filename=([plot_path 'Topo_RvRsus_inclDev_' ...
            file_in]);
        export_fig([filename '.png'],'-dpng');
    end
    
    %% Plot topo of the initial response around 100 Ms, all conditions
     h = figure(58);clf
    ft_topoplotER(cfg, GA_seqEvoked);
    cfg.xlim = [0.05 0.08];
    cfg.zlim = 'maxmin';
    title(['Onset response to all trials. Shown in timerange ' ...
        sprintf('%0.0f',cfg.xlim(1)*1000) ' : '  sprintf('%0.0f',cfg.xlim(2)*1000) ' ms '])
    if export
        figure(h)
        set(gcf,'Position',[353   359   560   420]);
        filename=([plot_path 'Topo_onset_' ...
            file_in]);
        export_fig([filename '.png'],'-dpng');
    end
    
end
