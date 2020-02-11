function [varargout] = rems1_analyze(what, varargin)
%% function [varargout] = rems1_analyze(what, varargin)
% Repetition Effects in upper limb Movement Sequences, experiment 1:
% analysis of behavioral data
%
% Operates with vararginoptions, e.g.:
% rems1_analyze('get_data', 'sn',[98]);
%
%
% usage | example calls:
%
%                           rems1_analyze('get_data');                  % extract and save data for each participant
%                           [D] = rems1_analyze('dataset');             % create specific dataset from subj
%                           [D] = rems1_analyze('re_MT');               % repetition effect on movement times (MT)
%                           [D] = rems1_analyze('re_RT');               % repetition effect on reaction times (RT)
%                           [D] = rems1_analyze('re_IRI');              % repetition effect on inter reach intervals (IRI)
%
%
% --
% gariani@uwo.ca - 2019.12.19

%% paths
path_to_data = '/Users/gariani/Documents/data/SequenceRepetition/rems1';
path_to_analyze = fullfile(path_to_data, 'analyze');
if ~exist(path_to_analyze, 'dir'); mkdir(path_to_analyze); end % if it doesn't exist already, create analyze folder

%% globals

% subjects
% incomplete:
% high-error:
subj = {'s98', 's97', 's96', 's95', 's94'};
ns = numel(subj);
subvec = zeros(1,ns);
for i = 1:ns; subvec(1,i) = str2double(subj{i}(2:3)); end

% colors
cbs_red = [213 94 0]/255;
cbs_blue = [0 114 178]/255;
cbs_yellow = [240 228 66]/255;
cbs_pink = [204 121 167]/255;
cbs_green = [0 158 115]/255;
blue = [49,130,189]/255;
lightblue = [158,202,225]/255;
red = [222,45,38]/255;
lightred = [252,146,114]/255;
green = [49,163,84]/255;
lightgreen = [161,217,155]/255;
orange = [253,141,60]/255;
yellow = [254,196,79]/255;
lightyellow = [255,237,160]/255;
purple = [117,107,177]/255;
lightpurple = [188,189,220]/255;
darkgray = [50,50,50]/255;
gray2 = [100,100,100]/255;
gray = [150,150,150]/255;
lightgray = [200,200,200]/255;
silver = [240,240,240]/255;
black = [0,0,0]/255;
white = [255,255,255]/255;

% plot defaults
fs = 20; %default font size for all figures
lw = 4; %4; %default line width for all figures
ms = 10; %12; %default marker size for all figures
rs = 0:3; %0:5; %[0,1,2,3,4,5,6]; %default repetition number subset (full range = 0:10; 0 = full set)
maxRep = 3; % default ceiling level for repetition number (0 = no ceiling)

% styles
style.reset;
style.custom({blue,lightblue,red,lightred,orange,yellow,lightyellow,purple,lightpurple,darkgray,gray,gray2,lightgray,green,lightgreen,black,silver,white,...
    cbs_red,cbs_yellow,cbs_blue,cbs_green,cbs_pink});
isrepsty = style.custom({lightgray, darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
lightgraysty = style.custom({lightgray}, 'markertype','none', 'linewidth',1, 'errorwidth',1, 'errorcap',0, 'linestyle','-');
darkgraysty = style.custom({darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
blacksty = style.custom({black}, 'markertype','none', 'linewidth',lw, 'linestyle','--','errorbars','plusminus', 'errorwidth',lw, 'errorcap',0);%, 'linestyle','none');
% graysty = style.custom({lightgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
% nm1sty = style.custom({cbs_pink, cbs_blue}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
% sffsty = style.custom({lightgray, gray, darkgray}, 'markersize',ms, 'linewidth',lw);
boxplotsty = style.custom({darkgray}, 'markertype','none', 'linewidth',lw);
% isrepstybox = style.custom({gray, black}, 'markersize',lw, 'linewidth',lw);

% legends
isrepleg = {'Switch', 'Repetition'};

%% types of analysis
switch (what)
    case 'get_data' % pre-analysis: extract and save relevant data for each participant
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        % main loop
        for s = 1:numel(sn)
            B = [];
            % get number of blocks for this subj, and filename for each block
            block = dir(fullfile(path_to_data, sprintf('s%02d/*.zip', sn(s))));
            nb = numel(block);
            for b = 1:nb
                fname = fullfile(path_to_data, sprintf('s%02d/%s', sn(s), block(b).name));
                data = sort_trials(zip_load(fname), 'execution');  % loads trials in execution order
                D = data.c3d;
                for t = 1:length(D)
                    % add trial info
                    T.TN(t,1) = t;
                    T.BN(t,1) = b;
                    T.SN(t,1) = sn(s);
                    T.seq_len(t,1) = D(t).TP_TABLE.SEQ_LEN(t);
                    T.seq_cue(t,:) = [D(t).TP_TABLE.TARGET_1(t), D(t).TP_TABLE.TARGET_2(t), D(t).TP_TABLE.TARGET_3(t), D(t).TP_TABLE.TARGET_4(t)];
                    if t>1
                        if all(T.seq_cue(t,:) == T.seq_cue(t-1,:)) && T.is_error(t-1,1) == 0
                            rn = rn + 1;
                            T.is_rep(t,1) = 1;
                            T.rep_num(t,1) = rn;
                        else
                            rn = 0;
                            T.is_rep(t,1) = 0;
                            T.rep_num(t,1) = rn;
                        end
                    else
                        rn = 0;
                        T.is_rep(t,1) = 0;
                        T.rep_num(t,1) = rn;
                    end
                    % add error info (tag error trials)
                    if ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'SEQ_END')))
                        T.is_error(t,1) = 0;
                        T.timing_error(t,1) = 0;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'GO_ANTICIP_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 1;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_ANTICIP_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 1;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_SLOW_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 1;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_WRONG_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 0;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TASK_RESET')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 0;
                    else
                        error('Unknown condition! Check event labels for this trial.');
                    end
                    % if correct trial
                    if T.is_error(t,1) == 0
                        T.prep_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'GO_SIGNAL')),3))*1000 - (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TARGETS_ON')),2))*1000;
                        % get reach times
                        switch T.seq_len(t,1)
                            case 1
                                T.reach1_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_1')),2))*1000;
                                T.reach2_time(t,1) = NaN;
                                T.reach3_time(t,1) = NaN;
                                T.reach4_time(t,1) = NaN;
                            case 2
                                T.reach1_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_1')),2))*1000;
                                T.reach2_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_2')),2))*1000;
                                T.reach3_time(t,1) = NaN;
                                T.reach4_time(t,1) = NaN;
                            case 3
                                T.reach1_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_1')),2))*1000;
                                T.reach2_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_2')),2))*1000;
                                T.reach3_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_3')),2))*1000;
                                T.reach4_time(t,1) = NaN;
                            case 4
                                T.reach1_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_1')),2))*1000;
                                T.reach2_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_2')),2))*1000;
                                T.reach3_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_3')),2))*1000;
                                T.reach4_time(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_4')),2))*1000;
                        end
                        % get mov onset
                        time_go_singal = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'GO_SIGNAL')),2))*1000;
                        time_reach_1 = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'REACH_1')),2))*1000;
                        handX_filt_vel = gradient(D(t).Right_HandX, 0.001);
                        handY_filt_vel = gradient(D(t).Right_HandY, 0.001);
                        TANVEL_filt = sqrt( (handX_filt_vel.^2) + (handY_filt_vel.^2) );
                        [max_vel, max_vel_time] = nanmax(TANVEL_filt(time_go_singal:time_reach_1));
                        time_stamp = find((TANVEL_filt(time_go_singal:time_go_singal+max_vel_time)) <= (0.05 * max_vel)); % Here I take the velocity and find the values above 5%
                        if isempty(time_stamp) || numel(time_stamp) < 100
                            T.is_error(t,1) = 1;
                            T.timing_error(t,1) = 1;
                            T.prep_time(t,1) = NaN;
                            T.reach1_time(t,1) = NaN;
                            T.reach2_time(t,1) = NaN;
                            T.reach3_time(t,1) = NaN;
                            T.reach4_time(t,1) = NaN;
                            T.RT(t,1) = NaN;
                            T.MT(t,1) = NaN;
                            T.TT(t,1) = NaN;
                        else
                            mov_onset = time_go_singal+time_stamp(end); % the first value above 5% peak velocity is your movement onset
                            T.RT(t,1) = mov_onset - time_go_singal;
                            % if RT is negative, it means cue anticipation, hence an error
                            if T.RT(t,1)<=50 || T.RT(t,1)>450
                                T.is_error(t,1) = 1;
                                T.timing_error(t,1) = 1;
                                T.prep_time(t,1) = NaN;
                                T.reach1_time(t,1) = NaN;
                                T.reach2_time(t,1) = NaN;
                                T.reach3_time(t,1) = NaN;
                                T.reach4_time(t,1) = NaN;
                                T.RT(t,1) = NaN;
                                T.MT(t,1) = NaN;
                                T.TT(t,1) = NaN;
                            end
                            % calc mov times
                            T.MT(t,1) = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'SEQ_END')),2))*1000 - mov_onset;
                            T.TT(t,1) = T.RT(t,1) + T.MT(t,1);
                        end
                    else % error trial
                        T.prep_time(t,1) = NaN;
                        T.reach1_time(t,1) = NaN;
                        T.reach2_time(t,1) = NaN;
                        T.reach3_time(t,1) = NaN;
                        T.reach4_time(t,1) = NaN;
                        T.RT(t,1) = NaN;
                        T.MT(t,1) = NaN;
                        T.TT(t,1) = NaN;
                    end
                    % add points info
                    T.points(t,1) = 0;
                end
                % add blocks
                B = addstruct(B, T);
            end
            % save this subj's data in a .mat file
            save(fullfile(path_to_data, sprintf('rems1_s%02d.mat',sn(s))), '-struct', 'B');
        end
        
    case 'dataset' % pre-analysis: create specific dataset from selected participants
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        % main loop
        ds = [];
        for s = 1:numel(sn)
            fprintf('\n%s\n\n', subj{s});
            S = load(fullfile(path_to_data, sprintf('rems1_%s.mat', subj{s}))); % load data structure for each subject
            % add subjects
            ds = addstruct(ds, S);
        end
        % save
        save(fullfile(path_to_analyze, sprintf('rems1_allsub_N=%02d.mat', numel(sn))), '-struct', 'ds'); % save all_data.mat file
        % out
        varargout={ds}; %return main structure
        
    case 're_MT' % repetition effect on MT
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        if nargin>1 % load single subj data
            D = load( fullfile(path_to_data, sprintf('rems1_s%02d.mat', sn)));
        else % load group data
            D = load(fullfile(path_to_analyze, sprintf('rems1_allsub_N=%02d.mat', numel(sn))));
        end
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.rep_num, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Rep Eff MT - subj %02d',sn)); else; figure('Name',sprintf('Rep Eff MT - group (N=%d)',numel(sn))); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep Eff MT
        T = tapply(D, {'SN', 'is_rep', 'seq_len'}, ...
            {D.MT, 'nanmedian', 'name', 'MT'}, ...
            'subset', D.is_error==0);
        
        for ip =  1:4
            subplot(3,4,ip); title(sprintf('Seq Len: %d', ip)); hold on;
            plt.box(T.is_rep, T.MT, 'style',darkgraysty, 'subset',T.seq_len==ip); hold on;
            plt.scatter(T.is_rep+1, T.MT, 'split',T.SN, 'subset',T.seq_len==ip);
            xticklabels({'Switch', 'Repetition'}); ylabel('Movement time (ms)'); set(gca,'fontsize',fs); axis square;
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Delta rep
        T = tapply(D, {'SN', 'seq_len'}, ...
            {D.MT, 'nanmedian', 'name', 'MT'}, ...
            {D.MT, 'nanmedian', 'name', 'MTs', 'subset',D.is_rep==0}, ...
            {D.MT, 'nanmedian', 'name', 'MTr', 'subset',D.is_rep==1}, ...
            'subset',D.is_error==0);
        
        subplot(3,4,5:8);
        plt.box(T.seq_len, ((T.MTs-T.MTr)./T.MT)*100, 'style',darkgraysty);
        hold on;
        plt.scatter(T.seq_len, ((T.MTs-T.MTr)./T.MT)*100, 'split',T.SN, 'regression','none');
        xlabel('Sequence length'); ylabel('Repetition difference (% of MT)');  set(gca,'fontsize',fs); %axis square;
        drawline(0, 'dir','horz', 'linestyle',':');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num
        if maxRep > 0; D.rep_num(D.rep_num >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'rep_num', 'seq_len'}, ...
            {D.MT, 'nanmedian', 'name', 'MT'}, ...
            'subset', D.is_error==0);
        
        for ip =  1:4
            subplot(3,4,8+ip);
            plt.scatter(T.rep_num, T.MT, 'split',T.SN, 'subset',T.seq_len==ip, 'regression','none');
            hold on;
            plt.line(T.rep_num, T.MT, 'style',darkgraysty, 'subset',T.seq_len==ip);
            xlabel('Repetition number'); ylabel('MT (ms)'); set(gca,'fontsize',fs); axis square;
        end
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 're_RT' % repetition effect on RT
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        if nargin>1 % load single subj data
            D = load( fullfile(path_to_data, sprintf('rems1_s%02d.mat', sn)));
        else % load group data
            D = load(fullfile(path_to_analyze, sprintf('rems1_allsub_N=%02d.mat', numel(sn))));
        end
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.rep_num, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Rep Eff RT - subj %02d',sn)); else; figure('Name',sprintf('Rep Eff RT - group (N=%d)',numel(sn))); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep Eff RT
        T = tapply(D, {'SN', 'is_rep', 'seq_len'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.is_error==0);
        
        for ip =  1:4
            subplot(3,4,ip); title(sprintf('Seq Len: %d', ip)); hold on;
            plt.box(T.is_rep, T.RT, 'style',darkgraysty, 'subset',T.seq_len==ip); hold on;
            plt.scatter(T.is_rep+1, T.RT, 'split',T.SN, 'subset',T.seq_len==ip);
            xticklabels({'Switch', 'Repetition'}); ylabel('Reaction time (ms)'); set(gca,'fontsize',fs); axis square;
        end
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Delta rep
        T = tapply(D, {'SN', 'seq_len'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            {D.RT, 'nanmedian', 'name', 'RTs', 'subset',D.is_rep==0}, ...
            {D.RT, 'nanmedian', 'name', 'RTr', 'subset',D.is_rep==1}, ...
            'subset',D.is_error==0);
        
        subplot(3,4,5:8);
        plt.box(T.seq_len, ((T.RTs-T.RTr)./T.RT)*100, 'style',darkgraysty);
        hold on;
        plt.scatter(T.seq_len, ((T.RTs-T.RTr)./T.RT)*100, 'split',T.SN, 'regression','none');
        xlabel('Sequence length'); ylabel('Repetition difference (% of RT)');  set(gca,'fontsize',fs); %axis square;
        drawline(0, 'dir','horz', 'linestyle',':');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep num
        if maxRep > 0; D.rep_num(D.rep_num >= maxRep) = maxRep; end % put a ceiling to nReps
        T = tapply(D, {'SN', 'rep_num', 'seq_len'}, ...
            {D.RT, 'nanmedian', 'name', 'RT'}, ...
            'subset', D.is_error==0);
        
        for ip =  1:4
            subplot(3,4,8+ip);
            plt.scatter(T.rep_num, T.RT, 'split',T.SN, 'subset',T.seq_len==ip, 'regression','none');
            hold on;
            plt.line(T.rep_num, T.RT, 'style',darkgraysty, 'subset',T.seq_len==ip);
            xlabel('Repetition number'); ylabel('RT (ms)'); set(gca,'fontsize',fs); axis square;
        end
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    case 're_IRI' % repetition effect on MT
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        if nargin>1 % load single subj data
            D = load( fullfile(path_to_data, sprintf('rems1_s%02d.mat', sn)));
        else % load group data
            D = load(fullfile(path_to_analyze, sprintf('rems1_allsub_N=%02d.mat', numel(sn))));
        end
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.rep_num, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Rep Eff IRI - subj %02d',sn)); else; figure('Name',sprintf('Rep Eff IRI - group (N=%d)',numel(sn))); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep Eff IRI overall
        T = tapply(D, {'SN', 'is_rep'}, ...
            {D.reach2_time - D.reach1_time, 'nanmedian', 'name', 'IRI1'}, ...
            {D.reach3_time - D.reach2_time, 'nanmedian', 'name', 'IRI2'}, ...
            {D.reach4_time - D.reach3_time, 'nanmedian', 'name', 'IRI3'}, ...
            'subset', D.is_error==0 & D.seq_len>1);
        
        for i = 1:3
            T.IRI(:,i) = eval( sprintf('T.IRI%d', i));
            T = rmfield(T,sprintf('IRI%d', i));
            T.IRInum(:,i) = ones(size(T.SN, 1), 1) * i;
            T.SN(:,i) = T.SN(:,1);
            T.is_rep(:,i) = T.is_rep(:,1);
        end
        T.IRI = reshape(T.IRI, size(T.IRI, 1) * size(T.IRI, 2), 1);
        T.IRInum = reshape(T.IRInum, size(T.IRInum, 1) * size(T.IRInum, 2), 1);
        T.SN = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        T.is_rep = reshape(T.is_rep, size(T.is_rep, 1) * size(T.is_rep, 2), 1);
        
        subplot(121); title('Across sequence length');
        plt.line([T.IRInum], T.IRI, 'split',T.is_rep, 'leg',isrepleg);
        xlabel('IRI number'); ylabel('Inter Reach Interval (ms)'); set(gca,'fontsize',fs); axis square;
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep Eff IRI by seq len
        T = tapply(D, {'SN', 'is_rep', 'seq_len'}, ...
            {D.reach2_time - D.reach1_time, 'nanmedian', 'name', 'IRI1'}, ...
            {D.reach3_time - D.reach2_time, 'nanmedian', 'name', 'IRI2'}, ...
            {D.reach4_time - D.reach3_time, 'nanmedian', 'name', 'IRI3'}, ...
            'subset', D.is_error==0 & D.seq_len>1);
        
        for i = 1:3
            T.IRI(:,i) = eval( sprintf('T.IRI%d', i));
            T = rmfield(T,sprintf('IRI%d', i));
            T.IRInum(:,i) = ones(size(T.SN, 1), 1) * i;
            T.SN(:,i) = T.SN(:,1);
            T.is_rep(:,i) = T.is_rep(:,1);
            T.seq_len(:,i) = T.seq_len(:,1);
        end
        T.IRI = reshape(T.IRI, size(T.IRI, 1) * size(T.IRI, 2), 1);
        T.IRInum = reshape(T.IRInum, size(T.IRInum, 1) * size(T.IRInum, 2), 1);
        T.SN = reshape(T.SN, size(T.SN, 1) * size(T.SN, 2), 1);
        T.is_rep = reshape(T.is_rep, size(T.is_rep, 1) * size(T.is_rep, 2), 1);
        T.seq_len = reshape(T.seq_len, size(T.seq_len, 1) * size(T.seq_len, 2), 1);
        
        subplot(122); title('Divided by sequence length');
        plt.line([T.seq_len T.IRInum], T.IRI, 'split',T.is_rep, 'leg',isrepleg);
        xlabel('IRI number'); ylabel('Inter Reach Interval (ms)'); set(gca,'fontsize',fs); axis square;
        
        plt.match('y');
        
        % out
        D.T=T; %incorporate the sub-structures as fields of main structure
        varargout={D}; %return main structure
        
    otherwise
        error('no such case!')
end