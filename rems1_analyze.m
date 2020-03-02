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
%subj = {'s98', 's97', 's96', 's95', 's94', 's93'};
subj = {'s01', 's02', 's03', 's04', 's05'};
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
%isrepsty = style.custom({lightgray, darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
%lightgraysty = style.custom({lightgray}, 'markertype','none', 'linewidth',1, 'errorwidth',1, 'errorcap',0, 'linestyle','-');
darkgraysty = style.custom({darkgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
%blacksty = style.custom({black}, 'markertype','none', 'linewidth',lw, 'linestyle','--','errorbars','plusminus', 'errorwidth',lw, 'errorcap',0);%, 'linestyle','none');
% graysty = style.custom({lightgray}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
% nm1sty = style.custom({cbs_pink, cbs_blue}, 'markersize',ms, 'linewidth',lw, 'errorbars','shade');
% sffsty = style.custom({lightgray, gray, darkgray}, 'markersize',ms, 'linewidth',lw);
%boxplotsty = style.custom({darkgray}, 'markertype','none', 'linewidth',lw);
% isrepstybox = style.custom({gray, black}, 'markersize',lw, 'linewidth',lw);

% legends
isrepleg = {'Switch', 'Repetition'};

%% types of analysis
switch (what)
    case 'get_data' % pre-analysis: extract and save relevant data for each participant
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        %-------------------------------------------------------------------------------------------------------------------------------------
        % main loop
        for s = 1:numel(sn)
            B = [];
            
            %-------------------------------------------------------------------------------------------------------------------------------------
            % get number of blocks for this subj, and filename for each block
            block = dir(fullfile(path_to_data, sprintf('s%02d/*.zip', sn(s))));
            nb = numel(block);
            for b = 1:nb
                T = [];
                fname = fullfile(path_to_data, sprintf('s%02d/%s', sn(s), block(b).name));
                data = sort_trials(zip_load(fname), 'execution');  % loads trials in execution order
                D = data.c3d; D = KINARM_add_hand_kinematics(D);
                % Information added after a few pilot subjs.
                % If absent, create but leave empty
                if ~isfield(D, 'ACH1'); D(1).ACH1 = []; end
                if ~isfield(D, 'ACH2'); D(1).ACH2 = []; end
                if ~isfield(D, 'ACH3'); D(1).ACH3 = []; end
                targets = [D(1).TARGET_TABLE.X_GLOBAL(1:9), D(1).TARGET_TABLE.Y_GLOBAL(1:9)].*.01;   % targets are the same for each trial
                
                for t = 1:length(D)
                    % Information added after a few pilot subjs.
                    % If empty, fill with NaNs
                    if isempty(D(t).ACH1); D(t).ACH1 = NaN; end
                    if isempty(D(t).ACH2); D(t).ACH2 = NaN; end
                    if isempty(D(t).ACH3); D(t).ACH3 = NaN; end
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % add trial info
                    T.TN(t,1) = t;
                    T.BN(t,1) = b;
                    T.SN(t,1) = sn(s);
                    T.seq_len(t,1) = D(t).TP_TABLE.SEQ_LEN(t);
                    T.seq_cue(t,:) = [D(t).TP_TABLE.TARGET_1(t), D(t).TP_TABLE.TARGET_2(t), D(t).TP_TABLE.TARGET_3(t), D(t).TP_TABLE.TARGET_4(t)];
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % tag same sequence repetitions
                    if t>1
                        % every trial but the first
                        if all(T.seq_cue(t,:) == T.seq_cue(t-1,:)) %&& T.is_error(t-1,1) == 0
                            % is repetition
                            rn = rn + 1;
                            T.is_rep(t,1) = 1;
                            T.rep_num(t,1) = rn;
                            T.swc_is_first(t,1) = -1;
                            T.swc_same_len(t,1) = -1;
                        else
                            % is switch
                            rn = 0;
                            T.is_rep(t,1) = 0;
                            T.rep_num(t,1) = rn;
                            % check if it's the first time this sequence
                            % was presented in the block (first switch).
                            % The logic is to check whether repetition
                            % effect is due to novelty (first switch) or
                            % to switching (not first switch).
                            if any( all(T.seq_cue(t,:) == T.seq_cue(1:t-1,:),2) )
                                T.swc_is_first(t,1) = 0;
                            else
                                T.swc_is_first(t,1) = 1;
                            end
                            % check if this switch was preceded by another
                            % sequence of the length. The logic is to
                            % check whether repetition effect is due to
                            % switching sequence type, or switching
                            % sequence length.
                            if T.seq_len(t,1) == T.seq_len(t-1,1)
                                T.swc_same_len(t,1) = 1;
                            else
                                T.swc_same_len(t,1) = 0;
                            end
                        end
                        
                        %-------------------------------------------------------------------------------------------------------------------------------------
                        % look for error trials using error count (ACH3)
                        if ~isnan(D(t).ACH3)
                            if D(t).ACH3(end) == D(t-1).ACH3(end)
                                T.is_error(t,1) = 0;
                            else
                                T.is_error(t,1) = 1;
                            end
                        end
                        
                        %-------------------------------------------------------------------------------------------------------------------------------------
                        % detect shared targets (swc_same_tgt) across
                        % different consecutive sequences (switches)
                        prev_cue = T.seq_cue(t-1, :);   prev_cue(prev_cue > 10) = NaN;
                        this_cue = T.seq_cue(t, :);     this_cue(this_cue > 10) = NaN;
                        T.swc_same_tgt(t,1) = polyval( find(this_cue == prev_cue), 10);
                        
                    else
                        %-------------------------------------------------------------------------------------------------------------------------------------
                        % first trial, always switch, always first, no same
                        rn = 0;
                        T.is_rep(t,1) = 0;
                        T.rep_num(t,1) = rn;
                        T.swc_is_first(t,1) = 1;
                        T.swc_same_len(t,1) = 0;
                        if ~isnan(D(t).ACH3)
                            if D(t).ACH3(end) == 0
                                T.is_error(t,1) = 0;
                            else
                                T.is_error(t,1) = 1;
                            end
                        end
                        T.swc_same_tgt(t,1) = 0;
                    end
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % calculate path length (distance between targets) for each segment and whole sequence
                    T.dist_all(t,:) = nan(1,numel(T.seq_cue(t,:))); % preallocate with NaNs
                    dist = sqrt( sum( ( diff( targets([1,T.seq_cue(t,T.seq_cue(t,:)<10)], :) ) ).^2, 2) ); % calculate distance for successive reaches (always starting from 1, which is home target)
                    T.dist_all(t,1:numel(dist)) = dist'; % store distance info
                    T.dist1(t,1) = T.dist_all(t,1);
                    T.dist2(t,1) = T.dist_all(t,2);
                    T.dist3(t,1) = T.dist_all(t,3);
                    T.dist4(t,1) = T.dist_all(t,4);
                    T.dist_sum(t,1) = nansum(T.dist_all(t,:));
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % add error info (tagging of different error types)
                    if ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'SEQ_END')))
                        T.is_error(t,1) = 0;
                        T.timing_error(t,1) = 0;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'GO_ANTICIP_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 1;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_WRONG_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 0;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_ANTICIP_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 2;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TGT_SLOW_ERROR')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 3;
                    elseif ~isempty(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TASK_RESET')))
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 0;
                    else
                        error('Unknown condition! Check event labels for this trial.');
                    end
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % add movement timing info
                    pt = (round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'GO_SIGNAL')),3)) *1000 ...
                        -(round(D(t).EVENTS.TIMES(contains(D(t).EVENTS.LABELS, 'TARGETS_ON')),2))*1000;
                    if isempty(pt); T.prep_time(t,1) = NaN; else; T.prep_time(t,1) = pt; end
                    % decide whether to use photodiode, or not, for RT
                    % calculation (for pilot data ? SN>90 ? it's not a
                    % good idea because use of photodiode changed
                    % between participants
                    if sn(s) > 90; diode = 0; else; diode = 1; end
                    M = staygo(D(t), t, diode);
                    T = addstruct(T, M); % combine structures to save staygo info
                    % calc reach times
                    T.reach1(t,1) = M.En1 - M.MO;
                    T.reach2(t,1) = M.En2 - M.Ex1;
                    T.reach3(t,1) = M.En3 - M.Ex2;
                    T.reach4(t,1) = M.En4 - M.Ex3;
                    T.reach_all(t,:) = [T.reach1(t,1), T.reach2(t,1), T.reach3(t,1), T.reach4(t,1)];
                    T.reach_sum(t,1) = nansum(T.reach_all(t,:));
                    % calc dwell times
                    T.dwell1(t,1) = M.Ex1 - M.En1;
                    T.dwell2(t,1) = M.Ex2 - M.En2;
                    T.dwell3(t,1) = M.Ex3 - M.En3;
                    T.dwell4(t,1) = M.End - M.En4;
                    T.dwell_all(t,:) = [T.dwell1(t,1), T.dwell2(t,1), T.dwell3(t,1), T.dwell4(t,1)];
                    T.dwell_sum(t,1) = nansum(T.dwell_all(t,:));
                    % calc overall mov times
                    T.RT(t,1) = M.MO - M.GS;
                    T.MT(t,1) = M.End - M.MO;
                    T.TT(t,1) = M.End - M.GS;
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % if RT is negative, or impossibly fast, it means
                    % cue anticipation, hence flag as error. Same for
                    % reach times longer than 1.5s, bad trials
                    if T.RT(t,1)<=50 || any(T.reach_all(t,:)>1500)
                        T.is_error(t,1) = 1;
                        T.timing_error(t,1) = 1;
                    end
                    
                    %-------------------------------------------------------------------------------------------------------------------------------------
                    % add points info
                    T.points_trial(t,1) = D(t).ACH1(end); % how many points for this trial
                    T.points_block(t,1) = D(t).ACH2(end); % cumulative count of points in this block
                    T.errors_block(t,1) = D(t).ACH3(end); % cumulative count of errors in this block
                end
                
                %-------------------------------------------------------------------------------------------------------------------------------------
                % add blocks
                B = addstruct(B, T);
            end
            
            %-------------------------------------------------------------------------------------------------------------------------------------
            % save this subj's data in a .mat file
            save(fullfile(path_to_data, sprintf('rems1_s%02d.mat',sn(s))), '-struct', 'B');
        end
        
    case 'dataset' % pre-analysis: create specific dataset from selected participants
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        
        % main loop
        ds = []; % dataset structure
        for s = 1:numel(sn)
            fprintf('\n%s\n\n', subj{s});
            S = load(fullfile(path_to_data, sprintf('rems1_%s.mat', subj{s}))); % load data structure for each subject
            % add subjects
            ds = addstruct(ds, S);
        end
        
        % save structure in .mat file
        fname = sprintf('rems1_ds_N=%02d', numel(sn));
        save(fullfile(path_to_analyze, [fname '.mat']), '-struct', 'ds'); % save all_data.mat file
        
        % export also as .txt file for later analysis in Python
        dsave(fullfile(path_to_analyze, [fname '.txt']), ds)
        
        % out
        varargout={ds}; %return main structure
        
    case 're_MT' % repetition effect on MT
        sn = subvec;
        vararginoptions(varargin, {'sn'});
        if nargin>1 % load single subj data
            D = load( fullfile(path_to_data, sprintf('rems1_s%02d.mat', sn)));
        else % load group data
            D = load(fullfile(path_to_analyze, sprintf('rems1_ds_N=%02d.mat', numel(sn))));
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
            D = load(fullfile(path_to_analyze, sprintf('rems1_ds_N=%02d.mat', numel(sn))));
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
            D = load(fullfile(path_to_analyze, sprintf('rems1_ds_N=%02d.mat', numel(sn))));
        end
        % select repetitions of interest
        if numel(rs)>1; D = getrow(D, ismember(D.rep_num, rs)); end
        
        % open figure
        if nargin>1; figure('Name',sprintf('Rep Eff IRI - subj %02d',sn)); else; figure('Name',sprintf('Rep Eff IRI - group (N=%d)',numel(sn))); end
        set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');
        
        %-------------------------------------------------------------------------------------------------------------------------------------
        % Rep Eff IRI overall
        T = tapply(D, {'SN', 'is_rep'}, ...
            {D.En2 - D.En1, 'nanmedian', 'name', 'IRI1'}, ...
            {D.En3 - D.En2, 'nanmedian', 'name', 'IRI2'}, ...
            {D.En4 - D.En3, 'nanmedian', 'name', 'IRI3'}, ...
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
            {D.En2 - D.En1, 'nanmedian', 'name', 'IRI1'}, ...
            {D.En3 - D.En2, 'nanmedian', 'name', 'IRI2'}, ...
            {D.En4 - D.En3, 'nanmedian', 'name', 'IRI3'}, ...
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