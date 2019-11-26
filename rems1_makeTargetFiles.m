function [varargout] = rems1_makeTargetFiles(varargin)
% function [varargout] = rems1_makeTargetFiles(varagin)
% Creates .tgt files for subject(s) s and block(s) b
%
% example calls:
%               [G] = rems1_makeTargetFiles([99], [1:12]);
%               [G] = rems1_makeTargetFiles([1:20], [1:12]);
%
% --
% gariani@uwo.ca - 2019.11.19

%% make target files for all subjects and blocks per subject
G = struct();
for s = varargin{1}
    S = struct();
    for b = varargin{2}
        fprintf(1, '\nsubj: %d   block: %d\n', s, b);
        [B] = rems1_target(s, b); % B=block (all trials)
        S = addstruct(S, B); % S=subject (all blocks)
    end
    fprintf(1, '\nSwitch: %03d trials (%2.0f%%)  Repetition: %02d trials (%2.0f%%)\n', sum(S.is_rep==0), (sum(S.is_rep==0)/numel(S.TN))*100, sum(S.is_rep==1), (sum(S.is_rep==1)/numel(S.TN))*100);
    fprintf(1, 'Seq Len=1:  %02d trials    Seq Len=2: %02d trials    Seq Len=3: %02d trials    Seq Len=4: %02d trials\n', sum(S.seq_len==1),  sum(S.seq_len==2),sum(S.seq_len==3),sum(S.seq_len==4));
    fprintf(1, '            (%2.0f%%)                    (%2.0f%%)                    (%2.0f%%)                    (%2.0f%%)\n', (sum(S.seq_len==1)/numel(S.TN))*100,  (sum(S.seq_len==2)/numel(S.TN))*100, (sum(S.seq_len==3)/numel(S.TN))*100, (sum(S.seq_len==4)/numel(S.TN))*100);
    G = addstruct(G, S); % G=group (all subjects)
end
varargout{1} = G;
end

function [varargout] = rems1_target(s, b)
% function [varargout] = rems1_target(s, b)
% This function generates .tgt files, one per block
%
% inputs: vector of subject numbers (s), vector of block numbers (b)
% output: saved filename (fn), block structure (B)

%% define target folder
targetFolder = '../../../../robotcode/projects/SequenceRepetition/rems1/target'; %save tgt files in the right relative path
if ~exist(targetFolder,'dir'); mkdir(targetFolder); end %create target folder if it doesn't already exist

%% experimental details
n_seq_len = 4; % how many sequence length type do you want?
n_seq_exe = 4; % how many exemplars per seq type?
n_seq_trl = 4; % how many trials per seq type?
n_seq = n_seq_len * n_seq_exe; % how many sequences in total?
n_seq_items = 4; % how many elements per sequence (max)?
n_trials = n_seq_len * n_seq_exe * n_seq_trl; % how many trials in a block?
p_rep = 0.50; % probability of same sequence repetition
p_swc = (1 - p_rep) / (n_seq - 1); % probability of switch to each other sequence
T = eye(n_seq); T(T==1) = p_rep; T(T==0) = p_swc; % transition matrix
p1 = ones(1,n_seq) * (1/n_seq); % initial probablity (the same for all sequences)
z = pickRandSeq(p1); % initialize the first random pick

%% pick sequences
allseq = permn(1:8, 4); % all possible sequences taken 4 at a time, with repetitions
seq_pool = zeros(1,size(allseq,2)); % remove "runs" (e.g. 1,2,3, ... 3,2,1, ...)
rowCount = 1;
for r = 1:size(allseq, 1)
    isRun = zeros(1,size(allseq, 2) - 2);
    for c = 1:size(allseq, 2) - 2
        if allseq(r, c)==(allseq(r, c+1) - 1) && allseq(r, c)==(allseq(r, c+2) - 2) % flag descending runs
            isRun(1, c) = 1;
        elseif allseq(r, c)==(allseq(r, c+1) + 1) && allseq(r,c)==(allseq(r, c+2) + 2) % flag ascending runs
            isRun(1, c) = 1;
        elseif any(allseq(r, c)==(allseq(r, c+1:end))) || any(allseq(r, c+1)==(allseq(r, c+2:end))) % flag same-number repetitions
            isRun(1, c) = 1;
        else
            isRun(1, c) = 0;
        end
    end
    if all(isRun==0)
        seq_pool(rowCount,:) = allseq(r,:); % pool of similarly difficult sequences to pick from
        rowCount = rowCount + 1;
    else
        continue
    end
end
randind = randperm( size(seq_pool, 1)); % get random indices
seq_pool = seq_pool(randind, :); % randomize order of sequence pool


%% add sequence length factor
seq_cue = zeros([n_seq_exe, n_seq_items, n_seq_len]); % initialize sequences
c = 1; % initialize counter
L = struct(); % initialize structure
for sl = 1 : n_seq_len
    for se = 1 : n_seq_exe
        f = 0; % initialize while loop flag
        while f==0
            seq_cue(se, :, sl) = seq_pool(c,:); % which sequence from pool?
            seq_idx = randind(c); % which sequence index from pool?
            % check if fisrt element has already appeared for this seq
            % length type. If so, keep searching. If not, pick the seq
            if se>1 && any(seq_cue(se, 1, sl) == seq_cue(1:se-1, 1, sl))
                c = c+1;
            else
                f = 1;
            end
        end
        % remove excess cues in shorter sequences (zero padding)
        seq_cue(se, 1+sl:end, sl) = 0;
        % store info in structure
        L1.seq = seq_cue(se, :, sl);
        L1.idx = seq_idx;
        L1.exe = se;
        L1.len = sl;
        L = addstruct(L, L1);
        c = c+1;
    end
end
% sanity check: how do sequences look? [L.exe, L.len, L.idx, L.seq]

%% fill in dataframe structure B for this block
B = struct(); % initialize block structure
all_trials = []; % initialize empty trial list
for t = 1:n_trials
    B.TN(t,1) = t; % trial number
    B.SN(t,1) = s; % subject number
    B.seq_idx(t,1) = L.idx(find(z,1)); % which sequence index?
    B.seq_len(t,1) = L.len(find(z,1)); % which sequence length?
    B.seq_cue(t,:) = L.seq(find(z,1), :); % which sequence cue?
    if (t > 1) && (B.seq_idx(t, 1) == B.seq_idx(t-1, 1)) % check if this is a repetition or not
        % repetition trial
        B.is_rep(t,1) = 1; % is this a repetition? (0=no; 1=yes)
        B.rep_num(t,1) = B.rep_num(t-1, 1) + 1; % what repetition number is this?
    else
        % switch trial
        B.is_rep(t,1) = 0; % is this a repetition? (0=no; 1=yes)
        B.rep_num(t,1) = 0; % what repetition number is this?
    end
    B.prep_time(t,1) = 2500; % fixed preparation time (in ms)
    B.iti(t,1) = 500; % fixed inter-trial-interval duration (in ms)
    p = T * z; % update probabilities
    z = pickRandSeq(p); % update pick
    % put an upper limit to repetition number (i.e. avoid cases in which
    % you have 6 or more repetitions in a row
    if B.rep_num(t,1) > 4
        while L.idx(find(z,1)) == B.seq_idx(t,1)
            z = pickRandSeq(p); % update pick
        end
    end
    
    % convert numbers to strings for relevant info
    seq_cue_str = num2str(B.seq_cue(t,:), '''%d'', ''%d'', ''%d'', ''%d''');
    prep_time_str = num2str(B.prep_time(t,1), '''%d''');
    iti_str = num2str(B.iti(t,1), '''%d''');
    % form trial info string
    this_trial_str = sprintf('[%s, %s, %s]', seq_cue_str, prep_time_str, iti_str);
    % append list of trial info for .dtp conversion
    if t==1; all_trials = this_trial_str; else; all_trials = [all_trials ', ' this_trial_str]; end %#ok<AGROW>
end
trial_list = sprintf('[%s]', all_trials);
% sanity check: how many repetitions for each condition of interest?
%figure; subplot(131); plt.hist(B.is_rep); subplot(132); plt.hist(B.seq_len); subplot(133); plt.hist(B.rep_num);

%% save structure B as a target file (.tgt) and return output data structure B
outfname = fullfile(targetFolder, sprintf('rems1_s%02d_b%02d', s, b));
dsave(sprintf('%s.tgt', outfname), B);
B.fn = outfname;
varargout{1} = B;

%% export target files in right format for exoskeleton (.dtp)
% read in template .dtp file
tmp = xmlread('template.dtp');

% replace relevant info about trial protocol
textNode = tmp.getElementsByTagName('tptable').item(0).item(0);
textNode.set('Data', trial_list);

% save out
xmlwrite(sprintf('%s.dtp', outfname), tmp);
end

function [z] = pickRandSeq(p)
% function [z] = pickRandSeq(p)
%
% input: p, probability vector of length nSeq (probability of each sequence being selected)
% output: z, indicator variable of length nSeq (which sequence has been selected: 1=yes, 0=no)

% make p a column vector
p = p(:);

% initialize indicator variable z
z = zeros(numel(p), 1);

% define cumulative probability
cump = cumsum(p);

% pick random number from continuous uniform distribution with range [0 1]
r = unifrnd(0,1);

% check in which class the number falls
i = find(r <= cump, 1);

% select and return the appropriate class
z(i) = 1;
end