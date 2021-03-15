function trial_plot(data, subj_num, block_num, trial_num)
% This is a function to draw some plots of one or more trials that use the
% same visual targets. The inputs are the summary data structure (D), and
% which subj, block, and trial numbers you want to plot (SN, BN, TN).
%
%usage: trial_plot(D, 1, 1, 5:8); 
%       trial_plot(D, 1, 1, 1:4);

% Select trials you want to plot
D = getrow(data, ismember(data.SN, subj_num) & ismember(data.BN, block_num) & ismember(data.TN, trial_num));

% open figure
figure('Name','Trial plot', 'color',[1 1 1]); 
set(gcf, 'Units','normalized', 'Position',[0.1,0.1,0.8,0.8], 'Resize','off', 'Renderer','painters');        
axis image; 
ax = gca;
ax.Color = [1 1 1];

% Look up the active targets which are specific to this sequence
targs = D.seq_cue(1, D.seq_cue(1,:)<22);                                                 %target = 22 means it was not presented
viscircles(D.tgt_loc{1},D.tgt_rad{1}, 'Color',[.8 .8 .8], 'LineWidth',16, 'LineStyle',':');                                     %plot all targets (grey)
hold on;

% Plot active target(s)
for i = 1:length(targs)
    viscircles(D.tgt_loc{1}(targs(i),:),D.tgt_rad{1}(targs(i)),'Color','r', 'LineWidth',16);                    %active target is black
end

%leg = {'Switch';'Rep 1';'Rep 2';'Rep 3'};
leg = cell(1);
% Loop through trials
for t = 1:numel(trial_num)
        
    %Now plot the kinematics, starting at movement onset and ending when they
    %have been in the target for long enough
    onsett = D.onsets(t,~isnan(D.onsets(t,:)));
    offsett = D.offsets(t,~isnan(D.offsets(t,:)));
    %plot(D.kin_rhx{t}(onsett(1):offsett(end)), D.kin_rhy{t}(onsett(1):offsett(end)), 'o', 'markersize',16, 'markeredgecolor','none', 'linewidth',1, 'markerfacecolor',[.8 .8 .8]-0.2*(t-1)); hold on;
    plot(D.kin_rhx{t}(onsett(1):offsett(end)), D.kin_rhy{t}(onsett(1):offsett(end)), 'linewidth',16, 'color',[.8 .8 .8]-0.2*(t-1)); hold on;
    %plot(D.kin_rhx{t}(offsett(1):offsett(1)+150), D.kin_rhy{t}(offsett(1):offsett(1)+150), 'o', 'markersize',16, 'markeredgecolor','none', 'linewidth',1, 'markerfacecolor','w'); hold on;
    %plot(D.kin_rhx{t}(offsett(2):offsett(2)+150), D.kin_rhy{t}(offsett(2):offsett(2)+150), 'o', 'markersize',16, 'markeredgecolor','none', 'linewidth',1, 'markerfacecolor','w'); hold on;
    %plot(D.kin_rhx{t}(offsett(3):offsett(3)+150), D.kin_rhy{t}(offsett(3):offsett(3)+150), 'o', 'markersize',16, 'markeredgecolor','none', 'linewidth',1, 'markerfacecolor','w'); hold on;
    
%     plot(D.kin_rhx{t}(onsett(2)+150), D.kin_rhy{t}(onsett(2)+150), '+', 'markersize',20, 'color','g', 'linewidth',3); hold on;
%     plot(D.kin_rhx{t}(onsett(3)+150), D.kin_rhy{t}(onsett(3)+150), '+', 'markersize',20, 'color','g', 'linewidth',3); hold on;
%     plot(D.kin_rhx{t}(onsett(4)+150), D.kin_rhy{t}(onsett(4)+150), '+', 'markersize',20, 'color','g', 'linewidth',3); hold on;

%     %Now plot the kinematics, starting at movement onset and ending when they
    %     %have been in the target for long enough
    %     for i = 1:length(targs)
    %         plot(D.kin_rhx{t}(D.onsets(t,i):D.offsets(t,i)), D.kin_rhy{t}(D.onsets(t,i):D.offsets(t,i)), 'LineWidth',2); hold on;
    %     end
    
    leg{t} = num2str(t);
end
legend(leg);
legend boxoff
set(gca,'fontsize',28);
set(gca,'XColor','none');
set(gca,'YColor','none');