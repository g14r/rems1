function traj_plot(data,index,onset_time,end_time)
% This is a function to draw some plots of one or more trials that use the
% same targets (e.g. repetitions). It doesn't require more than one index
% because index is only used to look up the active targets for this token.
% It requires the same number of movement onset and movement end times
% (from staygo 'MO' and staygo 'End') because these define how much of the
% kinematics are plotted.
%
%usage: traj_plot([data],index,[move_time],[end_time])

targets = [data(1).TARGET_TABLE.X_GLOBAL(1:9),data(1).TARGET_TABLE.Y_GLOBAL(1:9)].*.01; %plotting the global targets (constant per participant)
radii = data(1).TARGET_TABLE.Logical_radius(1:9)*.01;                                   %radii of these targets

%Now we look up the active targets which are specific to this sequence
Targs = [data(1).TP_TABLE.TARGET_1(index),data(1).TP_TABLE.TARGET_2(index),data(1).TP_TABLE.TARGET_3(index),data(1).TP_TABLE.TARGET_4(index)];
cmap = cool(5);                                                                         %colors used to differentiate reps

figure()

viscircles(targets,radii,'Color','0.75,0.75,0.75');                                     %plot all targets (grey)
hold on

for i = 1:length(Targs)                                                                 %plot active target(s)
    if Targs(i) < 22                                                                    %target = 22 means it was not presented
        viscircles(targets(Targs(i),:),radii(Targs(i)),'Color','k');                    %active target is black
    end
end

%Now plot the kinematics, starting at movement onset and ending when they
%have been in the target for long enough
for i = 1:length(data)
    plot(data(i).Right_HandX(onset_time(i):end_time(i)),data(i).Right_HandY(onset_time(i):end_time(i)),'LineWidth',2,'Color',cmap(i,:))
end

%Altering some figure properties so it's not quite as hideous
leg = {'Rep 0';'Rep 1';'Rep 2';'Rep 3';'Rep 4'};
legend(leg(1:length(data)))
xlim([targets(1,1)-.2 targets(1,1)+.2])
ylim([targets(1,2)-.2 targets(1,2)+.2])

axis square
set(gca,'XColor','none')
set(gca,'YColor','none')
end
