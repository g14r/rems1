function traj_3dplot(data,index,dwell,times)

%Visualize a reaching trajcetory. The velocity (in m/s) is plotted on the Z
%axis. Dwell is the dwell time is ms. Times is the structure from staygo.

targets = [data(1).TARGET_TABLE.X_GLOBAL(1:9),data(1).TARGET_TABLE.Y_GLOBAL(1:9)].*.01; %plotting the global targets (constant per participant)
radii = data(1).TARGET_TABLE.Logical_radius(1:9)*.01;  

figure()
viscircles(targets,radii,'Color','0.75,0.75,0.75');                                     %plot all targets (grey)
hold on

Targs = [data.TP_TABLE.TARGET_1(index),data.TP_TABLE.TARGET_2(index),data.TP_TABLE.TARGET_3(index),data.TP_TABLE.TARGET_4(index)];
active = sum(Targs<22);                                                                 %number of targts used


for i = 1:length(Targs)                                                                 %plot active target(s)
    if Targs(i) < 22                                                                    %target = 22 means it was not presented
        viscircles(targets(Targs(i),:),radii(Targs(i)),'Color','k');                    %active target is black
    end
end

vtot = sqrt(data.Right_HandXVel.^2+data.Right_HandYVel.^2);                             %total velocity
plot3(data.Right_HandX(times.MO:times.End),data.Right_HandY(times.MO:times.End),vtot(times.MO:times.End),'Color','k','LineWidth',2)


for i = 1:active-1                                                                      %Plot the parts while participant was not allowed
                                                                                        %to leave the target in red.
    st = eval(strcat('times.En', num2str(i)));
    plot3(data.Right_HandX(st:st+dwell),data.Right_HandY(st:st+dwell),vtot(st:st+200),'Color','r','LineWidth',2)

end

axis square
set(gca,'XColor','none')
set(gca,'YColor','none')
set(gca,'ZColor','none')