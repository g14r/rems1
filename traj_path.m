function kdata = traj_path(data,index,dwell,times)

%This function returns the path length covered during a trial while the
%participant was not waiting to be released from a waypoint. It also
%returns a vector of maximum velocities for each reach segment; NAN is
%returned for segments which were not executed in that trial (e.g. segment
%3 in a two-reach trial is NAN). The units of path length are meters and
%the units of velocity are M/S.

%Path length is computed by integrating velocity from
%KINARN_add_hand_kinematics because dt is uniform from this source but not
%from X_Position and Y_position.

%Inputs are on trial, index of that trial, dwell time (this could be
%retrieved from the file data...) and times from staygo.

Targs = [data.TP_TABLE.TARGET_1(index),data.TP_TABLE.TARGET_2(index),data.TP_TABLE.TARGET_3(index),data.TP_TABLE.TARGET_4(index)];

active = sum(Targs<22);                                                                 %number of targts used

vtot = sqrt(data.Right_HandXVel.^2+data.Right_HandYVel.^2);                             %compute total velocity

plen = trapz(vtot(times.MO:eval(strcat('times.En',num2str(active)))))*0.001;            %total path length until enters last target

vmax = nan(4,1);
vmax(1) = max(vtot(times.MO:times.En1));

pdel = 0;                                                                               %Removing distance traveled when participant
for i = 1:active-1                                                                      %had entered the target and wasn't released
    
    st = eval(strcat('times.En', num2str(i)));
    pdel = pdel + trapz(vtot(st:st+dwell))*.001;
    vmax(i+1) = max(vtot(st+dwell:eval(strcat('times.En', num2str(i+1)))));

end

kdata.pathlength = plen-pdel;
kdata.vmax = vmax;

end
