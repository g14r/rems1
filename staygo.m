function etimes = staygo(data,index)
% function etimes = staygo(data,index)
%
%   This is a function to retreive the time people spend moving between
%   targets and stopped at a given target before starting to move. For
%   sequences of various lenghts, it will look up target position and
%   compute the time participants were in that target, then return the
%   time they left each target. It is built to take an entire trial,
%   e.g D = zipload(filename); D = KINARM_add_hand_kinematics(D);
%   D = D.c3d; times = staygo(D(15),15);
%
%   Returns:
%   GS:    Go signal (same as GO_SIGNAL)
%   MO:    Movement Onset (for calculating Reaction Time)
%   En1-4: Time of target n acquisition (same as REACH_N)
%   Ex1-4: Time when participant left each target
%   End:   Time when trial is finished (same as SEQ_END)
%
%
%   Note: this code should run on sequences of lenght = 1 but it won't
%   be informative. Similarly, exit time for the last reach doesn't
%   mean anything, but the script will return a value at this point for
%   when the participant left, e.g. via servo
%
%   Right now the code will not run on error trials.

go = round(data.EVENTS.TIMES(contains(data.EVENTS.LABELS, 'GO_SIGNAL'))*1000);      % Look up the go signal; in the

targets = [data.TARGET_TABLE.X_GLOBAL(1:9),data.TARGET_TABLE.Y_GLOBAL(1:9)].*.01;   % targets are the same for each trial
radii = data.TARGET_TABLE.Logical_radius(1:9)*.01;                                  % initialize according to logical radius

rindex = contains(data.EVENTS.LABELS, 'REACH');                                     % Get Label indices of all reaches performed during trial
rtimes = round(data.EVENTS.TIMES(rindex)*1000);                                     % Look up time each event code was produced & convert to samples

handX_filt_vel = gradient(data.Right_HandX, 0.001);                                 % X and Y components of velocity
handY_filt_vel = gradient(data.Right_HandY, 0.001);
tanvel_filt = sqrt( (handX_filt_vel.^2) + (handY_filt_vel.^2) );                    % compute total velocity
vel = tanvel_filt(go:rtimes(1));                                                    % compute velocity in first reach interval

[~,y] = max(vel > 0.05*max(vel));                                                   % find index where v > 5% of max velocity
mo = y + go;                                                                        % get time w.r.t. beginning of trial

Targs = zeros(4,1);

Targs(1) = data.TP_TABLE.TARGET_1(index);                       %Looking up the target order
Targs(2) = data.TP_TABLE.TARGET_2(index);                       %
Targs(3) = data.TP_TABLE.TARGET_3(index);                       %If target ID > 9 it wasn't shown
Targs(4) = data.TP_TABLE.TARGET_4(index);


E = nan(4,1);                                                   %Entry times (same as rtimes)
L = nan(4,1);                                                   %Leave times (computed here)

for i = 1:length(rtimes)                                        %Loop through targets for reach sequences of length > 2
    
    pos = [data.Right_HandX(rtimes(i)+5:end)-targets(Targs(i),1),data.Right_HandY(rtimes(i)+5:end)-targets(Targs(i),2)];    %position w.r.t. target
    dist = sqrt(sum(pos.^2,2));                                 %distance from target center; when this is > logical radius, we are out of the target
    [~,I] = (min(dist < radii(i)));                             %Find the index of time participant exceeds logical radius
    L(i) = I + 5 + rtimes(i);                                   %Add back time so that we get a time w.r.t. beginning of the trial
    E(i) = rtimes(i);
    
end

etimes.GS = go;
etimes.MO = mo;
etimes.En1 = E(1);                                              %Build the output structure. NAN will be
etimes.En2 = E(2);                                              %retained if that target wasn't used
etimes.En3 = E(3);                                              %there are a maximum of 3 times for a 4-reach sequence
etimes.En4 = E(4);
etimes.Ex1 = L(1);
etimes.Ex2 = L(2);
etimes.Ex3 = L(3);
etimes.Ex4 = L(4);                                              %Exit for target 4 never means anything!
etimes.End = round(data.EVENTS.TIMES(contains(data.EVENTS.LABELS, 'SEQ_END'))*1000);