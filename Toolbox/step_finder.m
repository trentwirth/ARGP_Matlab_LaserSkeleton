function [steps_RL] = step_finder(w)
%  STEP_FINDER function to find steps based on density of points
% specifically geared towards finding steps from qualisys marker data

qual_dim_mar_fr = w.skellyPos;
qualiSkellyNames = w.skellyLabels;
fps = mean(w.fps);
walks = w.trials;



%%% Using Markers doesn't work because there are floaters, need to use
%%% skeleton data

rFootID = find(strcmp('RightFoot', qualiSkellyNames));
lFootID = find(strcmp('LeftFoot', qualiSkellyNames));

% get the average point of the forefoot

rFoot = squeeze(qual_dim_mar_fr(:,rFootID,:));
lFoot = squeeze(qual_dim_mar_fr(:,lFootID,:));

figure(1212);
plot3(rFoot(1,:),rFoot(2,:),rFoot(3,:),'.-r');
hold on;
plot3(lFoot(1,:),lFoot(2,:),lFoot(3,:),'.-b');
axis equal

figure(121233);
plot(rFoot(1,:),rFoot(2,:),'.-r');
hold on;
plot(lFoot(1,:),lFoot(2,:),'.-b');
axis equal

figure(4455)
hist3(rFoot(1:2,:)','Ctrs',{-2000:50:2000 -4000:50:10000},'CDataMode','auto','FaceColor','interp')
hold on;
hist3(lFoot(1:2,:)','Ctrs',{-2000:50:2000 -4000:50:10000},'CDataMode','auto','FaceColor','interp')
colorbar
view(2)
axis equal