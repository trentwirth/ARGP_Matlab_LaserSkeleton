function [qualStartTimeRelativeToPupil] = timeSyncPupilQualisys(headVecs, norm_pos_x, norm_pos_y, qualTimestamps_fromZero, pupTimestamps_fromZero)

headXhat = headVecs.xHat;

figure(47328)
clf
plot(headXhat)
title('Click on the start and end of the VOR calibration section (qualisys)')
[qual_VORstartEndFrame, ~] = ginput(2);

figure(47328)
clf
plot(norm_pos_x, 'r.-','DisplayName','Norm_pos_x')
hold on 
plot(norm_pos_y,'b.-','DisplayName','Norm_pos_y')
title('Click on the start and end of the VOR calibration section (pupilLabs)')
legend
[pup_VORstartEndFrame, ~] = ginput(2);


qual_VORstartEndFrame = round(qual_VORstartEndFrame);
qual_VORframes = qual_VORstartEndFrame(1):qual_VORstartEndFrame(2);
qual_clipStartTime = qualTimestamps_fromZero(qual_VORstartEndFrame(1));

pup_VORstartEndFrame = round(pup_VORstartEndFrame);
pup_VORframes = pup_VORstartEndFrame(1):pup_VORstartEndFrame(2);
pup_clipStartTime = pupTimestamps_fromZero(pup_VORstartEndFrame(1));

figure(785)
plot(pupTimestamps_fromZero(pup_VORframes), norm_pos_x(pup_VORframes), 'r.-','DisplayName','Norm_pos_x')
hold on 
% plot(pupTimestamps_fromZero(pup_VORframes), norm_pos_y(pup_VORframes), 'b.-','DisplayName','Norm_pos_y')
plot(qualTimestamps_fromZero(qual_VORframes), headXhat(qual_VORframes,:),'k.-','DisplayName','headXhat')

pupX = detrend(butterLowZero(4,5,300,norm_pos_x(pup_VORframes)));
headX = detrend(butterLowZero(4,5,300,headXhat(qual_VORframes,1)));

pupX = pupX-pupX(1);
headX = headX-headX(1);


plot(abs(pupX))
hold on 
plot(abs(headX))

[pupPks, pupLocs] =findpeaks(abs(pupX),'MinPeakHeight',.05);
[qualPks, qualLocs] =findpeaks(abs(headX),'MinPeakHeight',.05);


qualStartTimeRelativeToPupil = mean((qualTimestamps_fromZero(qualLocs) - pupTimestamps_fromZero(pupLocs)')) + (qual_clipStartTime - pup_clipStartTime);

% pupX = norm_pos_x(pup_VORframes);
% pupY = norm_pos_y(pup_VORframes);
% 
% headX =  headXhat(qual_VORframes,1);
% headY =  headXhat(qual_VORframes,2);
% headZ =  headXhat(qual_VORframes,3);
% 
% pupX = detrend(pupX./max(pupX));
% pupY = detrend(pupY./max(pupY));
% 
% headX = detrend(headX./max(headX));
% headY = detrend(headY./max(headY));
% headZ = detrend(headZ./max(headZ));
% 
% [r_XX, lags_XX] = xcorr(pupX, headX);
% [r_XY, lags_XY] = xcorr(pupX, headY);
% [r_XZ, lags_XZ] = xcorr(pupX, headZ);
% 
% 
% [r_YX, lags_YX] = xcorr(pupY, headX);
% [r_YY, lags_YY] = xcorr(pupY, headY);
% [r_YZ, lags_YZ] = xcorr(pupY, headZ);
% 
% figure(542378)
% plot(lags_XX, r_XX, 'DisplayName','XX')
% hold on 
% plot(lags_XY, r_XY, 'DisplayName','XY')
% plot(lags_XZ, r_XZ, 'DisplayName','XZ')
% 
% 
% plot(lags_YX, r_YX, 'DisplayName','YX')
% hold on 
% plot(lags_YY, r_YY, 'DisplayName','YY')
% plot(lags_YZ, r_YZ, 'DisplayName','YZ')

keyboard