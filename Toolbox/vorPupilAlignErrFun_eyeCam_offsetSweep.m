function [camAlignErr] = vorPupilAlignErrFun_eyeCam_offsetSweep(vData, camAlignGuess_xyz_a)
%VORPUPILALIGNERRFUN Summary of this function goes here
%   Detailed explanation goes here

% camQuatGuess = quaternion.eulerangles('123',camAlignGuess_xyz_a(1),camAlignGuess_xyz_a(2),camAlignGuess_xyz_a(3));
% camRotMatGuess = camQuatGuess.RotationMatrix;

camRotMatGuess = eul2rotm(camAlignGuess_xyz_a(1:3), 'XYZ');

vorFrames = vData.vorFrames;

offset = round(vData.initialOffsetGuess*camAlignGuess_xyz_a(4));

%%unpack vData

%%%Qualisys Data
headRotMat_row_col_fr = vData.headRotMat_row_col_fr;

eyeballCenterXYZ = vData.eyeballCenterXYZ(vorFrames,:);

qual_dim_mar_fr = vData.qual_dim_mar_fr(:,:,vorFrames);

calibPoint = vData.calibPoint;

%%%Pupil Data
eye_sphCenCam_x = vData.eye_sphCenCam_x(vorFrames-offset);
eye_sphCenCam_y = vData.eye_sphCenCam_y(vorFrames-offset);
eye_sphCenCam_z = vData.eye_sphCenCam_z(vorFrames-offset);

eye_pupCircCen_x = vData.eye_pupCircCen_x(vorFrames-offset);
eye_pupCircCen_y = vData.eye_pupCircCen_y(vorFrames-offset);
eye_pupCircCen_z = vData.eye_pupCircCen_z(vorFrames-offset);

confidence = vData.confidence(vorFrames-offset);

%Recipe for a gaze vector!
gazeXYZ = [eye_pupCircCen_x eye_pupCircCen_y eye_pupCircCen_z] ...  Take your "PupilCircleCenter" (in 3D Eye camera coordinate system, units are mm)
    -[eye_sphCenCam_x  eye_sphCenCam_y  eye_sphCenCam_z];        %Subtract EyeSphereCenter (in eye camera coordiates) >> Origin is now at the center of the EyeSphere in camera coords

%normalize its length
%normalize its length
for ll = 1:length(gazeXYZ)
    gazeXYZ(ll,:) = gazeXYZ(ll,:)/norm(gazeXYZ(ll,:));
end

%multiply it by your desired length ;)
calibDist = pdist([eyeballCenterXYZ(1,:); calibPoint]); %myboy pythag
gazeXYZ = gazeXYZ*calibDist;




%%
%%%%%%% This part's important - Rotate gaze vector by [this guess at the proper alignment matrix], prior to resituating  the origin on on the eyeball
%%%%


for rr = 1:length(gazeXYZ)
    
    thisET_frame_unrot = camRotMatGuess * [gazeXYZ(rr,1); gazeXYZ(rr,2); gazeXYZ(rr,3)];
    thisETframe = headRotMat_row_col_fr(:,:,rr) * thisET_frame_unrot;
    
    headOrVec(rr,:) =     headRotMat_row_col_fr(:,:,rr) * [2e3; 0;0];
    
    gazeXYZ(rr,:) = thisETframe;
    
end




% add the eyeball center (in shadow/world coordiates) to translate origin of gaze vector onto the shadow eyeball
gazeXYZ(:,1) = gazeXYZ(:,1)+ eyeballCenterXYZ(:,1);
gazeXYZ(:,2) = gazeXYZ(:,2)+ eyeballCenterXYZ(:,2);
gazeXYZ(:,3) = gazeXYZ(:,3)+ eyeballCenterXYZ(:,3);

gazeXYZ(abs(sum(diff(gazeXYZ)'))>40,:) = nan;

confThresh = 0.8;% was .8, bumped it up to .9
gazeXYZ(confidence < confThresh,:) = nan;


%% %%%%% Integrate "correctAlignmentError" code into this guy?

%%%%%%%% Dont need to do this, because we don't have an alignment error for
%%%%%%%% the qualisys data. yet. hopefully, ever.

%         corrAlignTheta = camAlignGuess_xyz_a(4);
% corrAlignTheta = vData.corrAlignTheta;
% 
% % center g on comXYZ (comXYZ reference frame)
% g_z(:,1) = gazeXYZ(:,1)-comXYZ(:,1);
% g_z(:,2) = gazeXYZ(:,3)-comXYZ(:,3);
% 
% %center camXYZ on comXYZ (comXYZ reference frame
% c_z(:,1) = eyeballCenterXYZ(:,1)-comXYZ(:,1);
% c_z(:,2) = eyeballCenterXYZ(:,3)-comXYZ(:,3);
% 
% [gTheta, gRho] = cart2pol(g_z(:,1), g_z(:,2));
% [g_z(:,1), g_z(:,2)] = pol2cart(gTheta-corrAlignTheta, gRho); %rotate by -CorrAlignTheta
% 
% [cTheta, cRho] = cart2pol(c_z(:,1), c_z(:,2));
% [c_z(:,1), c_z(:,2)] = pol2cart(cTheta-corrAlignTheta, cRho); %rotate by -CorrAlignTheta
% 
% % revert gaze to right coord system
% gazeXYZ(:,1) = g_z(:,1)+comXYZ(:,1);
% gazeXYZ(:,3) = g_z(:,2)+comXYZ(:,3);
% 
% % revert camXYZ to right coord system
% eyeballCenterXYZ(:,1) = c_z(:,1)+comXYZ(:,1);
% eyeballCenterXYZ(:,3) = c_z(:,2)+comXYZ(:,3);



%% calculate error for this alignment

%%%error as - minimizing distance between gazeXYZ and calib point (inherits
%%%error from calib point estimation)
for cc = 1:length(gazeXYZ)
    % error for each frame == distance between gaze/ground intersection and the calibration point
    err(cc) = sqrt( (calibPoint(1) - gazeXYZ(cc,1))^2 + (calibPoint(2) - gazeXYZ(cc,2))^2 +(calibPoint(3) - gazeXYZ(cc,3))^2 );
    
    
    %     err(cc) = sqrt( (calibPoint(1) - gazeGroundIntersectionXYZ(cc,1))^2 + (calibPoint(2) - gazeGroundIntersectionXYZ(cc,2))^2 +(calibPoint(3) - gazeGroundIntersectionXYZ(cc,3))^2 );
% tried error by distance, didn't change much. Average distance though was
% 10cm, bad? probably? Idk. TDW 10/17/2021
%    err(cc) = pdist([gazeXYZ(cc,:); calibPoint]);
    
end


% %%%error as - minimizing velocity of end of the gaze vector (from Rakshit
% %%%Kothari's version of this algorithm)
% err = nansum(abs(diff(gazeXYZ)),2);
% 
camAlignErr = nansum(err)/length(err);

% Trent and Jon quick try: it didn't work!
% camAlignErr = sum(nansum(abs(diff(gazeXYZ))))/length(gazeXYZ);

% if vData.dataType == 2
%     camAlignErr = camAlignErr + abs(corrAlignTheta- vData.rCorrAlignTheta); %if this is the left eye, add a small cost to having the alignment correction value be different from the right eye
% end
%% debugs ar the bets bugs
if vData.plotDebug
    
    %%%%% make sphere thingy fr eyeball guys
    sphRes = 20;
    r =12; %your eyeballs have a radius of 12 mm
    [th, phi] = meshgrid(linspace(0, 2*pi, sphRes), linspace(-pi, pi, sphRes));
    [x1,y1,z1] = sph2cart(th, phi, r);
    
    normScale = 1;
    plotSkel = true;
    
%     lLeg = [2 3 4 5 6 7 5];
%     rLeg = [2 8 9 10 11 12 10];
%     tors = [2 13 14 15 26 27 28];
%     lArm = [15 16 17 26 17 18 19 20];
%     rArm = [15 21 22 26 22 23 24 25];
        lLeg = [15 16 17 18];
        rLeg = [19 20 21 22];
        lArm = [ 7 8 9 10 ];
        rArm = [11 12 13 14];
        tors = [1 2 3 4 5 6];
    
     if vData.dataType == 2
         figureAdd = 1;
     else 
         figureAdd = 0;
     end
    
    figure(987+figureAdd);clf
    

        frames = 1;
    
    for ii = frames%:length(eye_sphCenCam_x);
        subplot(1,2,1)
        cla
        
        %%eyeball centers in shadow coordinats(not to be confused with "rEye_sphCen_x", which are in pupil camera coords)
        ex = eyeballCenterXYZ(ii,1);
        ey = eyeballCenterXYZ(ii,2);
        ez = eyeballCenterXYZ(ii,3);
        
        
        
        %%pull out the l and r eye sphere centers for this frame
        cx = eye_sphCenCam_x(ii);
        cy = eye_sphCenCam_y(ii);
        cz = eye_sphCenCam_z(ii);
        
        
        
        
        % right eye
        e1 =  mesh(x1+ex, y1+ey, z1+ez);
        e1.EdgeColor = 'k';
        e1.EdgeAlpha = 0.1;
        hold on
        
        
        %%% Plot circular patch for pupil - centered on pupilNorm (code jacked from - https://www.mathworks.com/matlabcentral/fileexchange/26588-plot-circle-in-3d)
        thisPupCenter = [eye_pupCircCen_x(ii)-cx eye_pupCircCen_y(ii)-cy eye_pupCircCen_z(ii)-cz]*1e3 ;
        thisPupNormal = thisPupCenter*normScale;
        
        %%%%
        switch vData.dataType
            case 1
                gzCol = 'r';
            case 2
                gzCol = 'b';
            case 3
                gzCol = 'k';
        end
        
        plot3([0+ex nanmean(gazeXYZ(:,1))],...
            [0+ey nanmean(gazeXYZ(:,2))],...
            [0+ez nanmean(gazeXYZ(:,3))],'-s','LineWidth',2,'Color',gzCol)
        
        
        gz1 = plot3(gazeXYZ(:,1),...
            gazeXYZ(:,2),...
            gazeXYZ(:,3),'o','Color',gzCol);


        
        %         plot3(gazeGroundIntersectionXYZ(:,1),...
        %             gazeGroundIntersectionXYZ(:,2),...
        %             gazeGroundIntersectionXYZ(:,3),'m.')
        %
        %         plot3(gazeGroundIntersectionXYZ(ii,1),...
        %             gazeGroundIntersectionXYZ(ii,2),...
        %             gazeGroundIntersectionXYZ(ii,3),'kp','MarkerFaceColor','r','MarkerSize',12)
        %
        %%%plot yer calibration point
        plot3(calibPoint(1), calibPoint(2), calibPoint(3), 'mh', 'MarkerFaceColor','m','MarkerSize',24)
        
        %%%Plotcherself a skeleetoon
        plot3(qual_dim_mar_fr(1,1:22,ii),qual_dim_mar_fr(2,1:22,ii),qual_dim_mar_fr(3,1:22,ii),'ko','MarkerFaceColor','k','MarkerSize',4)
        hold on
        
        plot3(qual_dim_mar_fr(1,lLeg,ii),qual_dim_mar_fr(2,lLeg,ii),qual_dim_mar_fr(3,lLeg,ii),'b','LineWidth',2)
        plot3(qual_dim_mar_fr(1,lArm,ii),qual_dim_mar_fr(2,lArm,ii),qual_dim_mar_fr(3,lArm,ii),'b','LineWidth',2)
        plot3(qual_dim_mar_fr(1,rLeg,ii),qual_dim_mar_fr(2,rLeg,ii),qual_dim_mar_fr(3,rLeg,ii),'r','LineWidth',2)
        plot3(qual_dim_mar_fr(1,rArm,ii),qual_dim_mar_fr(2,rArm,ii),qual_dim_mar_fr(3,rArm,ii),'r','LineWidth',2)
        plot3(qual_dim_mar_fr(1,tors,ii),qual_dim_mar_fr(2,tors,ii),qual_dim_mar_fr(3,tors,ii),'g','LineWidth',2)
        
%         plot3(shadow_fr_mar_dim(ii,lLeg,1),shadow_fr_mar_dim(ii,lLeg,2),shadow_fr_mar_dim(ii,lLeg,3),'b','LineWidth',2)
%         plot3(shadow_fr_mar_dim(ii,rLeg,1),shadow_fr_mar_dim(ii,rLeg,2),shadow_fr_mar_dim(ii,rLeg,3),'r','LineWidth',2)
%         plot3(shadow_fr_mar_dim(ii,tors,1),shadow_fr_mar_dim(ii,tors,2),shadow_fr_mar_dim(ii,tors,3),'g','LineWidth',2)
%         plot3(shadow_fr_mar_dim(ii,lArm,1),shadow_fr_mar_dim(ii,lArm,2),shadow_fr_mar_dim(ii,lArm,3),'b','LineWidth',2)
%         plot3(shadow_fr_mar_dim(ii,rArm,1),shadow_fr_mar_dim(ii,rArm,2),shadow_fr_mar_dim(ii,rArm,3),'r','LineWidth',2)
        
        
        
        bx =   qual_dim_mar_fr(1,1,ii);
        by =   qual_dim_mar_fr(2,1,ii);
        bz =   qual_dim_mar_fr(3,1,ii);
        
        grndSize = 3e3;
        g_x = meshgrid(0-grndSize:100:0+grndSize)+calibPoint(1);
        g_y = meshgrid(0-grndSize:100:0+grndSize)'+calibPoint(2);
        g_z = ones(size(g_x))+calibPoint(3);
        
        surface(g_x, g_y, g_z,'FaceColor','none','EdgeColor','k')
        
        
        
        axis equal
        title(num2str(ii))
        %     set(gca,'CameraUpVector',[0 1 0])
        xlabel('x');ylabel('y'); zlabel('z');
        axis([-grndSize+bx grndSize+bx -grndSize+by grndSize+by -grndSize+bz grndSize+bz])
        
        a = gca;
%         a.CameraTarget = [cx cy cz]; %point figure 'camera' at COM
%         a.CameraPosition = a.CameraTarget + [-1800 1800 2000]; %set camera position
%         a.CameraViewAngle = 80;
%         a.CameraUpVector = [ 0 1 0];
        %         a.Position = [0 0 1 1];
        
        hold off
        
        subplot(122)
        hold on;
        ylim([-1000 1000]);
        title('plot 122')
        plot(gazeXYZ(:,1)-calibPoint(1),'b.-');
        plot(gazeXYZ(:,2)-calibPoint(2),'r.-');
        plot(gazeXYZ(:,3)-calibPoint(3),'g.-');
        
%         rx = refline(0,calibPoint(1)); rx.Color = 'r';
%         ry = refline(0,calibPoint(2)); ry.Color = 'r';
%         rz = refline(0,calibPoint(3)); rz.Color = 'r';

        drawnow
    end
end