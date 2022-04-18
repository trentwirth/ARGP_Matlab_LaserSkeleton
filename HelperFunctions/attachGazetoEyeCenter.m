% the eye structure, eyeball center, camera alignment euler gotten from VOR
% calibration, and the calibration point.
function [GazeXYZ, calibDist] = attachGazetoEyeCenter(EyeData, qual_data, EyeballCenterXYZ, CamAlignEuler, pupLength, calibPoint, vorFrames, gazeMult)


% start at the offset, and count up to the length of the pupil data
    for iFr = 1:pupLength
        headRotMat_row_col_fr(:,:,iFr) = quat2rotm(squeeze(qual_data.Skeletons.OrganizedRotationData(:,6,iFr))'); % 6 = head
    end

Eye_sphCenCam_x = EyeData.sphere_center_x;
Eye_sphCenCam_y = EyeData.sphere_center_y;
Eye_sphCenCam_z = EyeData.sphere_center_z;

Eye_pupCircCen_x = EyeData.circle_3d_center_x;
Eye_pupCircCen_y = EyeData.circle_3d_center_y;
Eye_pupCircCen_z = EyeData.circle_3d_center_z;


qual_dim_mar_fr = qual_data.Skeletons.PositionData(1:pupLength); % this data is cleaned up top and saved in a different variable


confidence = EyeData.model_confidence;
%Recipe for a gaze vector!
GazeXYZ = [Eye_pupCircCen_x Eye_pupCircCen_y Eye_pupCircCen_z] ...  Take your "PupilCircleCenter" (in 3D Eye camera coordinate system, units are mm)
    -[Eye_sphCenCam_x  Eye_sphCenCam_y  Eye_sphCenCam_z];        %Subtract EyeSphereCenter (in eye camera coordiates) >> Origin is now at the center of the EyeSphere in camera coords

%normalize its length
for ll = 1:length(GazeXYZ)
    GazeXYZ(ll,:) = GazeXYZ(ll,:)/norm(GazeXYZ(ll,:));
end


%multiply it by your desired length ;)
calibDist = pdist([EyeballCenterXYZ(vorFrames(1),:); calibPoint])*gazeMult; %myboy pythag;; multiply by 1.25 to make the gaze vector a little longer
GazeXYZ = GazeXYZ*calibDist;

%%%%%%% This part's important - Rotate gaze vector by the rotation matrix, which we get from rCamAlignEuler
%%%%

% get the rotation matrix

GoodPupilRotM = eul2rotm(CamAlignEuler, 'XYZ');


for rr = 1:length(GazeXYZ)

    %     thisET_frame_unrot = camRotMatGuess * [gazeXYZ(rr,1); gazeXYZ(rr,2); gazeXYZ(rr,3)];
    %     thisETframe = headRotMat_row_col_fr(:,:,rr) * thisET_frame_unrot;

    thisET_frame_CalibRot = GoodPupilRotM * [GazeXYZ(rr,1); GazeXYZ(rr,2); GazeXYZ(rr,3)];
    thisETframe = headRotMat_row_col_fr(:,:,rr) * thisET_frame_CalibRot;

    headOrVec(rr,:) =     headRotMat_row_col_fr(:,:,rr) * [2e3; 0;0];

    rGazeXYZ_rotated(rr,:) = thisETframe;

end

% add the eyeball center (in shadow/world coordiates) to translate origin of gaze vector onto the shadow eyeball
rGazeXYZ_rotated(:,1) = rGazeXYZ_rotated(:,1)+ EyeballCenterXYZ(:,1);
rGazeXYZ_rotated(:,2) = rGazeXYZ_rotated(:,2)+ EyeballCenterXYZ(:,2);
rGazeXYZ_rotated(:,3) = rGazeXYZ_rotated(:,3)+ EyeballCenterXYZ(:,3);

% rGazeXYZ_rotated(abs(sum(diff(rGazeXYZ_rotated)'))>40,:) = nan;

confThresh = 0;
rGazeXYZ_rotated(confidence < confThresh,:) = nan;

GazeXYZ = rGazeXYZ_rotated;