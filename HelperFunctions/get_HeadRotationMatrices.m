function [rEyeballCenterXYZ, lEyeballCenterXYZ, hCenXYZ, headXhat, headYhat, headZhat] = get_HeadRotationMatrices(markers_fr_mar_dim, markerLabels, headRotMat_row_col_fr, debug)


%% identify 'head' markers, based on standard Qualisys AIM model (for now - 3 March 2020)

headL_xyz = squeeze(markers_fr_mar_dim(:,~cellfun(@isempty, strfind(markerLabels,'HeadL')),:));
headTop_xyz = squeeze(markers_fr_mar_dim(:,~cellfun(@isempty, strfind(markerLabels,'HeadTop')),:));
headR_xyz = squeeze(markers_fr_mar_dim(:,~cellfun(@isempty, strfind(markerLabels,'HeadR')),:));
headFront_xyz = squeeze(markers_fr_mar_dim(:,~cellfun(@isempty, strfind(markerLabels,'HeadFront')),:));

numFrames = length( headR_xyz);

%% define head axes
clear headOrigin_xyz
hCenXYZ(:,1) = mean([headL_xyz(:,1) headR_xyz(:,1)],2); %average of the HeadL and HeadR markers ... after all these years, I still haven't figured out a better way to do this....
hCenXYZ(:,2) = mean([headL_xyz(:,2) headR_xyz(:,2)],2);
hCenXYZ(:,3) = mean([headL_xyz(:,3) headR_xyz(:,3)],2);

%% set head points to origin
z_headL_xyz = headL_xyz - hCenXYZ; % zero head markers to the origin (maybe unnecessary, but conceptually easier for a math-dope like me :P )
z_headR_xyz = headR_xyz - hCenXYZ;
z_headTop_xyz = headTop_xyz - hCenXYZ;
z_headFront_xyz = headFront_xyz - hCenXYZ;


%% define head unit vectors  (this is the part where I pretend to understand linear algebra, based entirely on the YouTube series by 3b1b)
headXhat = normr(z_headFront_xyz);            %+x head axis points from head origin to headFront
headYhat = normr(z_headL_xyz);                %+y head axis points from head origin to headL
headZhat = cross( headXhat, headYhat,2); %+z is the cross product of x_hat and y_hat

%% define eyeballs in head reference frame

distToHeadFront = pdist([0 0 0; z_headFront_xyz(1,:)]); %distance from origin to headFront
distToHeadL = pdist([0 0 0; z_headL_xyz(1,:)]); %distance from origin to headL
distToHeadTop = pdist([0 0 0; z_headTop_xyz(1,:)]); %distance from origin to headTop

z_lEyeballCenterXYZ = headXhat*distToHeadFront*.9 + headYhat*distToHeadL*.5 + headZhat*distToHeadTop*-.5;
z_rEyeballCenterXYZ = headXhat*distToHeadFront*.9 + headYhat*distToHeadL*-.5 + headZhat*distToHeadTop*-.5;

lEyeballCenterXYZ = z_lEyeballCenterXYZ+hCenXYZ;
rEyeballCenterXYZ = z_rEyeballCenterXYZ+hCenXYZ;

%% debug plot
if debug
    figure(27834)
    for fr = [1:20:numFrames]
        
        hx = hCenXYZ(fr,1);
        hy = hCenXYZ(fr,2);
        hz = hCenXYZ(fr,3);
        
        clf
        subplot(121)
        
        
        marX = markers_fr_mar_dim(fr, :, 1);
        marY = markers_fr_mar_dim(fr, :, 2);
        marZ = markers_fr_mar_dim(fr, :, 3);
        
        hold on
        h_mar = plot3(marX, marY, marZ);
        h_mar.LineStyle = 'none';
        h_mar.Marker = '.';
        h_mar.Color = 'k';
        
        plot3(headL_xyz(fr,1), headL_xyz(fr,2), headL_xyz(fr,3), 'bo','DisplayName','HeadL')
        plot3(headR_xyz(fr,1), headR_xyz(fr,2), headR_xyz(fr,3), 'ro','DisplayName','HeadR')
        plot3(headTop_xyz(fr,1), headTop_xyz(fr,2), headTop_xyz(fr,3), 'go','DisplayName','HeadTop')
        plot3(headFront_xyz(fr,1), headFront_xyz(fr,2), headFront_xyz(fr,3), 'o','DisplayName','HeadFront')
        
        plot3(hCenXYZ(fr,1), hCenXYZ(fr,2), hCenXYZ(fr,3), 'mp','DisplayName','HeadOrigin')
        
        headVecScale  = 100;
        plot3([0 headXhat(fr,1)*headVecScale]+hx, [0 headXhat(fr,2)*headVecScale]+hy,  [0 headXhat(fr,3)*headVecScale]+hz, 'r-','LineWidth',2)
        plot3([0 headYhat(fr,1)*headVecScale]+hx, [0 headYhat(fr,2)*headVecScale]+hy,  [0 headYhat(fr,3)*headVecScale]+hz, 'g-','LineWidth',2)
        plot3([0 headZhat(fr,1)*headVecScale]+hx, [0 headZhat(fr,2)*headVecScale]+hy,  [0 headZhat(fr,3)*headVecScale]+hz, 'b-','LineWidth',2)
        
        plot3(rEyeballCenterXYZ(fr,1), rEyeballCenterXYZ(fr,2), rEyeballCenterXYZ(fr,3), 'kh','MarkerSize',12,'MarkerFaceColor','r')
        plot3(lEyeballCenterXYZ(fr,1), lEyeballCenterXYZ(fr,2), lEyeballCenterXYZ(fr,3), 'kh','MarkerSize',12,'MarkerFaceColor','b')
        
        axis equal
        view(3)
        title(['Frame#: ' num2str(fr)])
        %     legend
        
        %%% head zero-ed to origin
        subplot(122)
        
        plot3(z_headL_xyz(fr,1), z_headL_xyz(fr,2), z_headL_xyz(fr,3), 'bo','DisplayName','z_headL')
        hold on
        plot3(z_headR_xyz(fr,1), z_headR_xyz(fr,2), z_headR_xyz(fr,3), 'ro','DisplayName','z_headR')
        plot3(z_headTop_xyz(fr,1), z_headTop_xyz(fr,2), z_headTop_xyz(fr,3), 'go','DisplayName','z_headTop')
        plot3(z_headFront_xyz(fr,1), z_headFront_xyz(fr,2), z_headFront_xyz(fr,3), 'o','DisplayName','z_headFront')
        plot3(0,0,0,'kp')
        
        headVecScale  = 30;
        plot3([0 headXhat(fr,1)]*headVecScale, [0 headXhat(fr,2)]*headVecScale, [0 headXhat(fr,3)]*headVecScale, 'r-','LineWidth',2)
        plot3([0 headYhat(fr,1)]*headVecScale, [0 headYhat(fr,2)]*headVecScale, [0 headYhat(fr,3)]*headVecScale, 'g-','LineWidth',2)
        plot3([0 headZhat(fr,1)]*headVecScale, [0 headZhat(fr,2)]*headVecScale, [0 headZhat(fr,3)]*headVecScale, 'b-','LineWidth',2)
        
        plot3(z_rEyeballCenterXYZ(fr,1), z_rEyeballCenterXYZ(fr,2), z_rEyeballCenterXYZ(fr,3), 'kh','MarkerSize',12,'MarkerFaceColor','r')
        plot3(z_lEyeballCenterXYZ(fr,1), z_lEyeballCenterXYZ(fr,2), z_lEyeballCenterXYZ(fr,3), 'kh','MarkerSize',12,'MarkerFaceColor','b')
        
        axis equal
        view(3)
        xlim([-200 200])
        ylim([-200 200])
        zlim([-200 200])
        grid on 
        %     legend
        
        %%% Here
        headXVect = [ 100; 0; 0 ];
        headYVect = [ 0; 100; 0 ];
        headZVect = [ 0; 0; 100 ];
        
        plot3([0 headXVect(1)],[0 headXVect(2)],[0 headXVect(3)],'r--o');
        plot3([0 headYVect(1)],[0 headYVect(2)],[0 headYVect(3)],'g--o');
        plot3([0 headZVect(1)],[0 headZVect(2)],[0 headZVect(3)],'b--o');

        thisHeadXVectRot = headRotMat_row_col_fr(:,:,fr)*headXVect;
        thisHeadYVectRot = headRotMat_row_col_fr(:,:,fr)*headYVect;
        thisHeadZVectRot = headRotMat_row_col_fr(:,:,fr)*headZVect;
        
%         rX = vrrotvec(headXhat(1),thisHeadXVectRot)

        %plot again

        plot3([0 thisHeadXVectRot(1)],[0 thisHeadXVectRot(2)],[0 thisHeadXVectRot(3)],'r-o');
        plot3([0 thisHeadYVectRot(1)],[0 thisHeadYVectRot(2)],[0 thisHeadYVectRot(3)],'g-o');
        plot3([0 thisHeadZVectRot(1)],[0 thisHeadZVectRot(2)],[0 thisHeadZVectRot(3)],'b-o');

%         dr.drawEulerRotation(gca, euld);
 
        %%% end
        pause(.1)
        drawnow
        
    end
end