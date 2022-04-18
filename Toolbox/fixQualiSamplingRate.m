function[qual_data] = fixQualiSamplingRate(qual_data, prefRate,unixStartTime)
% FIXQUALISAMPLINGRATE Summary of this function goes here
% resample qualisys data at a constant rate
% prefRate = desired sampling rate in Hz

qual_trajData  = qual_data.Trajectories.Labeled.Data;
% Type isn't used, and it's not a double, so we're not going to worry about
% it....
% qual_trajType = qual_data.Trajectories.Labeled.Type;
qual_skelPos = qual_data.Skeletons.PositionData;
qual_skelRot  = qual_data.Skeletons.RotationData;

trial_length_seconds = qual_data.Frames/qual_data.FrameRate;

% 

startTime = unixStartTime;
endTime = unixStartTime + trial_length_seconds;

fakeTimestamp = (startTime-startTime):1/qual_data.FrameRate:(endTime-startTime);
fakeTimestamp = fakeTimestamp(1:end-1); % end-1 sort of issue? TDW: had to change this depending on the recorded condition...

qual_data.Timestamp_Resampled = (startTime-startTime):1/prefRate:(endTime-startTime);
% numFrames = length(desTimestamp);

%
[s1, ~, ~] = size(qual_trajData);

if length(fakeTimestamp) ~= size(qual_trajData,3)
print('UH OH Something is wrong!')
end

for count = 1:s1

    temp_Data = squeeze(qual_trajData(count,:,:));

    temp_resamp = resample(temp_Data', fakeTimestamp, prefRate, 'pchip')';
    
    resamp_qual_trajData(count,:,:) = temp_resamp;

end

qual_data.Trajectories.Labeled.Data_Resampled = resamp_qual_trajData;

[s2, ~, ~] = size(qual_skelPos);

for count2 = 1:s2

    temp_Data2 = squeeze(qual_skelPos(count2,:,:));

    temp_resamp2 = resample(temp_Data2', fakeTimestamp, prefRate, 'pchip')';
    
    resamp_qual_skelData(count2,:,:) = temp_resamp2;

end

qual_data.Skeletons.PositionData_Resampled = resamp_qual_skelData;

[s3, ~, ~] = size(qual_skelRot);

for count3 = 1:s3

    temp_Data3 = squeeze(qual_skelRot(count3,:,:));

    temp_resamp3 = resample(temp_Data3', fakeTimestamp, prefRate, 'pchip')';
    
    resamp_qual_skelRot(count3,:,:) = temp_resamp3;

end

qual_data.Skeletons.RotationData_Resampled = resamp_qual_skelRot;

