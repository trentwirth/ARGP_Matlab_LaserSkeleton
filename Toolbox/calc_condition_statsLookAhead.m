function [condition_struct] = calc_condition_statsLookAhead(condition_struct)
%UNTITLED2 Summary of this function goes here
% calculate mean lookahead distance for each conditoin

for cTrial = 1:condition_struct.numTrials

    % grab stats from look dist 
    condition_struct.trials(cTrial).meanLookDist = mean(condition_struct.trials(cTrial).lookDist)/1000;
    condition_struct.trials(cTrial).medianLookDist = median(condition_struct.trials(cTrial).lookDist)/1000;
    condition_struct.trials(cTrial).sdLookDist = std(condition_struct.trials(cTrial).lookDist)/1000;


end


end