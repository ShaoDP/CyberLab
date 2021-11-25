function [ ang_out ] = ang_round( ang )
% limit angle to -180 ~ 180
%   Detailed explanation goes here
num = length(ang);
for i = 1:num
    if(ang(i) > 180)
        ang_out(i) = ang(i) - 360;
    else
        ang_out(i) = ang(i);
    end
end
% ang_out = ang;
end

