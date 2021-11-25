function [ ang_out ] = rad_round( ang )
% limit angle to 0 ~ 2pi
%   Detailed explanation goes here
num = length(ang);
for i = 1:num
    while ang(i) > 2*pi
        ang(i) = ang(i) - 2*pi;
    end
    while ang(i) < -2*pi
        ang(i) = ang(i) + 2*pi;
    end
end
ang_out = ang;
end

