function [ uo ] = sat( u )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dp=50;dn=-50;
if u>dp
    uo=dp;
elseif u<dn
      uo=dn;
else
    uo=u;
end

end

