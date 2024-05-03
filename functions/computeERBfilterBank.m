function [S_erb] = computeERBfilterBank(S, fb)
%COMPUTEERBFILTERBANK Summary of this function goes here
%   Detailed explanation goes here

S_erb = max(fb.*repmat(S', size(fb, 1), 1), [], 2);
S_erb(S_erb<0) = 0;

end

