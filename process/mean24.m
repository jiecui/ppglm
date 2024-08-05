% Author: Timothée Proix
% License: GPL-3.0-only

function [Y] = mean24(X)
	% obtain the mean on consecutive non-overlapping windows of 24 samples
	% hours -> days
	% X = (time, whatever) with time multiple of 24

	for iDay = 1:size(X, 1)/24
		Y(iDay, :) = mean(X((iDay-1)*24 + 1:iDay*24, :), 1);
	end