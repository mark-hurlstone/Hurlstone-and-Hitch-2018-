function y = nonnans(y)
% NONNANS non-NaN matrix elements
%  
% DESCRIPTION 
% Returns a column vector of the non-NaN elements of y.
%  
% SYNTAX 
% y = NONNANS(y); 
% y   - vector of data
%
% EXAMPLES 
% y = [NaN NaN 3 4 5 6 NaN];
% z = nonnans(y);
%
% ......................................................................... 

y = y(~isnan(y(:)));
