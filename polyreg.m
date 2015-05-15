function [r, a, c, s] = polyreg(x,y,p)
%POLYREG Polynomial regression
% Copyright Anton Tkachev 2015
%   Inputs
%   x, y - dataset vectors
%   p - polynomial order
%   Outputs
%   r - regression points
%   a - polynomial coefficients array
%   c - correlation between the regression line and the points of dataset
%   s - residual standard deviation calculated from regression line

n = length(x);
X = zeros(n,p+1);
for i=1:n
     X(i,1) = 1;
    for j=1:p
        X(i,j+1) = x(i)^j;
    end
end

a = (X'*X)\X.'*y;
r = X*a;
s = std(y - r);
c = 1 - (s/std(y))^2;

end

