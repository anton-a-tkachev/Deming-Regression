function [Xr,Yr,B,COR,STD] = deming(X,Y,delta)
%DEMING Returns Deming regression
% Copyright Anton Tkachev 2015
%   X, Y - input data vectors
%   delta - X and Y variances ratio (if unknown assume delta = 1)
%   Xr, Yr - regression line data points
%   B - regression line coefficients (slope and intercept)
%   COR - X to Y variables correlation coefficient for Deming regression
%   STD - standard deviation of residues

n = length(X);
Xav = sum(X)/n;
Yav = sum(Y)/n;
Sxx = sum((X - Xav).^2)/(n-1);
Syy = sum((Y - Yav).^2)/(n-1);
Sxy = sum((X - Xav).*(Y - Yav))/(n-1);

B = zeros(2,1);
B(2) = (Syy - delta*Sxx + sqrt((Syy - delta*Sxx)^2 + 4*delta*Sxy^2))/(2*Sxy);
B(1) = Yav - B(2)*Xav;

Xr = X + B(2)/(B(2)^2 + delta)*(Y - B(1) - B(2)*X);
Yr = B(1) + B(2)*Xr;

D = -(B(2)*X - Y + B(1))/sqrt(B(2)^2 + 1);
RSS = sum(D.^2);
TSS = sum((X - Xav).^2 + (Y - Yav).^2);
COR = 1 - RSS/TSS;
STD = sqrt(RSS/(n-1));
end
