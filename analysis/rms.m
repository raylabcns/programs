% rms returns the root mean square value.
function y = rms(X)

lengthX = length(X);

y = sqrt(sum(X.^2)/lengthX);

end