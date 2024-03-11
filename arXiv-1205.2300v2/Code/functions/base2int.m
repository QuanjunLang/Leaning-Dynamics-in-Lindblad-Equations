function ints = base2int(N,b)
% convert the array N from a form with D columns containing the b-ary
% expansion of an integer to the standard form of the integer. Outputs a
% column vector.
% 
% Steve Flammia, Sep 2011

D = size(N,2);
base = ones(D,1);

for j=1:D; base(j) = b^(D-j); end;

ints = N*base;

end