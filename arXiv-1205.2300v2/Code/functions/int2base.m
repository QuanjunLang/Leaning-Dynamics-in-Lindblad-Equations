function dts=int2base(n,b,ndigits)
% convert the integer n into an array 
% ndigits long containing the b-ary expansion of n.

for j=1:ndigits
   dts(ndigits-j+1)=n-b*floor(n/b);
   n=floor(n/b);
end

end