function I = matrix_pencil(x, p, Ts)

N = length(x);
Y = hankel(x(1:end-p),x(end-p:end)) ;


Y1 = Y(:,1:end-1); % Y matrix without the last column
Y2 = Y(:,2:end);   % Y matrix without the first column

l = eig(pinv(Y1)*Y2);  % eigenvalues (eq. 16)

% frequency and damping factor as Prony's method
alfa = log(abs(l))/Ts;
freq = atan2(imag(l), real(l))/(2*pi*Ts);

% complex residues (amplitudes and angle phases) as Prony's method
Z = zeros(N, p);
for i = 1:length(l)
    Z(:,i) = transpose(l(i).^(0:N-1));
end
rZ = real(Z);
iZ = imag(Z);

rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));

Z = rZ + 1i*iZ;
h = Z\x;

Amp = abs(h);
theta = atan2(imag(h), real(h));

f = @(x) 0*x;
for t = 1:p
    f = @(x) f(x) + h(t).*l(t).^(x/Ts);
end

%%
I.Amp = Amp;
I.theta = theta;
I.freq = freq;
I.alfa = alfa;
I.r = l;
I.h = h;
I.f = f;
% I.a = a;
end
