function I = polynomial_method(x, p, Ts, method)




CLASSIC = 0;
LS = 1;
TLS = 2;
ID = 'ID';
RKHS = 'RKHS';
RKHS_h0 = 'RKHS_h0';
ID_h0 = 'ID_h0';
LS_h0 = 'LS_h0';
CLASSIC_h0 = 'CLASSIC_h0';
N = length(x);
if strcmpi(method, 'classic')
    if N ~= 2*p
        disp ('ERROR: length of x must be 2*p samples in classical method.');
        Amp = [];
        alfa = [];
        freq = [];
        theta = [];
        
        I.Amp = Amp;
        I.theta = theta;
        I.freq = freq;
        I.alfa = alfa;
        
        return;
    else
        solve_method = CLASSIC;
    end
elseif strcmpi(method,'LS')
    solve_method = LS;
elseif strcmpi(method,'TLS')
    solve_method = TLS;
elseif strcmpi(method,'RKHS')
    solve_method = RKHS;
elseif strcmpi(method,'RKHS_h0')
    solve_method = RKHS_h0;
elseif strcmpi(method,'ID_h0')
    solve_method = ID_h0;
elseif strcmpi(method,'ID')
    solve_method = ID;
elseif strcmpi(method,'CLASSIC_h0')
    solve_method = CLASSIC_h0;
elseif strcmpi(method,'LS_h0')
    solve_method = LS_h0;
else
    disp('ERROR: error in parsing the arqument "method".');
    Amp = [];
    alfa = [];
    freq = [];
    theta = [];
    
    I.Amp = Amp;
    I.theta = theta;
    I.freq = freq;
    I.alfa = alfa;
    return;
    
end

%% step 1 : solve the coefficients of polynomial
T = toeplitz(x(p:N-1),x(p:-1:1));
switch solve_method
    case {CLASSIC, LS, CLASSIC_h0, LS_h0}
        a = -T\x(p+1:N);
    case TLS
        a = tls(T,-x(p+1:N));
    case{RKHS, ID}
        [~, a] = L_curve(T'*T, -T'*x(p+1:N), solve_method, 1);
    case{RKHS_h0, ID_h0}
        method = solve_method(1:end-3);
        [~, a] = L_curve(T'*T, -T'*x(p+1:N), method, 1);
end

%% check for indeterminate forms
indeterminate_form = sum(isnan(a) | isinf(a));
if (indeterminate_form)
    Amp = []; alfa = []; freq = []; theta = [];
    return;
end

%% step 2
c = transpose([1; a]);
r = roots(c);
alfa = log(abs(r)) / Ts;
freq = atan2(imag(r),real(r)) /(2*pi*Ts);

%In case alfa equals to +/-Inf the signal will not be recovered for n=0
%(Inf*0 = Nan). Making alfa = +/-realmax that indeterminace will be solved

alfa(isinf(alfa)) = realmax*sign(alfa(isinf(alfa)));



%% step 3
switch solve_method
    case CLASSIC
        len_vandermonde = p; % exact case(N = 2p) find h with p samples
    case {LS, ID, RKHS}
        len_vandermonde = N; % overdetermined case (N>2p) find h with N samples
    case {ID_h0, RKHS_h0}
        len_vandermonde = N; % overdetermined case (N>2p) find h with N samples
    case {LS_h0, CLASSIC_h0}
        len_vandermonde = N; % overdetermined case (N>2p) find h with N samples
    case TLS
        len_vandermonde = N;  %overdetermined case (N>2p) find h with N samples
end

%% play with the roots
n = length(r);
rr = r;
for i = 1:n
    if abs(r(i)) >=1
        rr(i) = r(i)/abs(r(i));
    end
end

r = rr;






%%
Z = zeros(len_vandermonde, p);
for i=1:length(r)
    Z(:,i) = transpose(r(i) .^(0:len_vandermonde-1));
end

rZ = real(Z);
iZ = imag(Z);
% here Inf values are substituted by realmax values
rZ(isinf(rZ))=realmax*sign(rZ(isinf(rZ)));
iZ(isinf(iZ))=realmax*sign(iZ(isinf(iZ)));

Z = rZ+1i*iZ;

switch solve_method
    case {CLASSIC, LS}
        h = Z\x(1:len_vandermonde) ;
    case {RKHS, ID}
        [~, h] = L_curve(Z'*Z, Z'*x(1:len_vandermonde), solve_method, 1);
    case {RKHS_h0, ID_h0}
        A = [Z;100*log(conj(r'))];
        b = [x(1:len_vandermonde); 0];
        method = solve_method(1:end-3);
        [~, h] = L_curve(A'*A, A'*b, method, 1);
    case {CLASSIC_h0, LS_h0}
        A = [Z;100*log(conj(r'))];
        b = [x(1:len_vandermonde); 0];
        %         [~, h] = L_curve(A'*A, A'*b, 'RKHS', 1) ;
        h = (A'*A)\(A'*b);
        a = 1;
    case TLS
        % if exists nan values SVD won't work
        indeterminate_form = sum(sum(isnan(Z) | isinf(Z)));
        if (indeterminate_form)
            Amp = []; alfa = []; freq = []; theta = [];
            return;
        else
            h = tls(Z,x(1:len_vandermonde));
        end
end
Amp = abs(h);
theta = atan2(imag(h), real(h)) ;


f = @(x) 0*x;
for t = 1:p
    f = @(x) f(x) + h(t).*r(t).^(x/Ts);
end

% f = @(x) (f(x));

I.Amp = Amp;
I.theta = theta;
I.freq = freq;
I.alfa = alfa;
I.r = r;
I.h = h;
I.f = f;
I.a = a;
I.c = c;

end
