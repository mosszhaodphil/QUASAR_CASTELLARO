function [aif,r_true,rm_true,DM_true,DM,sigma]=symASL_QUASAR(t,f,r_type,s,delay,snr,nr)
step=0.3;
delay_second = delay;
delay=floor(delay./.3);
lambda = 0.9;
T1t = 1.3;
T1b = 1.6;
A = 1;
alpha = 2*pi/360*35;
deltaTI = t(2)- t(1);
Dt = 0.25 ;
tau = 0.64;
lookloker_param = (1/T1t-log(cos(alpha).^(delay)));
p = 0.1;
s=exp(s);
m = exp(-t.*lookloker_param);
nTI = length(t);

f_snr = .01;
r_snr=exp(-t*f_snr/lambda);
%Residue function
switch r_type
    case 1 % fisiological exponential
        r=exp(-t*f/lambda);
        
    case 2 % gamma-variate
        r=t.*exp(-t*sqrt(f/lambda));
    case 3 % low dispersed exponential
        theta = 0.5;
        r=1./(f/lambda*theta -1) .* exp(-t./theta) - 1./(f/lambda*theta -1).*exp(-f/lambda.*t);
    case 4 % medium dispersed exponential
        theta = 0.05;
        r=1./(f/lambda*theta -1) .* exp(-t./theta) - 1./(f/lambda*theta -1).*exp(-f/lambda.*t);
    case 5 % high dispersed exponential
        theta = 0.1;
        r=1./(f/lambda*theta -1) .* exp(-t./theta) - 1./(f/lambda*theta -1).*exp(-f/lambda.*t);
end

%Arterial Input function
aif = zeros(length(t),1);
for k=1:length(t)
    if and(t(k) >= Dt , t(k) < (Dt + tau))
        aif(k) = (1 - gammainc(1+s*p,s*(t(k)-Dt)));
    elseif t(k) >= (Dt + tau)
        aif(k) = (gammainc(1+s*p,s*(t(k)-Dt - tau)) - gammainc(1+s*p,s*(t(k)-Dt)));
    end
end
aif = aif .* exp(-t/T1b)';

r_true=r;
rm_true = r.*m;
DM_snr = A*f_snr.*(deltaTI).*conv(aif,m.*r_snr);

switch r_type
    case 1 %  exponential
        DM_noiseless = A*f.*(deltaTI).*conv(aif,m.*r);
        DM_noiseless = DM_noiseless(1:floor(end/2+1));
    case 2 % gamma-variate
        DM_noiseless = A*f.*(deltaTI).*conv(aif,m.*r);
        DM_noiseless = DM_noiseless(1:floor(end/2+1));
    case 3 % low dispersed exponential
        DM_noiseless = A*f.*(deltaTI).*conv(aif,m.*r);
        DM_noiseless = DM_noiseless(1:floor(end/2+1));
    case 4 % medium dispersed exponential
        DM_noiseless = A*f.*(deltaTI).*conv(aif,m.*r);
        DM_noiseless = DM_noiseless(1:floor(end/2+1));
    case 5 % high dispersed exponential
        DM_noiseless = A*f.*(deltaTI).*conv(aif,m.*r);
        DM_noiseless = DM_noiseless(1:floor(end/2+1));
end



if delay == 0
    DM_true = DM_noiseless;
else
    
    if step == 0.3
        DM_true((delay+1):nTI) = DM_noiseless (1:(nTI-delay));
        DM_true(1:delay)=0;
        DM_true = DM_true.*exp(-delay_second/T1b);
    elseif step == 0.06
        DM_true((delay*5+1):nTI) = DM_noiseless (1:(nTI-delay*5));
        DM_true(1:delay*5)=0;
        DM_true = DM_true.*exp(-delay_second/T1b);
    else
        t(2) - t(1)
    end
end

%Noise level estimation
sigma = max(DM_snr)/snr;

%Montecarlo simulations
DM = zeros(nTI,nr);
for k=1:nr
    noise = sigma.*randn(1,nTI);
    DM(:,k) = DM_true + noise;
end


