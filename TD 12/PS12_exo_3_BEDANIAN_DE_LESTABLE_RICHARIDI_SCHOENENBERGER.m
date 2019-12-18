%% Credit Risk - Problem Set 12, Exercise 3

%% Parameters
br = 0.077; 
betar = -0.96; 
sigmar = 0.014;
bg = 0.012; 

betag = -0.2; 
sigmag = 0.14;
t = 0; 
T = 1;

%A)

N = 1000; 
M = 1000;
dt = (T-t)/N;
dWr = sqrt(dt)*randn(N,M);
dWg = sqrt(dt)*randn(N,M);

r0 = 0.01;
r = nan(N,M);
r(1,:) = r0;
for t = 2:1:N
    r(t,:) = r(t-1,:)+(br+betar*r(t-1,:))*dt+sigmar.*dWr(t,:);
end



gamma0 = 0.05;
gamma = nan(N,M);
gamma(1,:) = gamma0;
for t=2:1:N
    gamma(t,:) = gamma(t-1,:)+(bg+betag*gamma(t-1,:))*dt+sigmag*sqrt(gamma(t-1,:)).*dWg(t,:);
end

delta = 0.5;

% calculating integral of R over the different paths
integralR = nan(N,M);
for tao = 1:1:N
    integralR(tao,:) = (sum(r(1:tao,:),1)*dt+sum(gamma(1:tao,:),1)*dt);
end
I = mean(exp(-integralR(N,:)));
% calculating second integral
    integral2 = sum(gamma(1:N,:).*(exp(-integralR(1:N,:))),1)*dt;

J = mean(integral2);
p1 = I + (1-delta)*J;
fprintf('t = 0 price p_1(t,T) of a defaultable bond with maturity T = 1 under Recovery of Face Value (RF) is %f \n',p1);
% result = 0.941969
%B)

t = 0; T_cds = 5;
N = 1000; M = 1000;
dt = (T_cds-t)/N;
dWr = sqrt(dt)*randn(N,1);
dWg = sqrt(dt)*randn(N,1);

r0 = 0.05;
r = nan(N,M);
r(1,:) = r0;
for t = 2:1:N
    r(t,:) = r(t-1,:)+(br+betar*r(t-1,:))*dt+sigmar.*dWr(t,:);
end

%C)
gamma0 = 0.05;
gamma = nan(N,M);
gamma(1,:) = gamma0;
for t = 2:1:N
    gamma(t,:) = gamma(t-1,:)+(bg+betag*gamma(t-1,:))*dt+sigmag*sqrt(gamma(t-1,:)).*dWg(t,:);
end

%integral of R over the different paths
integralRspread = nan(N,M);
for tau = 1:1:N
    integralRspread(tau,:) = (sum(r(1:tau,:),1)*dt+sum(gamma(1:tau,:),1)*dt);
end
I=mean(exp(-integralR(N,:)));
% second integral
integral2spread=sum(gamma(1:N,:).*(exp(-integralR(1:N,:))),1)*dt;
J=mean(integral2);

premiumTimes = 0.5:0.5:5;
premiumIndices = premiumTimes/dt;

Numerator=delta*J;
Denominator = 0.5*sum(mean(exp(-integralRspread(premiumIndices,:)),2));
x0 = Numerator/Denominator;
fprintf('\n t = 0 fair swap spread x0 of CDS with semi-annual premium payments and maturity T = 5 years is %f \n',x0)
% result = 0.006280