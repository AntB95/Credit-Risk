%PS9 CRM BEDANIAN DE LESTABLE RICHARDI SCHOENENBERGER
clear all;
rng(420);
%Exhibits
cap_vol=[29.3 29.3 29.3 20.8 20.8 18.3 18.3 17.8 17.8 ... 
16.3 16.3 16.7 16.7 16.1 16.1 15.7 15.7 15.7 15.7]/100;
LIBOR=[4.228 2.791 3.067 3.067 3.728 3.728 4.051 4.051 4.199 4.199 ...
    4.450 4.450 4.626 4.626 4.816 4.816 4.960 4.960 5.088 5.088]/100;
beta_table=[-0.4 -0.3273 -0.2545 -0.1818 -0.1091 -0.0364 ...
       0.0364 0.1091 0.1818 0.2545 0.3273 0.4000];
beta=0.07;

%Parameters
%sufficient step as per slides
step=1/12;
delta=1/2;
Maturities=0:0.5:10;
T=2:1:10;
Reset=Maturities(2:end-1)/step+1;
M=20;

% a)
% We use formula on slide 509
vm=get_vm(beta,cap_vol,Maturities)

%b)
k=0.035;
%MC simulation under both mesure and both correlation specifications
LIBOR_Risk_Neutral_I=LiborMC(LIBOR(2:end),step,delta,Maturities,vm,beta,1,10000,0,1);
LIBOR_Risk_Neutral_II=LiborMC(LIBOR(2:end),step,delta,Maturities,vm,beta,M-1,10000,0,1);
LIBOR_Forward_I=LiborMC(LIBOR(2:end),step,delta,Maturities,vm,beta,1,10000,0,2);
LIBOR_Forward_II=LiborMC(LIBOR(2:end),step,delta,Maturities,vm,beta,M-1,10000,0,2);

MM_Risk_Neutral_I=MoneyMarket(LIBOR(1),LIBOR_Risk_Neutral_I,Reset);
MM_Risk_Neutral_II=MoneyMarket(LIBOR(1),LIBOR_Risk_Neutral_II,Reset);
MM_Forward_I=MoneyMarket(LIBOR(1),LIBOR_Forward_I,Reset);
MM_Forward_II=MoneyMarket(LIBOR(1),LIBOR_Forward_II,Reset);

BondPrice=cumprod(1./(1+delta*LIBOR));

Cap_Maturities=[3 5 7 9 11 13 15 17 19];
%Get caplet prices
for i=1:M-1
    Caplets1(i)=mean(delta*subplus(squeeze(LIBOR_Risk_Neutral_I(:,Reset(i),i))-k)./MM_Risk_Neutral_I(:,i+2),1);
    Caplets2(i)=mean(delta*subplus(squeeze(LIBOR_Risk_Neutral_II(:,Reset(i),i))-k)./MM_Risk_Neutral_II(:,i+2),1);
    Caplets3(i)=mean(delta*subplus(squeeze(LIBOR_Forward_I(:,Reset(i),i))-k).*MM_Forward_I(:,end)./MM_Forward_I(:,i+2),1).*BondPrice(end);
    Caplets4(i)=mean(delta*subplus(squeeze(LIBOR_Forward_II(:,Reset(i),i))-k).*MM_Forward_II(:,end)./MM_Forward_II(:,i+2),1)*BondPrice(end);
end

%Get cap prices
 for i=1:length(T)
     Caps1(i)=sum(Caplets1(1:Cap_Maturities(i)))*10^4;
     Caps2(i)=sum(Caplets2(1:Cap_Maturities(i)))*10^4;
     Caps3(i)=sum(Caplets3(1:Cap_Maturities(i)))*10^4;
     Caps4(i)=sum(Caplets4(1:Cap_Maturities(i)))*10^4;
 end

 %Caps1 and 2 are compute under the risk neutral probability and Caps3 and
 %4 are compute under the terminal forward mesure. Caps 1 and 3 are compute
 %using specification I and 2 and 4 with specification II.
Caps1
Caps2
Caps3
Caps4
 %c)
 %Swaptions prices
 for i=1:length(beta_table)
     vm_beta=get_vm(beta_table(i),cap_vol,Maturities);
     Swaptions1(i)=Swaption(LIBOR,step,delta,Maturities,vm_beta,beta_table(i),1,10000,0,1)*10^4;
     Swaptions2(i)=Swaption(LIBOR,step,delta,Maturities,vm_beta,beta_table(i),M-1,10000,0,1)*10^4;
     Swaptions3(i)=Swaption(LIBOR,step,delta,Maturities,vm_beta,beta_table(i),1,10000,0,2)*10^4;
     Swaptions4(i)=Swaption(LIBOR,step,delta,Maturities,vm_beta,beta_table(i),M-1,10000,0,2)*10^4;
 end
 %Swaptions1 and 2 are compute under the risk neutral probability and Swaptions3 and
 %4 are compute under the terminal forward mesure. Swaptions 1 and 3 are compute
 %using specification I and 2 and 4 with specification II.
 Swaptions1
 Swaptions2
 Swaptions3
 Swaptions4
 %d)
 %Run different correlation specification
 gammas=[0.1 1 2];
 for i=1:length(gammas)
     Swaptions1_gamma(i)=Swaption(LIBOR,step,delta,Maturities,vm,beta,M-1,10000,gammas(i),1)*10^4;
     Swaptions2_gamma(i)=Swaption(LIBOR,step,delta,Maturities,vm,beta,M-1,10000,gammas(i),2)*10^4;
 end
 %Swaptions1_gamma is compute under risk neutral mesure and Swaptions2_gamma under the terminal forward mesure.
 %using specification II.
 Swaptions1_gamma
 Swaptions2_gamma
function res=get_vm(beta,cap_vol,Tenors)
    res=sqrt(cap_vol.^2.*Tenors(2:end-1).*2.*beta./(1-exp(-2.*beta.*Tenors(2:end-1))))';
end

function res=MoneyMarket(L0,Libor,Reset)
M=20;
delta=1/2;
res=zeros(size(Libor,1),M+1);
res(:,1)=1;
res(:,2)=(1+delta*L0)*res(:,1);
for i=1:M-1
    res(:,i+2)=(1+delta*squeeze(Libor(:,Reset(i),i))).*res(:,i+1);
end
end

function res=LiborMC(L0,step,delta,Tenors,vm,beta,d,P,gamma,measure)
N=9.5/step;
M=20;
T=cumsum([zeros(M-1,1) ones(M-1,N)*step],2);
TOL=10^-2*step;
indM=bsxfun(@lt,T+TOL,Tenors(2:end-1)');
sigma=repmat(vm,1,N+1).*exp(-beta*bsxfun(@minus,Tenors(2:end-1)',T(:,:))).*indM;
sigma(sigma==0)=nan;
% We define the correlation Matrix
if d==1
    rho=1;
    L=ones(M-1,1);
else if gamma ~=0
     rho=exp(-gamma*abs(bsxfun(@minus,Tenors(2:end-1),Tenors(2:end-1)')));
     else rho=diag(ones(d,1));
    end
    L=chol(rho,'lower');
end
Z=randn(P,N,d);
H=nan(P,N+1,M-1);
H0=log(L0);
H(:,1,:)=repmat(H0,P,1,1); 
%Risk Neutral
if (measure==1)
for i=1:N
    lambdaProd=repmat(sigma(:,i),1,M-1).*rho.*repmat(sigma(:,1)',M-1,1);
    for m=1:M-1
        if indM(m,i)==1
          lambda_im=sigma(m,i)*L(m,:);
          MM=(delta*exp(H(:,i,:)))./(1+delta*exp(H(:,i,:))).*reshape(repmat(lambdaProd(m,:),P,1),P,1,[]);
          alpha=nansum(MM(:,:,1:m),3)-1/2*norm(lambda_im)^2;
          H(:,i+1,m)=H(:,i,m)+alpha*step+squeeze(Z(:,i,:))*lambda_im'*sqrt(step);
        end 
    end
end
%Forward Measure
else if (measure==2)
        for i=1:N
            lambdaProd=repmat(sigma(:,i),1,M-1).*rho.*repmat(sigma(:,1)',M-1,1);
            for m=1:M-1
                if indM(m,i)==1
                    lambda_im=sigma(m,i)*L(m,:);
                    MM=(delta*exp(H(:,i,:)))./(1+delta*exp(H(:,i,:))).*reshape(repmat(lambdaProd(m,:),P,1),P,1,[]);
                    alpha=-nansum(MM(:,:,m+1:end),3)-1/2*norm(lambda_im)^2;
                    H(:,i+1,m)=H(:,i,m)+alpha*step+squeeze(Z(:,i,:))*lambda_im'*sqrt(step);
                end 
            end
        end
    end
end
res=exp(H);
end  

function res=Swaption(L0,step,delta,Tenors,vm,beta,d,P,gamma,measure)
indT=4/step+1; 
T8=8;
Reset=Tenors(2:end-1)/step+1;
Libor=LiborMC(L0(2:end),step,delta,Tenors,vm,beta,d,P,gamma,measure);
% Bond Price at T0
BondPrice=cumprod(1./(1+delta*L0));
% We set the strike
K=(BondPrice(T8)-BondPrice(end))/(sum(BondPrice(T8+2:2:end)));
L_T8=squeeze(Libor(:,indT,T8:end));
% Bond Price at T8
BondT8=cumprod(1./(1+delta*L_T8),2);
RSwapT8=(1-BondT8(:,end))./(sum(BondT8(:,2:2:end),2));
Payoff=sum(BondT8(:,2:2:end),2).*subplus(RSwapT8-K);
MM=MoneyMarket(L0(1),Libor,Reset);
%Risk Neutral
if (measure==1)
    res=mean(Payoff./MM(:,9));
%Forward Measure
else if (measure==2)
    res=mean(Payoff.*MM(:,end)./MM(:,9)).*BondPrice(end);
    end
end
end
