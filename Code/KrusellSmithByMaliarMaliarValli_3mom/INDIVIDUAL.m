% The program for the article "Solving the incomplete markets model with
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)

function [kprime,c]  = INDIVIDUAL(prob,ur_b,ur_g,ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id,k,km,kvar,kskew,er_b,er_g,a,epsilon,l_bar,alpha,delta,gamma,beta,mu,km_max,km_min,kvar_max,kvar_min,kskew_max,kskew_min,kprime,B,B2,B3,criter_k,k_min,k_max,update_k);
%
%__________________________________________________________________________
%
% Auxilary matrices of transition probabilities (needed for computing the 
% expectation term in the Euler equation)  
%__________________________________________________________________________

prob_bu=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for a bad agg. state 
                           % and unemployed idios. state in the next period
prob_be=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for a bad agg. state 
                           % and employed idios. state in the next period
prob_gu=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for a good agg. state 
                           % and unemployed idios. state in the next period
prob_ge=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for a good agg. state 
                           % and employed idios. state in the next period

prob_bu(:,:,:,:,1,1)=prob(1,1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_bu(:,:,:,:,1,2)=prob(2,1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_bu(:,:,:,:,2,1)=prob(3,1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_bu(:,:,:,:,2,2)=prob(4,1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);

prob_be(:,:,:,:,1,1)=prob(1,2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_be(:,:,:,:,1,2)=prob(2,2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_be(:,:,:,:,2,1)=prob(3,2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_be(:,:,:,:,2,2)=prob(4,2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);

prob_gu(:,:,:,:,1,1)=prob(1,3)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_gu(:,:,:,:,1,2)=prob(2,3)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_gu(:,:,:,:,2,1)=prob(3,3)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_gu(:,:,:,:,2,2)=prob(4,3)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);

prob_ge(:,:,:,:,1,1)=prob(1,4)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_ge(:,:,:,:,1,2)=prob(2,4)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_ge(:,:,:,:,2,1)=prob(3,4)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);
prob_ge(:,:,:,:,2,2)=prob(4,4)*ones(ngridk,ngridkm,ngridkvar,ngridkskew);

%__________________________________________________________________________
% 
% Auxilary matrices (needed for computing interest rate, wage and wealth)  
%__________________________________________________________________________

kaux=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for individual capital
kaux(:,:,:,:,1,1)=repmat(k,[1,ngridkm,ngridkvar,ngridkskew]);
kaux(:,:,:,:,1,2)=repmat(k,[1,ngridkm,ngridkvar,ngridkskew]);
kaux(:,:,:,:,2,1)=repmat(k,[1,ngridkm,ngridkvar,ngridkskew]);
kaux(:,:,:,:,2,2)=repmat(k,[1,ngridkm,ngridkvar,ngridkskew]);

kmaux=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for the mean of capital 
                                                   % distribution (km)
kmaux(:,:,:,:,1,1)=repmat(km',[ngridk,1,ngridkvar,ngridkskew]);
kmaux(:,:,:,:,1,2)=repmat(km',[ngridk,1,ngridkvar,ngridkskew]);
kmaux(:,:,:,:,2,1)=repmat(km',[ngridk,1,ngridkvar,ngridkskew]);
kmaux(:,:,:,:,2,2)=repmat(km',[ngridk,1,ngridkvar,ngridkskew]);

kmaux2=kmaux;
kmaux2(:,:,:,:,1,1)=repmat(permute(kvar,[3,2,1]),[ngridk,ngridkm,1,ngridkskew]);
kmaux2(:,:,:,:,1,2)=repmat(permute(kvar,[3,2,1]),[ngridk,ngridkm,1,ngridkskew]);
kmaux2(:,:,:,:,2,1)=repmat(permute(kvar,[3,2,1]),[ngridk,ngridkm,1,ngridkskew]);
kmaux2(:,:,:,:,2,2)=repmat(permute(kvar,[3,2,1]),[ngridk,ngridkm,1,ngridkskew]);

kmaux3=kmaux;
kmaux3(:,:,:,:,1,1)=repmat(permute(kskew,[4,3,2,1]),[ngridk,ngridkm,ngridkvar,1]);
kmaux3(:,:,:,:,1,2)=repmat(permute(kskew,[4,3,2,1]),[ngridk,ngridkm,ngridkvar,1]);
kmaux3(:,:,:,:,2,1)=repmat(permute(kskew,[4,3,2,1]),[ngridk,ngridkm,ngridkvar,1]);
kmaux3(:,:,:,:,2,2)=repmat(permute(kskew,[4,3,2,1]),[ngridk,ngridkm,ngridkvar,1]);

aglabor=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for aggregate labor 
                                              % (denoted by L in the paper)
aglabor(:,:,:,:,1,:)=er_b*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_id);
aglabor(:,:,:,:,2,:)=er_g*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag);

agshock_aux=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for the aggregate 
                                                         % shock (a)
agshock_aux(:,:,:,:,1,:)=a(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_id);
agshock_aux(:,:,:,:,2,:)=a(2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_id);

idshock_aux=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id); % for the idiosyncratic 
                                                         % shock (epsilon)
idshock_aux(:,:,:,:,:,1)=epsilon(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag);
idshock_aux(:,:,:,:,:,2)=epsilon(2)*ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag);

%__________________________________________________________________________
%
% Interest rate, wage and wealth under given k and km
%__________________________________________________________________________

ones4=ones(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_ag); % a four-dimensional matrix 
                                                  % of ones
                                                  
irateaux=alpha*(agshock_aux.*(kmaux./aglabor/l_bar).^(alpha-1)); % a 100*4*2*2-
                                         % dimensional matrix of interest rates
                                         
wageaux=(1-alpha)*(agshock_aux.*(kmaux./aglabor/l_bar).^alpha); % a 100*4*2*2-
                                         % dimensional matrix of wages
                                         
wealth=irateaux.*kaux+(wageaux.*idshock_aux)*l_bar+mu*(wageaux.*(ones4-idshock_aux))+(1-delta)*kaux-mu*(wageaux.*(1-aglabor)./aglabor).*idshock_aux; 
   % a 100*4*2*2-dimensional matrix of individual wealth

%__________________________________________________________________________
% 
% A new mean of capital distribution (km') for the given vector of the ALM 
% coefficients B (denoted by b in the paper)
%__________________________________________________________________________

% Remind that (the ALM for a bad aggregate state is ln(km')=B(1)+B(2)*ln(km) 
% and that for a good aggregate state is ln(km')=B(3)+B(4)*ln(km)

kmprime=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id);

% kmprime(:,:,1,1)=exp(B(1)*ones(ngridk,ngridkm)+B(2)*log(kmaux(:,:,1,1)));
% kmprime(:,:,1,2)=exp(B(1)*ones(ngridk,ngridkm)+B(2)*log(kmaux(:,:,1,2)));
% kmprime(:,:,2,1)=exp(B(3)*ones(ngridk,ngridkm)+B(4)*log(kmaux(:,:,2,1)));
% kmprime(:,:,2,2)=exp(B(3)*ones(ngridk,ngridkm)+B(4)*log(kmaux(:,:,2,2)));
kmprime(:,:,:,:,1,1)=exp(B(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B(2)*log(kmaux(:,:,:,:,1,1))+B(3)*log(kmaux2(:,:,:,:,1,1))+B(4)*(kmaux3(:,:,:,:,1,1)));
kmprime(:,:,:,:,1,2)=exp(B(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B(2)*log(kmaux(:,:,:,:,1,2))+B(3)*log(kmaux2(:,:,:,:,1,2))+B(4)*(kmaux3(:,:,:,:,1,2)));
kmprime(:,:,:,:,2,1)=exp(B(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B(6)*log(kmaux(:,:,:,:,2,1))+B(7)*log(kmaux2(:,:,:,:,2,1))+B(8)*(kmaux3(:,:,:,:,2,1)));
kmprime(:,:,:,:,2,2)=exp(B(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B(6)*log(kmaux(:,:,:,:,2,2))+B(7)*log(kmaux2(:,:,:,:,2,2))+B(8)*(kmaux3(:,:,:,:,2,2)));
kmprime=(kmprime>=km_min).*(kmprime<=km_max).*kmprime+(kmprime<km_min)*km_min+(kmprime>km_max)*km_max; % restricting km' to be in [km_min,km_max] range

kvarprime=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id);
kvarprime(:,:,:,:,1,1)=exp(B2(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B2(2)*log(kmaux(:,:,:,:,1,1))+B2(3)*log(kmaux2(:,:,:,:,1,1))+B2(4)*(kmaux3(:,:,:,:,1,1)));
kvarprime(:,:,:,:,1,2)=exp(B2(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B2(2)*log(kmaux(:,:,:,:,1,2))+B2(3)*log(kmaux2(:,:,:,:,1,2))+B2(4)*(kmaux3(:,:,:,:,1,2)));
kvarprime(:,:,:,:,2,1)=exp(B2(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B2(6)*log(kmaux(:,:,:,:,2,1))+B2(7)*log(kmaux2(:,:,:,:,2,1))+B2(8)*(kmaux3(:,:,:,:,2,1)));
kvarprime(:,:,:,:,2,2)=exp(B2(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B2(6)*log(kmaux(:,:,:,:,2,2))+B2(7)*log(kmaux2(:,:,:,:,2,2))+B2(8)*(kmaux3(:,:,:,:,2,2)));
kvarprime=(kvarprime>=kvar_min).*(kvarprime<=kvar_max).*kvarprime+(kvarprime<kvar_min)*kvar_min+(kvarprime>kvar_max)*kvar_max;

kskewprime=zeros(ngridk,ngridkm,ngridkvar,ngridkskew,nstates_ag,nstates_id);
kskewprime(:,:,:,:,1,1)=exp(B3(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B3(2)*log(kmaux(:,:,:,:,1,1))+B3(3)*log(kmaux2(:,:,:,:,1,1))+B3(4)*(kmaux3(:,:,:,:,1,1)));
kskewprime(:,:,:,:,1,2)=exp(B3(1)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B3(2)*log(kmaux(:,:,:,:,1,2))+B3(3)*log(kmaux2(:,:,:,:,1,2))+B3(4)*(kmaux3(:,:,:,:,1,2)));
kskewprime(:,:,:,:,2,1)=exp(B3(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B3(6)*log(kmaux(:,:,:,:,2,1))+B3(7)*log(kmaux2(:,:,:,:,2,1))+B3(8)*(kmaux3(:,:,:,:,2,1)));
kskewprime(:,:,:,:,2,2)=exp(B3(5)*ones(ngridk,ngridkm,ngridkvar,ngridkskew)+B3(6)*log(kmaux(:,:,:,:,2,2))+B3(7)*log(kmaux2(:,:,:,:,2,2))+B3(8)*(kmaux3(:,:,:,:,2,2)));
kskewprime=(kskewprime>=kskew_min).*(kskewprime<=kskew_max).*kskewprime+(kskewprime<kskew_min)*kskew_min+(kskewprime>kskew_max)*kskew_max;

%__________________________________________________________________________
%
% Future interest rate and wage conditional on state
%__________________________________________________________________________

irate_b=alpha*a(1)*((kmprime./(er_b*ones4*l_bar)).^(alpha-1));  
                                      % under a bad future aggregate state
irate_g=alpha*a(2)*((kmprime./(er_g*ones4*l_bar)).^(alpha-1));  
                                      % under a good future aggregate state
                                      
wage_b=(1-alpha)*a(1)*((kmprime./(er_b*ones4*l_bar)).^(alpha)); 
                                      % under a bad future aggregate state
wage_g=(1-alpha)*a(2)*((kmprime./(er_g*ones4*l_bar)).^(alpha)); 
                                      % under a good future aggregate state

%__________________________________________________________________________
% 
% SOLVING THE INDIVIDUAL PROBLEM
%__________________________________________________________________________   

dif_k=1;
   while dif_k>criter_k
   
   % Bad aggregate state and unemployed idiosyncratic state 
   
     k2prime_bu=interpn(k,km,kvar,kskew,kprime(:,:,:,:,1,1),kprime,kmprime,kvarprime,kskewprime,'linear'); 
        % finding the individual policy function k''=k(k',km') by interpolating
        % the previously found policy function k'=k(k,km) in new points (k',km')
     cprime_bu=irate_b.*kprime+mu*(wage_b.*ones4)+(1-delta)*kprime-k2prime_bu; 
                                                  % future consumption (c')
     cprime_bu=(cprime_bu>0).*cprime_bu+(cprime_bu<=0)*10^-10; 
                      % if c'<=0, then set it to a very low positive number 
     muprime_bu=cprime_bu.^(-gamma); % marginal utility of future consumption
   
   % Bad aggregate state and employed idiosyncratic state
   
     k2prime_be=interpn(k,km,kvar,kskew,kprime(:,:,:,:,1,2),kprime,kmprime,kvarprime,kskewprime,'linear');
     cprime_be=irate_b.*kprime+wage_b.*(epsilon(2)*l_bar*ones4)+(1-delta)*kprime-mu*(wage_b.*((ur_b./(1-ur_b))*ones4))-k2prime_be;
     cprime_be=(cprime_be>0).*cprime_be+(cprime_be<=0)*10^-10;
     muprime_be=cprime_be.^(-gamma);
 
   % Good aggregate state and unemployed idiosyncratic state
   
     k2prime_gu=interpn(k,km,kvar,kskew,kprime(:,:,:,:,2,1),kprime,kmprime,kvarprime,kskewprime,'linear');
     cprime_gu=irate_g.*kprime+mu*(wage_g.*ones4)+(1-delta)*kprime-k2prime_gu;
     cprime_gu=(cprime_gu>0).*cprime_gu+(cprime_gu<=0)*10^-10;
     muprime_gu=cprime_gu.^(-gamma);
   
   % Good aggregate state and employed idiosyncratic state
   
     k2prime_ge=interpn(k,km,kvar,kskew,kprime(:,:,:,:,2,2),kprime,kmprime,kvarprime,kskewprime,'linear');
     cprime_ge=irate_g.*kprime+wage_g.*(epsilon(2)*l_bar*ones4)+(1-delta)*kprime-mu*(wage_g.*((ur_g./(1-ur_g))*ones4))-k2prime_ge;
     cprime_ge=(cprime_ge>0).*cprime_ge+(cprime_ge<=0)*10^-10;
     muprime_ge=cprime_ge.^(-gamma);
   
   % Expectation term in the Euler equation (condition (5) in the paper) 
   
    expec=(muprime_bu.*((1-delta)*ones4+irate_b)).*prob_bu+...
    (muprime_be.*((1-delta)*ones4+irate_b)).*prob_be+...
    (muprime_gu.*((1-delta)*ones4+irate_g)).*prob_gu+...
    (muprime_ge.*((1-delta)*ones4+irate_g)).*prob_ge;
   
   cn=(beta*expec).^(-1/gamma); % current consumption found from the Euler 
                    %equation (5) if the borrowing constraint is not binding
                    
   kprimen=wealth-cn; % new capital function (i.e., condition (5) in the paper)
    
   kprimen=(kprimen>=k_min).*(kprimen<=k_max).*kprimen+(kprimen<k_min)*k_min+(kprimen>k_max)*k_max; 
                    % restricting kprimen to be in [km_min,km_max] range 
   
  
   dif_k=max(max(max(max(max(max(abs(kprimen-kprime),[],6),[],5),[],4),[],3),[],2),[],1); % difference between the 
                                       % new and previous capital functions
   
   kprime=update_k*kprimen+(1-update_k)*kprime; % updating kprimen (condition 
                                       % (6) in the paper)
      
   end

% Consumption function

c=wealth-kprime; % follows from the budget constraint
