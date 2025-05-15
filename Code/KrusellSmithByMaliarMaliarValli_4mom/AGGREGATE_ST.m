% The program for the article "Solving the incomplete markets model withlinear
% aggregate uncertainty using the Krusell-Smith algorithm" from the special 
% JEDC issue edited by Den Haan, Judd and Juillard (2008)  
%
% Written by Lilia Maliar, Serguei Maliar and Fernando Valli (2008)

function [kmts,kvarts,kskewts,kkurtts,kcross_out]  = AGGREGATE_ST(T,idshock,agshock,km_max,km_min,kvar_max,kvar_min,kskew_max,kskew_min,kkurt_max,kkurt_min,kprime,km,kvar,kskew,kkurt,k,epsilon2,k_min,k_max,kcross,a2);

kmts=zeros(T,1);         % a time series of the mean of capital distribution 
kvarts=zeros(T,1);
kskewts=zeros(T,1);
kkurtts=zeros(T,1);
kmts_b=zeros(T,1);         
kvarts_b=zeros(T,1);
kskewts_b=zeros(T,1);
kkurtts_b=zeros(T,1);
kcross_out = kcross;

for t=1:T
  
   kmts(t)=mean(kcross); % find the t-th observation of kmts by computing 
                         % the mean of the t-th period cross-sectional 
                         % distribution of capital
   kmts_b(t)=kmts(t)*(kmts(t)>=km_min)*(kmts(t)<=km_max)+km_min*(kmts(t)<km_min)+km_max*(kmts(t)>km_max); % restrict kmts to be within [km_min, km_max]
   
   kvarts(t)=std(kcross);
   kvarts_b(t)=kvarts(t)*(kvarts(t)>=kvar_min)*(kvarts(t)<=kvar_max)+kvar_min*(kvarts(t)<kvar_min)+kvar_max*(kvarts(t)>kvar_max);
   kskewts(t)=(skewness(kcross));
   kskewts_b(t)=kskewts(t)*(kskewts(t)>=kskew_min)*(kskewts(t)<=kskew_max)+kskew_min*(kskewts(t)<kskew_min)+kskew_max*(kskewts(t)>kskew_max);
   kkurtts(t)=(kurtosis(kcross));
   kkurtts_b(t)=kkurtts(t)*(kkurtts(t)>=kkurt_min)*(kkurtts(t)<=kkurt_max)+kkurt_min*(kkurtts(t)<kkurt_min)+kkurt_max*(kkurtts(t)>kkurt_max);
   
   % To find kmts(t+1), we should compute a new cross-sectional distribution 
   % at t+1. For this purpose, we first find kprime by interpolation for 
   % realized kmts(t) and agshock(t) (kprimet below) and then use it to 
   % compute new kcross by interpolation (kcrossn below) given the previous 
   % kcross and the realized idshock(t)
   
   kprimet4=interpn(k,km,kvar,kskew,kkurt,a2,epsilon2,kprime,k, kmts_b(t),kvarts_b(t),kskewts_b(t),kkurtts_b(t),agshock(t),epsilon2,'spline');
      % a four-dimensional capital function at time t is obtained by fixing
      % known kmts(t) and agshock (t)
      
   kprimet=squeeze(kprimet4); % the size of kprimet4 is ngridk*1*1*nstates_id; 
                              % in kprimet, all singleton dimensions (i.e.,
                              % those with only one column per page) are removed
                              
   kcrossn=interpn(k,epsilon2,kprimet,kcross,idshock(t,:),'spline'); 
                              % given kcross and idiosyncratic shocks we
                              % compute kcrossn
                              
   kcrossn=kcrossn.*(kcrossn>=k_min).*(kcrossn<=k_max)+k_min*(kcrossn<k_min)+k_max*(kcrossn>k_max); % restrict kcross to be within [k_min, k_max]
                              
   kcross=kcrossn;
   if kmts(t)==kmts_b(t) && kvarts(t)==kvarts_b(t) && kskewts(t)==kskewts_b(t) && kkurtts(t)==kkurtts_b(t)
      kcross_out = kcrossn;
   end
   
end
end