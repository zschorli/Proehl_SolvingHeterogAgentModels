% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE September 2026
%
% DESCRIPTION
% This file implements the policy function iteration for the model with
% aggregate shocks.
%__________________________________________________________________________
function [c,h,a_prime,b,wealth,r,q_pr,diffs] ...
         = OptPol_AggShock_PFI(c,h,a_prime,wealth,r,q_pr,basic_grid,StaticParams,Poly,grid_N,distrGrid,idx,SymFuns)

diffs = 1e8.*ones(1,3);
iter=0; count = 0;
c_pr = zeros([size(a_prime,1),4,4,size(a_prime,3)]);
h_pr = zeros([size(a_prime,1),4,4,size(a_prime,3)]);
r_pr = zeros([size(a_prime,1),4,4]);
q = repmat(distrGrid(idx(:,1),1),[1,4,size(a_prime,3)]);
qmax = distrGrid(end,1)+0.*(distrGrid(end,1)-distrGrid(1,1));
qmin = distrGrid(1,1)-0.*(distrGrid(end,1)-distrGrid(1,1));
aux = c;
q_pr_new = q_pr;
options = optimoptions('fsolve','Display','off','FunctionTolerance',StaticParams.criter_a,...
           'SpecifyObjectiveGradient',true,'MaxIterations',150);
options2 = optimoptions('fsolve','Display','off','FunctionTolerance',StaticParams.criter_a,...
           'SpecifyObjectiveGradient',true);
sol_start = [r,q_pr];
update_agg = 0.7;
update = 0.1;
r_large = repmat(r(:,[1,1,2,2]),[1,1,size(c,3)]);
q_pr_large = repmat(q_pr(:,[1,1,2,2]),[1,1,size(c,3)]);
uq = max(eps,(1+r_large).*q-q_pr_large);
diff_r = 0.01;
diff_q_pr = 0.01;

innerAggr = @(c,h) c.^((StaticParams.eps-1)/StaticParams.eps)...
               +StaticParams.eta^(1/StaticParams.eps).*...
               h.^((StaticParams.eps-1)/StaticParams.eps);
Aggr = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1));
Aggr_derC = @(c,h) (innerAggr(c,h)).^(StaticParams.eps/(StaticParams.eps-1)-1)...
               .*c.^(-1/StaticParams.eps);
MargUt = @(g) g.^(-StaticParams.gamma);
MargUt_inv = @(g) g.^(-1./StaticParams.gamma);

F_uq = @(r,q,q_pr) max(eps,(1+r).*q-q_pr);
if 1/StaticParams.eps==StaticParams.gamma
    OptC_mult = @(r,q,q_pr) ones(size(q)); 
else
    OptC_mult = @(r,q,q_pr) (1+StaticParams.eta.*max(eps,F_uq(r,q,q_pr)).^(1-StaticParams.eps))...
                .^((1-StaticParams.eps*StaticParams.gamma)./(StaticParams.eps-1));
end
a_large = repmat(permute(StaticParams.aGrid_pol,[3,1,2]),[size(c,1),4,1]);
h_bound = 1/(1-StaticParams.theta).*a_large./q;
F_h = @(r,q,q_pr,c,a) min((StaticParams.eta.*(F_uq(r,q,q_pr)).^(-StaticParams.eps).*c)...
      ,(1/(1-StaticParams.theta).*a./q));
diff_grid=StaticParams.aGrid_pol(2:end)-StaticParams.aGrid_pol(1:end-1);
diff_grid=repmat(permute(diff_grid,[3,1,2]),[size(c,1),4,1]);
bool_nosolve = zeros(size(a_prime,1),1);

while (max(diffs(max(iter,1),:)) > StaticParams.criter_a) || iter<5
   parfor m=1:size(a_prime,1)
       a_prime_large = interp1(StaticParams.aGrid_pol,squeeze(a_prime(m,:,:))',StaticParams.aGrid)';   
       [~,weights_new] = calcAggregates(a_prime_large,squeeze(basic_grid(m,:,:,:,:)),StaticParams,Poly);
       c_prime = @(z_idx,z,a,w)...
		   interp1(StaticParams.aGrid_pol,...
           getCApprox(grid_N,distrGrid,squeeze(c(:,2*z(1)+z(2)+1,:)),...
           idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.a_min,a),'linear','extrap');
	   h_prime = @(z_idx,z,a,w)...
		   interp1(StaticParams.aGrid_pol,...
           getCApprox(grid_N,distrGrid,squeeze(h(:,2*z(1)+z(2)+1,:)),...
           idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'),...
           max(StaticParams.a_min,a),'linear','extrap');

	   r_prime = @(z_idx,z,w)...
		   max(-1+eps,getCApprox(grid_N,distrGrid,squeeze(r(:,z(1)+1)),...
           idx,squeeze(w((z_idx>2)+1,z(1)+1,:))'));
       for z_idx=1:4
		   c_pr(m,z_idx,:,:) = [c_prime(z_idx,[0,0],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[0,1],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,0],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                c_prime(z_idx,[1,1],squeeze(a_prime(m,z_idx,:))',weights_new)] ;
		   h_pr(m,z_idx,:,:) = [h_prime(z_idx,[0,0],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                h_prime(z_idx,[0,1],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                h_prime(z_idx,[1,0],squeeze(a_prime(m,z_idx,:))',weights_new);...
                                h_prime(z_idx,[1,1],squeeze(a_prime(m,z_idx,:))',weights_new)] ;
		   r_pr(m,z_idx,:) = [r_prime(z_idx,[0,0],weights_new);...
                              r_prime(z_idx,[0,1],weights_new);...
                              r_prime(z_idx,[1,0],weights_new);...
                              r_prime(z_idx,[1,1],weights_new)] ;
       end
   end
   for i=1:4 
       MargUt_pr = StaticParams.beta...
				.*sum(repmat(StaticParams.P(i,:),[size(a_prime,1),1,length(StaticParams.aGrid_pol)])...
                .*repmat(squeeze(1+r_pr(:,i,:)),[1,1,size(a_prime,3)])...
                .*squeeze(MargUt(Aggr(max(eps,c_pr(:,i,:,:)),...
                              max(eps,h_pr(:,i,:,:))))...
                .*Aggr_derC(max(eps,c_pr(:,i,:,:)),...
                            max(eps,h_pr(:,i,:,:)))),2);
       aux(:,i,:) = MargUt_inv(MargUt_pr);
       % c=aux.*MargUt_inv(1./OptC_mult)
   end
   aux = cummin(aux,3,"reverse");
   if count==0
       parfor m=1:size(a_prime,1)  
           [sol,~,flag] = fsolve(@(sol) solveEquCond(sol,squeeze(q(m,:,:)),squeeze(r_large(m,:,:)),squeeze(a_large(m,:,:)),squeeze(aux(m,:,:)),...
                        squeeze(h_bound(m,:,:)),squeeze(basic_grid(m,:,:,:,:)),qmin,qmax,StaticParams,Poly,SymFuns),...
                        sol_start(m,:)',options);
           if flag>0
               sol_start(m,:) = sol';
               q_pr_new(m,:) = max(qmin,min(qmax,real(sol(3:4)')));
               r_new(m,:) = max(-1+eps,sol(1:2)');
           elseif flag==0
               display([m,flag]);
               sol_start(m,:) = sol';
           else
               display([m,flag]);
               bool_nosolve(m) = 1;
           end
       end
       a = find(bool_nosolve);
       for i=1:sum(bool_nosolve)
           aux_idx = find(idx(:,1)==idx(a(i),1));
           aux_idx = find(aux_idx~=a(i));
           aux_idx = find(abs(a(i)-aux_idx)==min(abs(a(i)-aux_idx)));
           sol_start(a(i),:) = 0.5.*mean(sol_start(aux_idx,:),1);
       end
       diff_r = max(max(abs(r-r_new)));
       diff_q_pr = max(max(abs(q_pr-q_pr_new)));
       r =  (1-update_agg)*r+update_agg*r_new;
       q_pr = (1-update_agg)*q_pr+update_agg*q_pr_new;
       r_large = repmat(r(:,[1,1,2,2]),[1,1,size(c,3)]);
       q_pr_large = repmat(q_pr(:,[1,1,2,2]),[1,1,size(c,3)]);
       uq = max(eps,(1+r_large).*q-q_pr_large);
       wealth = repmat(permute(StaticParams.wage,[3,1,2]),[size(c,1),1,1])...
                +(1+r_large).*a_large;
   end

   % explicit solution of Euler equation for h unconstraint, the same holds
   % for h constraint when 1/eta=equal gamma which is our calibration
   c_new = aux.*MargUt_inv(1./OptC_mult(r_large,q,q_pr_large));
   c_new = max(eps,c_new);
   diff_c = mean(mean(mean(abs(c-c_new))));

   c = (1-update)*c +update*c_new;
   aux_h = F_h(r_large,q,q_pr_large,c,a_large);
   aux_a_prime = wealth-c-F_uq(r_large,q,q_pr_large).*aux_h;
   if count == 0
       aux_idx1 = find(aux_a_prime<=0);
       aux_idx = aux_idx1;
       c_max = [];
   else
       aux_idx1 = find(aux_a_prime<=0);
       aux_idx1 = setdiff(aux_idx1,aux_idx);
       aux_idx = [aux_idx;aux_idx1];
   end
   aux_c_max = zeros(size(aux_idx1));
   parfor i=1:length(aux_idx1)
       aux_c_max(i) = fsolve(@(sol) solveMaxCons(sol,wealth(aux_idx1(i)),uq(aux_idx1(i)),h_bound(aux_idx1(i)),SymFuns)...
                        ,0.9.*wealth(aux_idx1(i)),options2);
   end
   c_max = [c_max;aux_c_max];
   c(aux_idx) = min(c(aux_idx),c_max);
   h = F_h(r_large,q,q_pr_large,c,a_large);
   a_prime = max(0,wealth-c-F_uq(r_large,q,q_pr_large).*h);
   b = a_prime-q_pr_large.*h;

   iter = iter+1;
   diffs(iter,:) = [diff_c,diff_r,diff_q_pr];   
      
   if (diffs(iter,1)>StaticParams.criter_a) && (count<9) 
       count = count+1;
   else
       count =0;
   end
   bool_nosolve = zeros(size(a_prime,1),1);
end
display(iter); 
poolobj = gcp('nocreate');
delete(poolobj);
end