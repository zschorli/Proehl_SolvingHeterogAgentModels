% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Amsterdam
% DATE May 2025
%
% DESCRIPTION
% This file sets the base random vaiables and generates the corresponding 
% orthogonal polynomials for the polynomial chaos expansion.
%__________________________________________________________________________
function setPoly(StaticParams,folder,order,case_nr)

% Load the solution of the model without aggregate risk

switch case_nr
    case {1,5,6,13,14,19,20,21}
        count = 1;
    case {2,3,4,7}
        count = case_nr;
    otherwise
        count = 7;
end
preSol = load(strcat('res_case',num2str(count),'/preSol_PFI.mat'));

%--------------------------------------------------------------------------
% Set the base distributions indexed by n
%--------------------------------------------------------------------------
Poly = struct;
for n=1:2
    if n==1 % \xi^z
        Poly.grid{n} = [1,2,3];
        Poly.pdf{n} = [StaticParams.ur(2),StaticParams.ur(1)-StaticParams.ur(2),StaticParams.er(1)];  % base distribution
    elseif n==2 % \xi^k
        Poly.grid{n} = StaticParams.kGrid;
        Poly.pdf{n} = 2.*StaticParams.p(1:2)'*preSol.pdf_cond_b;%StaticParams.p'*[preSol.pdf_cond_b;preSol.pdf_cond_g]; 
    end
    Poly.pdf{n} = Poly.pdf{n}./repmat(sum(Poly.pdf{n},2),[1,length(Poly.grid{n})]);
    Poly.cdf{n} = min(1,cumsum(Poly.pdf{n},2)); 
    Poly.cdf{n}(:,(Poly.cdf{n}==max(Poly.cdf{n}))) = 1;
end

%--------------------------------------------------------------------------
% Computing the orthogonal polynomials w.r.t. the base distributions
%--------------------------------------------------------------------------
for M=0:order
Poly.phi_N = M;                           %  max. degree of the polynomial

% generate the polynomials of the base random variables separately
for n=1:2
    Poly.phi_coeffs{n} = zeros(Poly.phi_N+1,Poly.phi_N+1);
    Poly.xi{n} = repmat(Poly.grid{n},[Poly.phi_N+1,1])...
          .^(repmat((0:Poly.phi_N)',[1,length(Poly.grid{n})]));   
    for i=1:Poly.phi_N+1
        if i==1
           Poly.phi_coeffs{n}(1,1) = 1;
        elseif i==2
            alpha = Poly.grid{n}*Poly.pdf{n}';
            if i<=Poly.phi_N+1
               Poly.phi_coeffs{n}(2,1:2) = [-alpha,1];
            end
        else        
           if i<=Poly.phi_N+1
               alpha = ((Poly.grid{n}.*(Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}')...
                       /max(eps,((Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}');
               beta = (((Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}')...
                      /max(eps,((Poly.phi_coeffs{n}(i-2,:)*Poly.xi{n}).^2)*Poly.pdf{n}');
               Poly.phi_coeffs{n}(i,:) = [0,Poly.phi_coeffs{n}(i-1,1:end-1)]...
                     -alpha*Poly.phi_coeffs{n}(i-1,:)-beta*Poly.phi_coeffs{n}(i-2,:);
           end
        end
    end
    Poly.phi_coeffs{n}=max((Poly.phi_coeffs{n}>eps),(Poly.phi_coeffs{n}<-eps)).*Poly.phi_coeffs{n};
end
Poly1 = Poly.phi_coeffs{1}*Poly.xi{1};
Poly2 = Poly.phi_coeffs{2}*Poly.xi{2};

% Form the product to obtain the joint polynomial
Poly.total = zeros(Poly.phi_N+1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)));
for i=1:Poly.phi_N+1
    for j=1:i
    for l=1:i
        if (j-1)+(l-1)==(i-1)
        Poly_mult = repmat(Poly1(j,:)',[1,length(Poly2(1,:))]).*repmat(Poly2(l,:),[length(Poly1(1,:)),1]);
        Poly.total(i,:,:,:,:) =squeeze(Poly.total(i,:,:,:,:))+Poly_mult;
        end
    end
    end
end
clear alpha beta i n phi_N xi phi_coeffs Poly1 Poly2 Poly3 Poly_mult;
Poly.total_pdf = repmat(Poly.pdf{1}',[1,length(Poly.pdf{2})]).*repmat(Poly.pdf{2},[length(Poly.pdf{1}),1]);

%--------------------------------------------------------------------------
% Computing the transition probabilities for \xi^z
%--------------------------------------------------------------------------
T = @(j,l) StaticParams.P(2*(j-1)+(1:2),2*(l-1)+(1:2))...
           ./repmat(sum((StaticParams.P(2*(j-1)+(1:2),2*(l-1)+(1:2))),2),[1,2]);
Poly.T = zeros(2,2,3,3);
T1 = T(1,1); 
Poly.T(1,1,:,:) = [Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(1,1),T1(1,2);...
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(1,1),T1(1,2);...
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(2,1),T1(2,2)];
T1 = T(1,2);
Poly.T(1,2,:,:) = [T1(1,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(1,2);...
                   T1(1,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(1,2);...
                   T1(2,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(2,2)];
T1 = T(2,1);
Poly.T(2,1,:,:) = [Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(1,1),T1(1,2);...
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(2,1),T1(2,2);...
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(2,1),T1(2,2)];
T1 = T(2,2);
Poly.T(2,2,:,:) =  [T1(1,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(1,2);...
                    T1(2,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(2,2);...
                    T1(2,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(2,2)];
clear T1 T;

save(strcat(folder,'/Poly',num2str(M),'_PFI.mat'),'-struct','Poly');
end
end