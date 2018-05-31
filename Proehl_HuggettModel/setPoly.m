% This is code for the working paper "Approximating Equilibria with Ex-Post 
% Heterogeneity and Aggregate Risk" by Elisabeth Pröhl
%
% AUTHOR Elisabeth Pröhl, University of Geneva and Swiss Finance Institute
% DATE May 2018
%
% DESCRIPTION
% This file sets the base random vaiables and generates the corresponding 
% orthogonal polynomials for the polynomial chaos expansion.
%__________________________________________________________________________
function setPoly(StaticParams,folder,order)

% Use either the Proximal point algorithm (method = 0) or policy function
% iteration (method = 1)
for method = 1:-1:0
StaticParams.methodEqPFI = method;   

% Load the solution of the model without aggregate risk
if StaticParams.methodEqPFI == 0
	preSol = load(strcat(folder,'/preSol_PPA.mat'));
else
	preSol = load(strcat(folder,'/preSol_PFI.mat'));
end

%--------------------------------------------------------------------------
% Set the base distributions indexed by n
%--------------------------------------------------------------------------
Poly = struct;
for n=1:3
    switch n
        case 1 % \xi^z
        Poly.grid{n} = [1,2,3];
        Poly.pdf{n} = [StaticParams.ur(2),StaticParams.ur(1)-StaticParams.ur(2),StaticParams.er(1)];  % base distribution
        case 2 % \xi^a_1
        Poly.grid{n} = StaticParams.kGrid;
        Poly.pdf{n} = StaticParams.p([1,3])'*[preSol.pdf_cond_b(1,:);preSol.pdf_cond_g(1,:)]./sum(StaticParams.p([1,3]));
        case 3 % \xi^a_2
        Poly.grid{n} = StaticParams.kGrid;
        Poly.pdf{n} = StaticParams.p([2,4])'*[preSol.pdf_cond_b(2,:);preSol.pdf_cond_g(2,:)]./sum(StaticParams.p([2,4]));
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
for n=1:3
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
               alpha = ((Poly.grid{n}.*(Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}')/max(eps,((Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}');
               beta = (((Poly.phi_coeffs{n}(i-1,:)*Poly.xi{n}).^2)*Poly.pdf{n}')/max(eps,((Poly.phi_coeffs{n}(i-2,:)*Poly.xi{n}).^2)*Poly.pdf{n}');
               Poly.phi_coeffs{n}(i,:) = [0,Poly.phi_coeffs{n}(i-1,1:end-1)]-alpha*Poly.phi_coeffs{n}(i-1,:)-beta*Poly.phi_coeffs{n}(i-2,:);
           end
        end
    end
    Poly.phi_coeffs{n}=max((Poly.phi_coeffs{n}>eps),(Poly.phi_coeffs{n}<-eps)).*Poly.phi_coeffs{n};
end
Poly1 = Poly.phi_coeffs{1}*Poly.xi{1};
Poly2 = Poly.phi_coeffs{2}*Poly.xi{2};
Poly3 = Poly.phi_coeffs{3}*Poly.xi{3};

% Form the product to obtain the joint polynomial
Poly.total = zeros(Poly.phi_N+1,length(Poly.xi{1}(1,:)),length(Poly.xi{2}(1,:)),length(Poly.xi{3}(1,:)));
for i=1:Poly.phi_N+1
    for j=1:i
    for l=1:i
    for m=1:i
        if (j-1)+(l-1)+(m-1)==(i-1)
        Poly_mult = repmat(Poly1(j,:)',[1,length(Poly2(1,:)),length(Poly3(1,:))])...
                    .*repmat(Poly2(l,:),[length(Poly1(1,:)),1,length(Poly3(1,:))])...
                    .*repmat(permute(Poly3(m,:),[3,1,2]),[length(Poly1(1,:)),length(Poly2(1,:)),1]);
        Poly.total(i,:,:,:,:) =squeeze(Poly.total(i,:,:,:,:))+Poly_mult;
        end
    end
    end
    end
end
clear alpha beta i n phi_N xi phi_coeffs Poly1 Poly2 Poly3 Poly4 Poly_mult a;
Poly.total_pdf = repmat(Poly.pdf{1}',[1,length(Poly.pdf{2}),length(Poly.pdf{3})])...
                 .*repmat(Poly.pdf{2},[length(Poly.pdf{1}),1,length(Poly.pdf{3})])...
                 .*repmat(permute(Poly.pdf{3},[3,1,2]),[length(Poly.pdf{1}),length(Poly.pdf{2}),1]);
             
%--------------------------------------------------------------------------
% Computing the transition probabilities for \xi^z
%--------------------------------------------------------------------------
T = @(j,l) StaticParams.P(2*(j-1)+(1:2),2*(l-1)+(1:2))...
           ./repmat(sum((StaticParams.P(2*(j-1)+(1:2),2*(l-1)+(1:2))... 
                    ),2),[1,2]);
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
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(1,1),T1(1,2);...
                   Poly.pdf{1}(1:2)/sum(Poly.pdf{1}(1:2)).*T1(1,1),T1(1,2)];
T1 = T(2,2);
Poly.T(2,2,:,:) =  [T1(1,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(1,2);...
                    T1(1,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(1,2);...
                    T1(2,1),Poly.pdf{1}(2:3)/sum(Poly.pdf{1}(2:3)).*T1(2,2)];
clear T1 T;

if StaticParams.methodEqPFI == 0
	save(strcat(folder,'/Poly',num2str(M),'_PPA.mat'),'-struct','Poly');
else
	save(strcat(folder,'/Poly',num2str(M),'_PFI.mat'),'-struct','Poly');
end
end
end
end