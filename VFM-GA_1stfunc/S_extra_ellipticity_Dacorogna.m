clear;
clc;

G0 = 59.4

B = 121

Jmin = 0.1646

C1 = 0.8433

K10 =  -0.12

dK1 = 0.37667

X1p = 5.2

X2p = 0.604

C0 = 0.44667

p = 3.6

q = 4.8

C2 = 10.3333

C3 = 0.5

r = 1.0


  
% Initialize for a given K3
delta = eye(3,3);
%N = [1/sqrt(2) 0 0; 0 0 0; 0 0 -1/sqrt(2)]; % K3=0
N = [1/sqrt(6) 0 0; 0 1/sqrt(6) 0; 0 0 -2/sqrt(6)]; % K3 = -1
%N = [2/sqrt(6) 0 0; 0 -1/sqrt(6) 0; 0 0 -1/sqrt(6)]; % K3 = 1



K3 = 3*sqrt(6)*det(N);
Y = 3*sqrt(6)*N*N - sqrt(6)*eye(3,3) - 3*K3*N;

% Initialize K1 and K2
 K1v = -1.4:0.01:0.2;
 K2v = 0.01:0.01:1.2;

[K1, K2] = meshgrid(K1v,K2v);
stability = zeros(size(K1,1),size(K1,2));

%%
% Loop over pairs of K1 and K2
for ii = 1:size(K1,1)
    for jj = 1:size(K1,2)
        
        % Extract the strain tensor, Jacobian, and B-tensor
        E = (1/3)*K1(ii,jj)*delta + K2(ii,jj)*N;
        J = exp(K1(ii,jj));
        [eigvec,eigval] = eig(E);
        principal_stretches = zeros(3,3);
        for mm=1:3
            eigval(mm,mm) = exp(2*eigval(mm,mm));
            principal_stretches(mm,mm) = sqrt(eigval(mm,mm));
        end
        Bstretch = eigvec*eigval*eigvec';
        Fdefgrad = eigvec*principal_stretches*eigvec';
        Finv = inv(Fdefgrad);
        
        % Calculate X and its derivatives
        X = (1/2)*(X1p + X2p)*K1(ii,jj) + (dK1/2)*(X1p - X2p)*log((cosh((K1(ii,jj) - K10)/dK1))/(cosh((K10)/dK1))) + 1;
        dXdK1 = (1/2)*(X1p + X2p) + (1/2)*(X1p - X2p)*tanh((K1(ii,jj) - K10)/dK1);
        d2XdK12 = (1/2)*((X1p - X2p)/dK1)*(sech((K1(ii,jj) - K10)/dK1))^2;
        
        % Calculate L and its derivatives
        L = C0*(K2(ii,jj))^p + C1*(1+K3)*(K2(ii,jj))^q;
        dLdK2 = p*C0*(K2(ii,jj))^(p-1) + q*C1*(1+K3)*(K2(ii,jj))^(q-1);
        dLdK3 = C1*(K2(ii,jj))^q;
        d2LdK22 = p*(p-1)*C0*(K2(ii,jj))^(p-2) + q*(q-1)*C1*(1+K3)*K2(ii,jj)^(q-2);
        d2LdK2dK3 = q*C1*(K2(ii,jj))^(q-1);
        
        % Calculate f and its derivatives
        dfdK1 = (exp(C2*K1(ii,jj))-1)/C2 + C3*exp(K1(ii,jj))*(exp(-2*K1(ii,jj)) - ((1-Jmin)^2)/((exp(K1(ii,jj))-Jmin)^2));
        d2fdK12 = exp(C2*K1(ii,jj)) + C3*exp(K1(ii,jj))*(exp(-2*K1(ii,jj)) - ((1-Jmin)^2)/((exp(K1(ii,jj))-Jmin)^2)) + C3*exp(K1(ii,jj))*(-2*exp(-2*K1(ii,jj)) + 2*exp(K1(ii,jj))*((1-Jmin)^2)/((exp(K1(ii,jj))-Jmin)^3));
        
        % Calculate first derivatives of psi
        dpsidK1 = G0*(dXdK1*K2(ii,jj)^2) + B*dfdK1;
        dpsidK2 = G0*(2*X*K2(ii,jj) + dLdK2);
        dpsidK3 = G0*dLdK3;
        
        % Calculate the Kirchhoff stress
        TK = dpsidK1*delta + dpsidK2*N + dpsidK3*(1/K2(ii,jj))*Y;
        
        % Calculate second derivatives of psi
        d2psidK12 = G0*d2XdK12*K2(ii,jj)^2 + B*d2fdK12;
        d2psidK22 = G0*(2*X + d2LdK22);
        d2psidK32 = 0;
        d2psidK1dK2 = G0*2*dXdK1*K2(ii,jj);
        d2psidK1dK3 = 0;
        d2psidK2dK3 = G0*d2LdK2dK3;
        
        % Calculate the material contribution to the spatial tangents
        Dtan = zeros(3,3,3,3);
        for mm = 1:3
            for nn = 1:3
                for pp = 1:3
                    for qq = 1:3
                        Dtan(mm,nn,pp,qq) = Dtan(mm,nn,pp,qq) + d2psidK12*delta(mm,nn)*delta(pp,qq) + (d2psidK22 - dpsidK2*(1/K2(ii,jj)) - dpsidK3*(3*K3/(K2(ii,jj)^2)))*N(mm,nn)*N(pp,qq) + (d2psidK32/(K2(ii,jj)^2))*Y(mm,nn)*Y(pp,qq) + (d2psidK1dK2 - 2*sqrt(6)*(dpsidK3/(K2(ii,jj)^2)))*(delta(mm,nn)*N(pp,qq) + N(mm,nn)*delta(pp,qq)) + d2psidK1dK3*(1/K2(ii,jj))*(delta(mm,nn)*Y(pp,qq) + Y(mm,nn)*delta(pp,qq)) + (d2psidK2dK3*(1/K2(ii,jj)) - 3*dpsidK3*(1/(K2(ii,jj)^2)))*(N(mm,nn)*Y(pp,qq) + Y(mm,nn)*N(pp,qq)) + (dpsidK2*(1/K2(ii,jj)) - dpsidK3*(3*K3/(K2(ii,jj)^2)))*((1/2)*delta(mm,pp)*delta(nn,qq) + (1/2)*delta(mm,qq)*delta(nn,pp) - (1/3)*delta(mm,nn)*delta(pp,qq)) + (3*sqrt(6)/2)*dpsidK3*(1/(K2(ii,jj)^2))*(delta(mm,pp)*N(nn,qq) + delta(mm,qq)*N(nn,pp) + N(mm,pp)*delta(nn,qq) + N(mm,qq)*delta(nn,pp));
                    end
                end
            end
        end
        
        % Calculate the derivative of the log(B) with respect to B
        Ltan = dlnxdx(Bstretch);
        
        % Calculate the Btan
        Btan = zeros(3,3,3,3);
        for mm = 1:3
            for nn = 1:3
                for pp = 1:3
                    for qq = 1:3
                        Btan(mm,nn,pp,qq) = Btan(mm,nn,pp,qq) + delta(mm,pp)*Bstretch(nn,qq) + delta(nn,pp)*Bstretch(mm,qq);
                    end
                end
            end
        end
        
        % Calculate the spatial tangents
        Ctan = mprod4(mprod4(Dtan,Ltan),Btan)/(2*J);
        for mm = 1:3
            for nn = 1:3
                for pp = 1:3
                    for qq = 1:3
                        Ctan(mm,nn,pp,qq) = Ctan(mm,nn,pp,qq) - TK(mm,qq)*delta(nn,pp)/J;
                    end
                end
            end
        end
        
        % Calculate the material tangents
        Cmat = zeros(3,3,3,3);
        for mm = 1:3
            for nn = 1:3
                for pp = 1:3
                    for qq = 1:3
                        for rr = 1:3
                            for ss = 1:3
                                Cmat(mm,nn,pp,qq) = Cmat(mm,nn,pp,qq) + J*Finv(nn,rr)*Finv(qq,ss)*Ctan(mm,rr,pp,ss);
                            end
                        end
                    end
                end
            end
        end
        
        if (Cmat(1,1,1,1)<=0)
            %display('1')
            stability(ii,jj) = 1;
        end
        
        if (Cmat(2,2,2,2)<=0)
            %display('2')
            stability(ii,jj) = 1;
        end
        
        if (Cmat(3,3,3,3)<=0)
            %display('3')
            stability(ii,jj) = 1;
        end
        
        if (Cmat(1,2,1,2)<=0)
            %display('4')
            stability(ii,jj) = 1;
        end
        
        if (Cmat(1,3,1,3)<=0)
            %display('5')
            stability(ii,jj) = 1;
        end
        
        if (Cmat(2,3,2,3)<=0)
            %display('6')
            stability(ii,jj) = 1;
        end
        
        temp = Cmat(1,1,1,1)*Cmat(2,2,2,2) + Cmat(1,2,1,2)*Cmat(1,2,1,2) - (Cmat(1,1,2,2) + Cmat(1,2,2,1))^2 + 2*Cmat(1,2,1,2)*sqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2));
        if (temp<=0)
            %display('7')
            stability(ii,jj) = 1;
        end
        
        temp = Cmat(1,1,1,1)*Cmat(3,3,3,3) + Cmat(1,3,1,3)*Cmat(1,3,1,3) - (Cmat(1,1,3,3) + Cmat(1,3,3,1))^2 + 2*Cmat(1,3,1,3)*sqrt(Cmat(1,1,1,1)*Cmat(3,3,3,3));
        if (temp<=0)
            %display('8')
            stability(ii,jj) = 1;
        end
        
        temp = Cmat(2,2,2,2)*Cmat(3,3,3,3) + Cmat(2,3,2,3)*Cmat(2,3,2,3) - (Cmat(2,2,3,3) + Cmat(2,3,3,2))^2 + 2*Cmat(2,3,2,3)*sqrt(Cmat(2,2,2,2)*Cmat(3,3,3,3));
        if (temp<=0)
            %display('9')
            stability(ii,jj) = 1;
        end
        
        deltalist = [1 1 1];
        %
        tempM = zeros(3,3);
        tempM(1,1) = Cmat(1,1,1,1);
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1));
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1));
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2));
        tempM(2,2) = Cmat(2,2,2,2);
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2));
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3));
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3));
        tempM(3,3) = Cmat(3,3,3,3);
        %
        temp = tempM(1,2)*sqrt(Cmat(3,3,3,3)) + tempM(1,3)*sqrt(Cmat(2,2,2,2)) + tempM(2,3)*sqrt(Cmat(1,1,1,1)) + sqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3));
        %
        if(temp<0)
            stability(ii,jj) = 1;
        end
        
        deltalist = [1 1 -1];
        %
        tempM = zeros(3,3);
        tempM(1,1) = Cmat(1,1,1,1);
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1));
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1));
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2));
        tempM(2,2) = Cmat(2,2,2,2);
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2));
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3));
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3));
        tempM(3,3) = Cmat(3,3,3,3);
        %
        temp = tempM(1,2)*sqrt(Cmat(3,3,3,3)) + tempM(1,3)*sqrt(Cmat(2,2,2,2)) + tempM(2,3)*sqrt(Cmat(1,1,1,1)) + sqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3));
        %
        if(temp<0)
            stability(ii,jj) = 1;
        end
        
        deltalist = [1 -1 1];
        %
        tempM = zeros(3,3);
        tempM(1,1) = Cmat(1,1,1,1);
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1));
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1));
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2));
        tempM(2,2) = Cmat(2,2,2,2);
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2));
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3));
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3));
        tempM(3,3) = Cmat(3,3,3,3);
        %
        temp = tempM(1,2)*sqrt(Cmat(3,3,3,3)) + tempM(1,3)*sqrt(Cmat(2,2,2,2)) + tempM(2,3)*sqrt(Cmat(1,1,1,1)) + sqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3));
        %
        if(temp<0)
            stability(ii,jj) = 1;
        end
        
        deltalist = [-1 1 1];
        %
        tempM = zeros(3,3);
        tempM(1,1) = Cmat(1,1,1,1);
        tempM(1,2) = Cmat(1,2,1,2) + deltalist(1)*deltalist(2)*(Cmat(1,1,2,2) + Cmat(1,2,2,1));
        tempM(1,3) = Cmat(1,3,1,3) + deltalist(1)*deltalist(3)*(Cmat(1,1,3,3) + Cmat(1,3,3,1));
        tempM(2,1) = Cmat(2,1,2,1) + deltalist(2)*deltalist(1)*(Cmat(2,2,1,1) + Cmat(2,1,1,2));
        tempM(2,2) = Cmat(2,2,2,2);
        tempM(2,3) = Cmat(2,3,2,3) + deltalist(2)*deltalist(3)*(Cmat(2,2,3,3) + Cmat(2,3,3,2));
        tempM(3,1) = Cmat(3,1,3,1) + deltalist(3)*deltalist(1)*(Cmat(3,3,1,1) + Cmat(3,1,1,3));
        tempM(3,2) = Cmat(3,2,3,2) + deltalist(3)*deltalist(2)*(Cmat(3,3,2,2) + Cmat(3,2,2,3));
        tempM(3,3) = Cmat(3,3,3,3);
        %
        temp = tempM(1,2)*sqrt(Cmat(3,3,3,3)) + tempM(1,3)*sqrt(Cmat(2,2,2,2)) + tempM(2,3)*sqrt(Cmat(1,1,1,1)) + sqrt(Cmat(1,1,1,1)*Cmat(2,2,2,2)*Cmat(3,3,3,3));
        %
        if(temp<0)
            stability(ii,jj) = 1;
        end
        
    end
end

figure(3);
contourf(K1,K2,stability);

%plot(-0.21,0.3094,'*')
