function [stress1, lambda2, K1, K2] = S_subsub_solve_explicit(lambda1, dt, props, types, branches)

lambda2 = ones(length(lambda1), 1);
stress1 = zeros(length(lambda1), 1);
K2 = zeros(length(lambda1), 1);
K1 = zeros(length(lambda1), 1);



Fet1 = eye(3);
Fvt1 = eye(3);
Fet2 = eye(3);
Fvt2 = eye(3);
Fet3 = eye(3);
Fvt3 = eye(3);
Fet4 = eye(3);
Fvt4 = eye(3);

for i = 1:length(lambda1)
    ite2 = 0;
    if i == 1
        if (lambda1(i) < 1)
            lambda2_up = sqrt(1/lambda1(i))+0.05;
            lambda2_lo =0.85;

        else
            lambda2_lo =sqrt(1/lambda1(i))-0.1;

            lambda2_up = 1.2;

        end
        lambda2_c = 1;
    else
        if (lambda1(i) < 1)    
            lambda2_lo = lambda2_c*0.95;
            lambda2_up = lambda2_c*1.1;

        else
            lambda2_up = lambda2_c*1.05;
            lambda2_lo = lambda2_c*0.9;
        end
    
        ite2= 0;
        max_ite = 5;
        while (ite2<=max_ite)
            lambda2_mi = (lambda2_up + lambda2_lo)/2;
            F_up = [lambda1(i), 0, 0; 0, lambda2_up, 0; 0,0,lambda2_up];
            F_lo = [lambda1(i), 0, 0; 0, lambda2_lo, 0; 0,0,lambda2_lo];
            F_mi = [lambda1(i), 0, 0; 0,lambda2_mi, 0; 0,0, lambda2_mi];
            [~, ~, ~, ~, ~, ~,~,~,stress_mi]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_mi, dt, props, types, branches);
            [~, ~, ~, ~, ~, ~,~,~,stress_lo]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_lo, dt, props, types, branches);
            [~, ~, ~, ~, ~, ~,~,~,stress_up]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_up, dt, props, types, branches);
            if stress_up(2,2)*stress_mi(2,2)<=0
                break;
            elseif stress_lo(2,2)*stress_mi(2,2)<=0
                break;
            else 
                if ite2 == max_ite
                    display('totally out of range');
                    
                end
                if lambda1(i)<1
                    lambda2_lo = lambda2_lo*0.95;
                    lambda2_up = lambda2_up*1.1;
                else
                    lambda2_lo = lambda2_lo*0.9;
                    lambda2_up = lambda2_up*1.05;
                end
            end
            ite2 = ite2+1;
            

        end
    end
    

    
    
    stress2 = 1;
    ite = 0;
    while (abs(stress2) > 1e-5 && ite<2500)
        ite = ite +1;
        lambda2_mi = (lambda2_up + lambda2_lo)/2;
        F_up = [lambda1(i), 0, 0; 0, lambda2_up, 0; 0,0,lambda2_up];
        F_lo = [lambda1(i), 0, 0; 0, lambda2_lo, 0; 0,0,lambda2_lo];
        F_mi = [lambda1(i), 0, 0; 0,lambda2_mi, 0; 0,0, lambda2_mi];
        
        [~, ~, ~, ~, ~, ~,~,~,stress_mi]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_mi, dt, props, types, branches);
        [~, ~, ~, ~, ~, ~,~,~,stress_lo]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_lo, dt, props, types, branches);
        [~, ~, ~, ~, ~, ~,~,~,stress_up]  = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,F_up, dt, props, types, branches);
        
        if (stress_mi(2,2) * stress_up(2,2) <= 0)
            lambda2_lo = lambda2_mi;
        else
            if (stress_mi(2,2) * stress_lo(2,2) <= 0)
                lambda2_up = lambda2_mi;
                
            else
%                 lambda1(i)
%                 ite2
%                 lambda2_lo
%                 lambda2_up
%                 disp('bisection cannot find a solution');
%                 disp(stress_up(2,2))
%                 disp(stress_mi(2,2))
%                 disp(stress_lo(2,2))
%                 disp(stress_mi(3,3))
                break;
            end
            
            stress2 = stress_mi(3,3);
            
        end
    end
    
    lambda2(i) = lambda2_mi;
    lambda2_c  = lambda2_mi;
%    ppp = sprintf('Lamda1 %f, Lamda2 %f, stress %f ,stress22 %f ,stress33 %f ',lambda1(i),lambda2_mi,stress_mi(1,1),stress_mi(2,2),stress_mi(3,3));
%    display(ppp)
    stress1(i) = stress_mi(1,1);
    [Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4, ~] = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4, F_mi, dt, props, types, branches);
    [K1c, K2c] = Ks(F_mi);
%     X1 = props(19);
%     X2 = props(20);
%     K10 = props(17);
%     delta_K = props(18);
%     X_val(i) = X_fun(K1c,X1,X2,K10,delta_K)
    K2(i) = K2c;
    K1(i) = K1c;
    
end

end

function [Fetau1, Fvtau1,Fetau2, Fvtau2,Fetau3, Fvtau3,Fetau4, Fvtau4, stress] = stress_fun(Fet1, Fvt1,Fet2, Fvt2,Fet3, Fvt3,Fet4, Fvt4,  Ftau, dt, props, types, branches)


% initialize to be 0
c = {0, 0, 0,0,0,0,0,0,0,0,0,0};
[Fetau1, Fvtau1,Fetau2, Fvtau2,Fetau3, Fvtau3, Tetau1, Tetau2, Tetau3,Fetau4,Fvtau4,Tetau4] = c{:};

ratio = props(26);
Gneq1 = props(27);
Gneq2 = props(29);
Gneq3 = props(31);
Gneq4 = props(33);
Bneq1 = Gneq1 * ratio;
Bneq2 = Gneq2 * ratio;
Bneq3 = Gneq3 * ratio;
Bneq4 = Gneq4 * ratio;

t1 = props(28);
t2 = props(30);
t3 = props(32) ;
t4 = props(34) ;



if branches(1) == 1
    type = types(1);
    number = 1;
    [Fetau1, Fvtau1,Tetau1] = branch(number, type, Ftau,Fet1, Fvt1, Gneq1, Bneq1, Gneq1*t1, Bneq1*t1, dt, props);
end

if branches(2) == 1
    type = types(2);
    number = 2;
    [Fetau2, Fvtau2,Tetau2] = branch(number,type, Ftau,Fet2, Fvt2, Gneq2, Bneq2, Gneq2*t2, Bneq2*t2, dt, props);
end

if branches(3) == 1
    type = types(3);
    number = 3;
    [Fetau3, Fvtau3,Tetau3] = branch(number,type, Ftau,Fet3, Fvt3, Gneq3, Bneq3, Gneq3*t3, Bneq3*t3, dt, props);
end

if branches(4) == 1
    type= types(4);
    number = 4;
    [Fetau4, Fvtau4,Tetau4] = branch(number,type, Ftau,Fet4, Fvt4, Gneq4, Bneq4, Gneq4*t4, Bneq4*t4, dt, props);
end


Ttau = stress_eq(Ftau, props);

stress = Tetau1 + Tetau2 + Tetau3 + Tetau4 + Ttau ;
stress = inv(Ftau) * stress* det(Ftau);

end

function [Fetau, Fvtau,Tetau] = branch(number, type, Ftau,Fet, Fvt, Gneq, Kneq, eta, kappa,dt, props)
[Me, ~, ~,~ ] = stress_neq(number, Fet, Gneq, Kneq, props);

%shear thining is not used for poron, just give a arbitary value for m, gamma0_dot
m1 = 1;
gamma0_dot_1 = 1;

% shear thining
if type == 2   
    [Dv] = Dv_fun_shear_thinning(Me, m1, gamma0_dot_1);
end
% newtonian 
if type == 1
    [Dv] = Dv_fun_newtonian(Me, eta, kappa);
end



Fvtau = updateFv(Fvt, Dv, dt);
Fetau = Ftau /(Fvtau);
[~, Tetau,~,~] = stress_neq(number,Fetau, Gneq, Kneq, props);

end

function Fvtau = updateFv(Fvt, Dv, dt)
if sum(Dv.^2) == 0
    Fvtau = Fvt;
else
    
    d = Dv*dt;
    expd = eye(3);
    expd(1,1) = exp(d(1,1));
    expd(2,2) = exp(d(2,2));
    expd(3,3) = exp(d(3,3));
    Fvtau = expd * Fvt;
    
end
end

function [Dv] = Dv_fun_newtonian( Me, eta, kappa)

IDEN = eye(3);
Me_tr = (Me(1,1) + Me(2,2) + Me(3,3))/3;
Me_dev = Me - Me_tr*IDEN/3;
Dv = Me_dev/2/eta + Me_tr*IDEN/3/kappa;
end

function  Dv = Dv_fun_shear_thinning(Me, m1, gamma0_dot_1)

IDEN = eye(3);
Me_tr = (Me(1,1) + Me(2,2) + Me(3,3))/3;
Me_dev = Me - Me_tr*IDEN/3;
taue = (sum(sum(Me_dev.^2))/2).^0.5;


Tau0 = 1e5;
sigma0 = 6e5;
eps0_dot_1 = gamma0_dot_1;

gamma_dot_v = gamma0_dot_1*abs(taue/Tau0 )^(1/m1);

eps_dot_v = eps0_dot_1*abs(Me_tr/sigma0)^(1/m1);
%
if (abs(taue) < 1e-8)
    Dv = zeros(3);
else
    
    if (Me_tr > 0.)
        Dv =  Me_dev/taue*0.5*gamma_dot_v + IDEN/3*eps_dot_v;
        
    else
        Dv = Me_dev/taue*0.5*gamma_dot_v - IDEN/3*eps_dot_v;
        
    end
    
end

end


function [Me, Te, K1, K2] = stress_neq(number, F, Gneq, Kneq, props)

Jmin = props(15);
C1 = props(16);

K10 = props(17);
delta_K = props(18);
X1 = props(19);
X2 = props(20);
C0 = props(21);
p = props(22);
q = props(23);
C2 = props(24);
r = props(25);

C3neqs = [props(35),props(36),props(37),props(38)];
C3 = C3neqs(number);
lambda1 = F(1,1);
lambda2 = F(2,2);
lambda3 = F(3,3);
IDEN = eye(3);
E = eye(3);
for i = 1:3
    E(i,i) = log(F(i,i));
end

K1 = log(lambda1*lambda2*lambda3);
dev_E = E - 1/3*K1*IDEN;
K2 = sqrt((log(lambda1)-K1/3)^2+(log(lambda2)-K1/3)^2+(log(lambda3)-K1/3)^2);
if K2 == 0
    Phi = zeros(3);
else
    Phi = dev_E/K2;
end
if (lambda1>=1)
    K3 = 1;
else
    K3 = -1;
end
dfdK1 = dfdK1_fun(K1, Jmin, C2, C3, r);
[X, dXdK1] = X_fun(K1, X1, X2, K10, delta_K);
[dLdK2,~] = dL_fun(K2, K3, C0, C1, p, q);


dWdK1 = dfdK1*Kneq + K2^2*dXdK1*Gneq ;
dWdK2 = (2*K2*X + dLdK2)*Gneq;
%dWdK3 is 0 for simple compression/tension
Me = (dWdK1*IDEN + dWdK2*Phi);
Te = Me/exp(K1);

end




function T = stress_eq(F, props)
G0 = props(1);
B = props(2);
Jmin = props(3);
C1 = props(4);
K10 = props(5);
delta_K = props(6);
X1 = props(7);
X2 = props(8);
C0 = props(9);
p = props(10);
q = props(11) ;
C2 = props(12) ;
C3 = props(13);
r = props(14) ;


lambda1 = F(1,1);
lambda2 = F(2,2);
lambda3 = F(3,3);
IDEN = eye(3);
E = eye(3);
for i = 1:3
    E(i,i) = log(F(i,i));
end

K1 = log(lambda1*lambda2*lambda3);
dev_E = E - 1/3*K1*IDEN;
K2 = sqrt((log(lambda1)-K1/3)^2+(log(lambda2)-K1/3)^2+(log(lambda3)-K1/3)^2);
if K2 == 0
    Phi = zeros(3);
else
    Phi = dev_E/K2;
end
if (lambda1>=1)
    K3 = 1;
else
    K3 = -1;
end

Y = 3*sqrt(6)*Phi*Phi-sqrt(6)*eye(3)-3*K3*Phi;
dWdK3 = C1*(K2^q)*G0/K2;
dfdK1 = dfdK1_fun(K1, Jmin, C2, C3, r);
[X, dXdK1] = X_fun(K1,X1,X2,K10,delta_K);
[dLdK2,~] = dL_fun( K2,K3,C0,C1,p,q);

dWdK1 = dfdK1*B + K2^2*dXdK1*G0 ;
dWdK2 = (2*K2*X + dLdK2)*G0;


%dWdK3 is 0 for simple compression/tension

T = 1/exp(K1)*(dWdK1*IDEN + dWdK2*Phi);
%T = (dWdK1*IDEN + dWdK2*Phi);
end

function [X, dXdK1] = X_fun(K1,X1,X2,K10,delta_K)

dXdK1 = (X1+X2)/2 + (X1-X2)/2*tanh((K1 - K10)/delta_K);

X =  (X1+X2)/2*K1 + (X1-X2)/2*...
    delta_K*log(cosh((K1-K10)/delta_K));

cc =  (X1+X2)/2*0 + (X1-X2)/2*...
    delta_K*log(cosh((0-K10)/delta_K));
X = X - cc+1;

end

function [dLdK2,dLdK3]  = dL_fun(K2,K3,C0,C1,p,q)

dLdK2 = (C1*K2.^(q-1)*q).*(1+K3) + (C0*K2.^(p-1)*p);
dLdK3 = C1*K2.^q;

end

function dfdK1 = dfdK1_fun(K1, Jmin, C2, C3, r)

J = exp(K1);
dfdK1 = (exp(C2*K1)-1)/C2 + C3*(J^(-r)-(1-Jmin)^(r)/(J-Jmin)^r) * exp(K1);

end

function [K1, K2, Phi, K3] = Ks(F)

lambda1 = F(1,1);
lambda2 = F(2,2);
lambda3 = F(3,3);

K1 = log(det(F));
K2 = sqrt((log(lambda1)-K1/3)^2+(log(lambda2)-K1/3)^2+(log(lambda3)-K1/3)^2);

E = eye(3);
for i = 1:3
    E(i,i) = log(F(i,i));
end

dev_E = E - 1/3*K1*eye(3);
if K2 == 0
    Phi = zeros(3);
else
    Phi = dev_E/K2;
end

K3 = 3*6^0.5*det(Phi);

end



