clear

%Load data
d=load("p_acm.mat");

%STARTING STEADY STATE
d.T=1;
d.tau=0; %tau=0 means steady state

d.Niter=5000; %max iters
d.err_tol=1.0e-8; %tolerance
d.speed=0.4; %update speed

% Make sure: Shapes are correct
d.price=d.price_ss; %initial prices
d.L=ones(d.DJ,d.T)/d.DJ; %guess, total labor=1
d.V=(1/(1-d.bta))*ones(d.DJ,d.T); %guess

disp("Solving initial steady state")
rs=Solve_ACM(d);

%TRANSITION
d.T=30;
d.tau=1; %tau=1 means transition
%Make sure that (1) L(:,1) is the initial labor alloc (2) Shapes are correct
d.price=repmat(d.price_ft,1,d.T); %new prices
d.V=repmat(rs.V,1,d.T); %copy steady state
d.L=repmat(rs.L,1,d.T); %copy steady state

disp("Solving transition")

tic
rt=Solve_ACM(d);
toc

disp("Figures")
Lt=[repmat(rs.L,1,5) rt.L];
waget=[repmat(rs.wage,1,5) rt.wage];

T=d.T;
DJ=d.DJ;

t = (1:(T+5)) - 6;

figure;
plot(t, Lt', 'LineWidth', 1.5); 
grid on;
xlabel('Time');
ylabel('Labor');
title('Labor allocation');

figure;
plot(t, waget', 'LineWidth', 1.5); 
grid on;
xlabel('Time');
ylabel('Wage');
title('Real wages');


%------------------------

function r=Solve_ACM(d)

Niter=d.Niter;
err_tol=d.err_tol;
T=d.T;
tau=d.tau;
speed=d.speed;
price=d.price;
DJ=d.DJ;
nu=d.nu;
C=d.C;
theta=d.theta;
alpha=d.alpha;
psi=d.A;
bta=d.bta;
rho=d.rho;

L=d.L;
V=d.V;
L_star=L;
V_star=V;

err_big=true;
iter=0;
diag_mask_T=repmat(logical(eye(DJ)),1,1,T);
while err_big && iter<Niter

    iter = iter+1;
    
    %SOLUTION ALGORITHM BEGINS 
    %if tau=0 then i_tau=[1 2 3 .. T-1 T]
    %if tau=1 then i_tau=[2 3 .. T-1 T T] thus last period repeats   
    i_tau = [(tau+1):(tau-1)+T T]; %index to use later

    %Step 1: cpi and wages
    cpi = exp(sum(log(price).*theta)); %consumer price index
    wage = (price./cpi).*psi.*alpha.*(L.^(rho-1)).*(alpha.*(L.^rho)+(1-alpha)).^((1./rho)-1); %Eq. 14
    
    %Step 2: moving probabilities
    mij_num = (permute(bta*V(:,i_tau),[3,1,2])-C)./nu; %Eq. 10, two pieces
    mij = exp(mij_num)./sum(exp(mij_num),2); %Eq. 10, two pieces
    
    %Step 3: option values
    omega = -nu*log(reshape(mij(diag_mask_T),DJ,T)); %Eq. 7
    
    %Step 4: implied values
    V_star(:,1:T) = wage + bta.*V(:,i_tau)+omega; %Eq. 3 
    
    %Step 5: implied labor allocations
    L_star(:,tau+1:T) = reshape(sum(permute(L(:,1:T-tau),[1,3,2]).*mij(:,:,1:T-tau),1),DJ,T-tau); %labor flows

    %Step 6: check error and update
    err1 = 100*sum(abs(L(:)-L_star(:)))/T;
    err2 = sum(abs(V(:)-V_star(:)))/T;
    err = [err1 err2];

    L = (1-speed)*L+speed*L_star;
    V = (1-speed)*V+speed*V_star;

    if mod(iter,100)==0
        err_big=(err1 > err_tol || err2 > err_tol );
        disp(iter)
        disp(err)
    end

    %SOLUTION ALGORITHM ENDS

end %while

r.V=V;
r.L=L;
r.wage=wage;
r.mij=mij;
r.err=err;
r.iter=iter;
r.omega=omega;
r.cpi=cpi;

end %function