clear all;
format long;

eps=0.1;
c=-30;
lambda = linspace(0.83,1,20);

kappa=0.65;%works fine, relationship to alpha
threshold = eps^(kappa)/(1+eps^(kappa));
threshold2 = eps^(kappa-1)/(1+eps^(kappa-1));

%%%%%%%%%%%%%%%%%%%%%%%%
%% BEGIN COMPUTE SOLITON
%

%% Set soliton amplitude guess "amp0", interval [x0,xmax0], and initial # of points in mesh "npts0"
%
amp0 = 4.3373876;

x0 = 0;
xmax0 = 15;
npts0 = 20;

%% Set further parameters needed for extending soliton to [xmax0,xmax1]
%
deltaxmax = 5;
xmax1  = 750;
nmax = floor((xmax1-xmax0)/deltaxmax);
xmax = xmax0;

%% Create initial guess
%
guessinit = @(x)solit_guess(amp0,x);
xmesh = linspace(x0,xmax0,npts0);
sol = bvpinit(xmesh, guessinit);

%% Set singular term "S" and call bvp4c to calculate soliton on [x0,xmax0]
% 
S = [0 0; 0 -2];

solit_sys_handle = @(x,y)solit_sys(c,eps,x,y);
solit_jac_handle = @(x,y)solit_jac(c,eps,x,y);

solit_bc_handle = @(ya,yb)solit_bc(xmax0,ya,yb);
solit_bc_jac_handle = @(ya,yb)solit_bc_jac(xmax0,ya,yb);

options = bvpset('SingularTerm',S,'FJacobian',solit_jac_handle,'BCJacobian',solit_bc_jac_handle);    
sol = bvp4c(solit_sys_handle, solit_bc_handle, sol, options);

%% Iterate & extend to [x0,xmax1]
%
for i=1:nmax
    xmax = xmax + deltaxmax;
    sol = bvpinit(sol,[x0,xmax]);
    
    solit_sys_handle = @(x,y)solit_sys(c,eps,x,y);
    solit_jac_handle = @(x,y)solit_jac(c,eps,x,y);

    solit_bc_handle = @(ya,yb)solit_bc(xmax,ya,yb);
    solit_bc_jac_handle = @(ya,yb)solit_bc_jac(xmax,ya,yb);

    options = bvpset('SingularTerm',S,'FJacobian',solit_jac_handle,'BCJacobian',solit_bc_jac_handle);    
    sol = bvp4c(solit_sys_handle, solit_bc_handle, sol, options);
end

%% Plot soliton
%
% plot(sol.x, sol.y(1,:));
% axis([0 10 0 inf]);

%% Vectorize the soliton so that it may be made into a function.
%
r = sol.x';
u = sol.y(1,:)';
%return
%% Clear all un-necessary values
%
clear amp0 solit_bc_jac_handle deltaxmax guessinit i;
clear solit_jac_handle nmax npts0 options S solit_bc_handle;
clear solit_sys_handle x0 xmax0 xmax1 xmax xmesh sol;

%% END COMPUTE SOLITON
%
%%%%%%%%%%%%%%%%%%%%%%%%

eq2 = [1;0;0;0;0;0];%IC/Fp starting at sigma=0

t0=0;
tmax0 = 360;

S = zeros(6,6);
S(2,2) = -2;
S(3,3) = -2;
S(4,4) = -2;
S(5,5) = -4;
S(6,6) = 1;

[Q, Lambda] = eig(S);

%assigning values
lambda_u = Lambda(1,6);
xi_u = Q(:, 6);

%%

for h=1:length(lambda)
gamma = ((1-lambda(h))/eps^2);
Fpb=atan(-sqrt(gamma));
Fp1=[1;Fpb];

DF=jac2(Fp1,c,lambda(h),eps);
% Compute the eigenvectors and eigenvalues of the Jacobian
[V1, D] = eig(DF);

% Extract the center/stable eigenvector
if D(1,1) < D(2,2)    
    V1 = V1(:,1);
else
    V1 = V1(:,2);
end

scale = 10^(-4);
x0 = Fp1 - scale*V1;

options=odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, W0c] = ode45('Vop2',[0 -1E6], x0, options, flag, c, lambda(h), eps);


x1 = eq2 + scale*xi_u;

[t1, W0u2] = ode45('Vop',[t0,tmax0], x1, options, flag, c, lambda(h), eps, r, u);

cross_indices = find(diff(W0c(:,1) >= threshold) ~= 0);
cross_indices2 = find(diff(W0u2(:,6) >= threshold2) ~= 0);

Q = cos(W0c(:,2));
eta = sin(W0c(:,2));

re = -Q/sqrt(2) - eta;
im = -eta/sqrt(2) - Q;
top = complex(re,-im);%
bottom = complex(re,im);
det2 = top./bottom;
RE = real(det2);
IM = imag(det2);

for i=1:length(det2)
tht(i) = atan2(IM(i),RE(i));
end
xx=t;
newtht = unwrap(tht);%unwrapping angle values
EndP=abs(newtht(1)-newtht(end));

%For unstable
%Starting det^2 computation
re2 = W0u2(:,1)-W0u2(:,5);
im2 = W0u2(:,2)-W0u2(:,3);
top2 = complex(re2,im2);
bottom2 = complex(re2,-im2);
det21 = top2./bottom2;
RE2 = real(det21);
IM2 = imag(det21);

for i=1:length(det21)
tht2(i) = atan2(IM2(i),RE2(i));
end
xx2=t1;
newtht2 = unwrap(tht2);%unwrapping angle values


figure(1)
set(gcf, 'Color', 'w', 'Position', [100, 100, 600, 400])

LL=length(W0u2(:,6));
LC=length(W0c(:,1));

% Left subplot
subplot(1, 2, 1)
hold on
plot(W0u2(:,6), newtht2(1:LL)+ 6*pi,'LineWidth', 1) % plotting unwrap angle values
plot(W0u2(:,6), newtht2(1:LL)+ 4*pi,'LineWidth', 1) % plotting unwrap angle values
plot(W0u2(:,6), newtht2(1:LL)+ 2*pi,'LineWidth', 1) % plotting unwrap angle values
xline(1,'-', 'LineWidth', 1.5);
xlabel('$\sigma$', 'Interpreter', 'latex')
ylabel('$\det^2$', 'Interpreter', 'latex')
%xlim([0 1])
axis([0 1 -8 20])
set(gca, 'FontSize', 20, 'TickLabelInterpreter', 'latex')
hold off

% Right subplot  
subplot(1, 2, 2)
hold on
plot(W0c(:,1), newtht(1:LC)-2*pi,'r','LineWidth', 1)
xline(0,'-','LineWidth',1.5);
xlabel('$\tau$','Interpreter','latex')
ylabel('$\det^2$','Interpreter','latex')
set(gca, 'YAxisLocation', 'right')  % This moves the y-axis to the right
%xlim([0 1])
axis([0 1 -8 20])
set(gca,'FontSize',20,'TickLabelInterpreter','latex')
hold off

end

%print(1, '-djpeg', '-r300', 'JL2.jpg')
