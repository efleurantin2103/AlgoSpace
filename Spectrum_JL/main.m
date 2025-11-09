clear all;
format long;

eps=0.1;
c=-30;
%lambda = 0.83;
%lambda=0.991052631578947;
lambda = linspace(0.83,1,20);

alpha=0.75;
threshold = eps^(alpha);

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

%% Vectorize the soliton so that it may be made into a function.
%
r = sol.x';
u = sol.y(1,:)';
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

%%
fig = figure('Position', [100, 100, 800, 600]);
hold on

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
[t, W0c] = ode15s('Vop2',[1E7 0], x0, options, flag, c, lambda(h), eps, r, u);

cross_indices = find(diff(W0c(:,1) >= threshold) ~= 0);

Q = cos(W0c(:,2));
eta = sin(W0c(:,2));

re = -Q/sqrt(2) - eta;
im = -eta/sqrt(2) + Q;
top = complex(re,im);%
bottom = complex(re,-im);
det2 = top./bottom;
RE = real(det2);
IM = imag(det2);

for i=1:length(det2)
tht(i) = atan2(IM(i),RE(i));
end
xx=t;
newtht = unwrap(tht);%unwrapping angle values

EndP=abs(newtht(1)-newtht(end))
% newtht(1)
% newtht(end)
jumptau =  round(abs((newtht(1)-newtht(end))/(2*pi)))

%plot(W0c(1:cross_indices,1), newtht(1:cross_indices),'r','LineWidth', 3) %plotting unwrap angle values
plot(W0c(:,1), newtht(1:length(W0c)),'r','LineWidth', 3) %plotting unwrap angle values
xlabel('\tau', 'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'tex')
ylabel('\nu(\tau)', 'FontSize', 24, 'FontWeight', 'bold', 'Interpreter', 'tex')
%xline(threshold,'-');
set(gca, 'FontSize', 20, 'LineWidth', 1.5, 'Box', 'on')
set(gcf, 'color', 'w')
axis([0 1 -inf inf])

% Tight layout
    ax = gca;
    ax.XAxis.TickLength = [0.015 0.015];
    ax.YAxis.TickLength = [0.015 0.015];

% Save as high-resolution PNG
% print(fig, 'psi_det2_highfigure', '-dpng', '-r300')


drawnow
end
