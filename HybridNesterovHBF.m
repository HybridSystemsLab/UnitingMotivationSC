%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: HeavyBall.m
%--------------------------------------------------------------------------
% Project: Testing out parameters lambda and gamma for fast, oscillatory
% convergence globally and slow, smooth convergence locally.
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00

clear all

set(0,'defaultTextInterpreter','tex');

% global variables
global delta M gamma lambda c_0 c_10 z1Star zeta cTilde_0 cTilde_10 d_0 d_10 alpha d beta a 

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting the parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
setMinima();

% Nesterov constants
kappa = 1; % kappa >= 1, kappa = M/mu = 2/2
M = 2;
d = 1/(sqrt(kappa) + 1);
beta = (sqrt(kappa) - 1)/(sqrt(kappa) + 1);
a = d + beta/(2*kappa);

% Heavy Ball constants
gamma = 2/3; 
lambda = 40; 

% Uniting parameters for \mathcal{U}_0 and \mathcal{T}_{1,0}:
c_0 = 3000;  
c_10 = 402.83351565;  % 400 

alpha = 1;

% eps_0 has to be bigger than eps_10
eps_0 = 20;
eps_10 = 15;

cTilde_0 = eps_0*alpha
cTilde_10 = eps_10*alpha
d_0 = c_0 - gamma*((cTilde_0^2)/alpha)
d_10 = c_10 - (a^2)*(cTilde_10/alpha)^2 - (1/M)*((cTilde_10^2)/alpha)

delta = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setting intial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaVecNest = [0,0,0];
deltaVecHBF = [0,0,0];
deltaVecUniting = [0,0,0];
deltaVecUniting2 = [0,0,0];

lDeltaNest = 0;
lDeltaHBF = 0;
lDeltaUniting = 0;
lDeltaUniting2 = 0;

% initial conditions
z1_0 = 50;
z2_0 = 0;
q_0 = 1;

% Assign initial conditions to vector
x0 = [z1_0;z2_0;q_0];
x00 = [z1_0;z2_0;q_0];

% simulation horizon
TSPAN=[0 10]; % 800
TSPAN_HBF =[0 400];
JSPAN = [0 5000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

%% Simulate
[tNest,jNest,xNest] = HyEQsolver(@fNesterov,@gNesterov,@CNesterov,@DNesterov,...
    x0,TSPAN,JSPAN,rule,options);

% Find the L values for Nesterov:
lNest = zeros(1,length(tNest));
[deltaVecNest lNest lDeltaNest] = timeToConv(xNest,tNest);

lNestAvg = lNest.';
% This is the dotted line indicating the average value:
PO = polyfit(tNest,log10(lNestAvg(:,1)),1);
yO = polyval(PO,tNest);

[tHBF,jHBF,xHBF] = HyEQsolver(@fHBF,@gHBF,@CHBF,@DHBF,...
    x00,TSPAN_HBF,JSPAN,rule,options);

% Find the L values for HBF:
lHBF = zeros(1,length(tHBF));
[deltaVecHBF lHBF lDeltaHBF] = timeToConv(xHBF,tHBF);

% Finally, simulate the hybrid closed-loop heavy ball system
[tUniting,jUniting,xUniting] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

% Find the L values for Uniting:
lUniting = zeros(1,length(tUniting));
[deltaVecUniting lUniting lDeltaUniting] = timeToConv(xUniting,tUniting);

% Modifying c_10 for profile comparison:
c_10 = 402.834; % 
d_10 = c_10 - (a^2)*(cTilde_10/alpha)^2 - (1/M)*((cTilde_10^2)/alpha)
[tUniting2,jUniting2,xUniting2] = HyEQsolver(@fU,@gU,@CU,@DU,...
    x0,TSPAN,JSPAN,rule,options);

lUniting2 = zeros(1,length(tUniting2));
[deltaVecUniting2 lUniting2 lDeltaUniting2] = timeToConv(xUniting2,tUniting2);

%% Plots
figure(1) % position
clf
modificatorF{1} = '';
modificatorF{2} = 'LineWidth';
modificatorF{3} = 1.5;
modificatorJ{1} = '*--';
modificatorJ{2} = 'LineWidth';
modificatorJ{3} = 1.5;
subplot(3,1,1), plotHarc(tHBF,jHBF,xHBF(:,1),[],modificatorF,modificatorJ);
axis([0 10 -20 60])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
axes('Position',[0.15 0.8 0.15 0.08])
box on
plot(tHBF,xHBF(:,1),'LineWidth',1.5)
set(gca,'xtick',[0 200 400])
set(gca,'ytick',[-20 20 60])
grid on
axis([0 400 -20 60])
subplot(3,1,2), plotHarc(tNest,jNest,xNest(:,1),[],modificatorF,modificatorJ);
axis([0 10 -20 60])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
subplot(3,1,3), plotHarc(tUniting,jUniting,xUniting(:,1),[],modificatorF,modificatorJ);
axis([0 10 -20 60])
grid on
ylabel('z_1','Fontsize',16)
xlabel('t','Fontsize',16)
saveas(gcf,'Plots\UnitingMotivation','png')
saveas(gcf,'Plots\UnitingMotivation','epsc')

figure(2)
clf
semilogy(tUniting,lUniting,'LineWidth',1.5);
hold on
semilogy(tUniting2,lUniting2,'Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
semilogy(tHBF,lHBF,'Color',[0.4660 0.6740 0.1880],'LineWidth',1.5);
semilogy(tNest,lNest,'Color',[0.4940 0.1840 0.5560],'LineWidth',1.5,'LineStyle','--');
semilogy(tNest,10.^(yO(:,1)),'k--','LineWidth',1.5);
hold off
axis([0 10 10^(-20) 10^(6)]);
ylabel('L(z_1)-L^*','FontSize',20)
xlabel('t','FontSize',20)
legend({'Hybrid, tuning 1','Hybrid, tuning 2','Heavy ball','Nesterov','Nesterov, average'},'Location','southeast')
saveas(gcf,'Plots\Semilog','epsc')