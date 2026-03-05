close all
clc
clear all

global chiN_
chiN_ = NaN; % chiN_ keeps track of control switches

sim01 = 0; % 0 = pseudo-random trajectories, 1 = patrolling example

s = @(angle) sin(angle);
c = @(angle) cos(angle);
Rz = @(angle) [c(angle) -s(angle); s(angle) c(angle)];

gamma = pi/4; % half of the FoV angle
vM = 5; % maximum speed
n = 5; % number of agents
r0 = [0 0]'; % relay initial condition
gbi = [0 -1]'-r0; % bisector
gbi = gbi/norm(gbi);
gfov1 = Rz(-gamma)*gbi; % first FoV vector
gfov2 = Rz(gamma)*gbi; % second FoV vector

x10 = r0+20*gfov1; % agent 1 is initially located on first border
%x20 = r0+30*Rz(-pi/8)*gfov2; % worst case scenario
x20 = r0+27*Rz(pi/100)*gfov1;
x30 = r0+20*Rz(pi/8)*gfov1; % chattering dance
x40 = r0+6*Rz(-pi/6)*gbi; % circle
x50 = [r0+30*Rz(pi/8)*gfov1; -10]; % Lorenz strange attractor

v1 = 0.02*vM*Rz(-pi/2)*gfov1; % velocity v1
v2 = 0.03*vM*Rz(pi/3)*gfov2; % velocity v2
v2(2) = 0;
v3par = [0.9 0.02 1]'; % params for velocity v3
v4par = [0.4 0.4 0.5]; % params for velocity v4
beta = 1+rand;
sigma = beta+rand;
rho = 1+rand;
v5par = [sigma rho beta]'; % params for velocity 5 (sigma,rho,beta)

if sim01 
    vM = 5;
    e1 = [1 0]';
    p.vi1 = vM*e1;
    p.vi2 = vM*Rz(2/3*pi)*e1;
    p.vi3 = vM*Rz(4/3*pi)*e1;
    p.omega = 0.2;
    l = 10*vM;
    rr = (2*sqrt(3)/3)*l;
    p.rr = rr;
    m1 = 20;
    x10 = r0+m1*Rz(-pi/6)*gbi;
    x20 = r0+(m1+l)*Rz(-pi/6)*gbi;
    x30 = r0+(m1+2*l)*Rz(-pi/6)*gbi;
    x40 = x30+l*e1;
    x50 = x10-rr*[0 1]';
    x50 = [x50; 0];
end

Kr_crit = 4*vM/(3*sqrt(2)-2); % critical control gain, valid only for gamma = pi/2
if gamma ~= pi/4
    Kr_crit = vM/(s(gamma)^3); % this is a lower bound for critical gain that works in any case
end
Kr = 1*Kr_crit; % control gain

p.v1 = v1;
p.v2 = v2;
p.v3par = v3par;
p.v4par = v4par;
p.v5par = v5par;
p.vM = vM;
p.gbi = gbi;
p.gfov1 = gfov1;
p.gfov2 = gfov2;
p.Kr = Kr;
p.n = n;
p.epsl = 5;
p.vM = vM;

dt = 0.001;
tspan = (0:dt:30);
x0 = [x10; x20; x30; x40; x50; r0];

if sim01
    [t,y] = ode45(@(t,y) dynNagentsEx(t,y,p),tspan,x0);
else
    [t,y] = ode45(@(t,y) dynNagentsRand(t,y,p),tspan,x0);
end

tt = t;
ttL = length(tt);
x1 = y(1:ttL,1:2)';
x2 = y(1:ttL,3:4)';
x3 = y(1:ttL,5:6)';
x4 = y(1:ttL,7:8)';
x5 = y(1:ttL,9:10)';
r = y(1:ttL,12:13)';

% Fov constraint check + dancing around the bisector
Rz90 = [0 -1; 1 0]';
Pgbi = eye(2)-gbi*gbi';
for kt = 1:length(tt)
    gr1 = x1(:,kt)-r(:,kt);
    gr1 = gr1/norm(gr1);
    gr2 = x2(:,kt)-r(:,kt);
    gr2 = gr2/norm(gr2);
    gr3 = x3(:,kt)-r(:,kt);
    gr3 = gr3/norm(gr3);
    gr4 = x4(:,kt)-r(:,kt);
    gr4 = gr4/norm(gr4);
    gr5 = x5(:,kt)-r(:,kt);
    gr5 = gr5/norm(gr5);
    check1 = gfov1'*Rz90*gr1 >= 0 && gfov2'*Rz90*gr1 <= 0;
    check2 = gfov1'*Rz90*gr2 >= 0 && gfov2'*Rz90*gr2 <= 0;
    check3 = gfov1'*Rz90*gr3 >= 0 && gfov2'*Rz90*gr3 <= 0;
    check4 = gfov1'*Rz90*gr4 >= 0 && gfov2'*Rz90*gr4 <= 0;
    check5 = gfov1'*Rz90*gr5 >= 0 && gfov2'*Rz90*gr5 <= 0;
    if check1*check2*check3*check4*check5 == -1
        fprintf('FoV constraints have been violated at')
        time = t(kt)
    end
    chiN = side([gr1 gr2 gr3 gr4 gr5],gfov1,gfov2,n);
    if chiN ~= chiN_
        if kt > 1
            fprintf(strcat('old_chiN = ',num2str(chiN_),'\n'))
            fprintf(strcat('new_chiN = ',num2str(chiN),'\n\n'))
            time = t(kt)
        end
        chiN_ = chiN;
    end
end


figure
yellow = [255,215,10]/255;
cyan = [0,191,255]/255;
green = [50,205,50]/255;
grid on
hold on
lw = 2.5;
ftsz = 40;
plot(x1(1,:),x1(2,:),'b','linewidth',lw)
plot(x2(1,:),x2(2,:),'m','linewidth',lw)
plot(x3(1,:),x3(2,:),'color',yellow,'linewidth',lw)
plot(x4(1,:),x4(2,:),'color',cyan,'linewidth',lw)
plot(x5(1,:),x5(2,:),'color',green,'linewidth',lw)
plot(r(1,:),r(2,:),'r','linewidth',lw)
set(gca,'FontSize',ftsz);



L = 5;
magquiv = 3; 
if sim01
    magquiv = 10;
end
DeltaT = (tt(end)-tt(1))/L;
hc = 10^3;
ls = linspace(0,1);
seg1 = zeros(2,length(ls));
seg2 = zeros(2,length(ls));
for kl = 0:L
    kt = round(kl*DeltaT/dt);
    if kt < 1
        kt = 1;
    end
    x1_kt = x1(:,kt);
    x2_kt = x2(:,kt);
    x3_kt = x3(:,kt);
    x4_kt = x4(:,kt);
    x5_kt = x5(:,kt);
    r_kt = r(:,kt);
    fov1_kt = r_kt+gfov1*hc;
    fov2_kt = r_kt+gfov2*hc;
    for kk = 1:length(ls)
        seg1(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov1_kt;
        seg2(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov2_kt;
    end

    hfov1 = plot(seg1(1,:),seg1(2,:),':k','linewidth',lw/2);
    hfov2 = plot(seg2(1,:),seg2(2,:),':k','linewidth',lw/2);
    dr1_kt = norm(x1_kt-r_kt);
    hgr1 = quiver(r_kt(1),r_kt(2),x1_kt(1)-r_kt(1),x1_kt(2)-r_kt(2),...
        'color','blue','linewidth',lw,'AutoScaleFactor',magquiv/dr1_kt,...
        'MaxHeadSize',0.5);
    dr2_kt = norm(x2_kt-r_kt);
    hgr2 = quiver(r_kt(1),r_kt(2),x2_kt(1)-r_kt(1),x2_kt(2)-r_kt(2),...
        'color','magenta','linewidth',lw,'AutoScaleFactor',magquiv/dr2_kt,...
        'MaxHeadSize',0.5);
    dr3_kt = norm(x3_kt-r_kt);
    hgr3 = quiver(r_kt(1),r_kt(2),x3_kt(1)-r_kt(1),x3_kt(2)-r_kt(2),...
        'color',yellow,'linewidth',lw,'AutoScaleFactor',magquiv/dr3_kt,...
        'MaxHeadSize',0.5);
    dr4_kt = norm(x4_kt-r_kt);
    hgr4 = quiver(r_kt(1),r_kt(2),x4_kt(1)-r_kt(1),x4_kt(2)-r_kt(2),...
        'color',cyan,'linewidth',lw,'AutoScaleFactor',magquiv/dr4_kt,...
        'MaxHeadSize',0.5);
    dr5_kt = norm(x5_kt-r_kt);
    hgr5 = quiver(r_kt(1),r_kt(2),x5_kt(1)-r_kt(1),x5_kt(2)-r_kt(2),...
        'color',green,'linewidth',lw,'AutoScaleFactor',magquiv/dr5_kt,...
        'MaxHeadSize',0.5);
    
    hx1 = plot(x1(1,kt),x1(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','blue',...
        'MarkerFaceColor','blue');
    hx2 = plot(x2(1,kt),x2(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','magenta',...
        'MarkerFaceColor','magenta');
    hx3 = plot(x3(1,kt),x3(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor',yellow,...
        'MarkerFaceColor',yellow);
    hx4 = plot(x4(1,kt),x4(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor',cyan,...
        'MarkerFaceColor',cyan);
    hx5 = plot(x5(1,kt),x5(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor',green,...
        'MarkerFaceColor',green);
    hr = plot(r(1,kt),r(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor','red');
    text(2+r(1,kt),1+r(2,kt),strcat('$k = ',num2str(kl),'$'),'interpreter','latex','fontsize',ftsz)
end

axis equal
x_1 = [x1(1,:) x2(1,:) x3(1,:) x4(1,:) x5(1,:) r(1,:)]; 
x_2 = [x1(2,:) x2(2,:) x3(2,:) x4(2,:) x5(2,:) r(2,:)];
axx = 40;
if sim01
    axx = 100;
end
xlim([-10+min(x_1) axx+max(x_1)])
ylim([-1+min(x_2) 5+max(x_2)])
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz)
hplots = [hr hx1 hx2 hx3 hx4 hx5 hgr1 hgr2 hgr3 hgr4 hgr5 hfov1];
leg = legend(hplots,{...
    '$\mathbf{p}_{r}(t_{k})$',...
    '$\mathbf{p}_{1}(t_{k})$',...
    '$\mathbf{p}_{2}(t_{k})$',...
    '$\mathbf{p}_{3}(t_{k})$',...
    '$\mathbf{p}_{4}(t_{k})$',...
    '$\mathbf{p}_{5}(t_{k})$',...
    '$\mathbf{g}_{r1}(t_{k})$',...
    '$\mathbf{g}_{r2}(t_{k})$',...
    '$\mathbf{g}_{r3}(t_{k})$',...
    '$\mathbf{g}_{r4}(t_{k})$',...
    '$\mathbf{g}_{r5}(t_{k})$',...
    '$\mathbf{h}_{FoVj}(t_{k}), ~j=1,2$',...
    },'Interpreter','latex','location','northeast');
xlabel('$x$ [m]','interpreter','latex','fontsize',ftsz)
ylabel('$y$ [m]','interpreter','latex','fontsize',ftsz)




