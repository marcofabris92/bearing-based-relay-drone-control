close all
clc
clear all

sim01 = 0; % 0 = worst case scenario, 1 = chattering dance

global chi2_
chi2_ = NaN; % chi2_ keeps track of control switches

ftsz = 40;

s = @(angle) sin(angle);
c = @(angle) cos(angle);
Rz = @(angle) [c(angle) -s(angle); s(angle) c(angle)];

gamma = pi/4; % half of the FoV angle
vM = 5; % maximum speed

r0 = [0 0]'; % relay initial condition
gbi = [0 -1]'-r0; % bisector
gbi = gbi/norm(gbi);
gfov1 = Rz(-gamma)*gbi; % first FoV vector
gfov2 = Rz(gamma)*gbi; % second FoV vector
x10 = r0+20*gfov1; % agent 1 is initially located on first border
x20 = r0+30*Rz(-pi/8)*gfov2; % worst case scenario
if sim01
    x20 = r0+30*gfov1; % chattering dance
end

kv1 = 1; % = 1 for worst case scenario, = 0.1 chattering dance
if sim01
    kv1 = 0.1;
end
v1 = kv1*vM*Rz(-pi/2)*gfov1; % velocity v1
v2 = vM*Rz(pi/2)*gfov2; % velocity v2

Kr_crit = 4*vM/(3*sqrt(2)-2); % critical control gain, valid only for gamma = pi/2
if gamma ~= pi/4
    Kr_crit = vM/(s(gamma)^3); % this is a lower bound for critical gain that works in any case
end
Kr = 1*Kr_crit; % control gain

p.n = 2;
p.v1 = v1;
p.v2 = v2;
p.gbi = gbi;
p.gfov1 = gfov1;
p.gfov2 = gfov2;
p.Kr = Kr;
p.epsl = 5;
p.vM = vM;

dt = 0.0001;
tspan = (0:dt:30);
if sim01
   [t,y] = ode45(@(t,y) dyn2agentsdance(t,y,p),tspan,[x10; x20; r0]);
else
    [t,y] = ode45(@(t,y) dyn2agents(t,y,p),tspan,[x10; x20; r0]);
end

tt = t;
ttL = length(tt);
x1 = y(1:ttL,1:2)';
x2 = y(1:ttL,3:4)';
r = y(1:ttL,5:6)';
% Fov constraint check + dancing around the bisector
Rz90 = [0 -1; 1 0]';
Pgbi = eye(2)-gbi*gbi';
for kt = 1:length(tt)
    gr1 = x1(:,kt)-r(:,kt);
    gr1 = gr1/norm(gr1);
    gr2 = x2(:,kt)-r(:,kt);
    gr2 = gr2/norm(gr2);
    check1 = gfov1'*Rz90*gr1 >= 0 && gfov2'*Rz90*gr1 <= 0;
    check2 = gfov1'*Rz90*gr2 >= 0 && gfov2'*Rz90*gr2 <= 0;
    if check1*check2 == 0
        fprintf('FoV constraints have been violated at')
        time = t(kt)
    end
    chi2 = sign(gr1'*Pgbi*gr2);
    if chi2 ~= chi2_
        if kt > 1
            fprintf(strcat('old_chi2 = ',num2str(chi2_),'\n'))
            fprintf(strcat('new_chi2 = ',num2str(chi2),'\n\n'))
            time = t(kt)
        end
        chi2_ = chi2;
    end
end


figure
grid on
hold on
lw = 2.5;
tp = 1:60000;
% tp is used instead : in the plot to stop at t1 
plot(x1(1,:),x1(2,:),'b','linewidth',lw)
plot(x2(1,:),x2(2,:),'m','linewidth',lw)
plot(r(1,:),r(2,:),'r','linewidth',lw)
set(gca,'FontSize',ftsz);



L = 5;
magquiv = 20; 
if sim01
    magquiv = 5;
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
    r_kt = r(:,kt);
    fov1_kt = r_kt+gfov1*hc;
    fov2_kt = r_kt+gfov2*hc;
    for kk = 1:length(ls)
        seg1(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov1_kt;
        seg2(:,kk) = ls(kk)*r_kt+(1-ls(kk))*fov2_kt;
    end
    hfov1 = plot(seg1(1,:),seg1(2,:),':k','linewidth',lw/1.5);
    hfov2 = plot(seg2(1,:),seg2(2,:),':k','linewidth',lw/1.5);
    dr1_kt = norm(x1_kt-r_kt);
    hgr1 = quiver(r_kt(1),r_kt(2),x1_kt(1)-r_kt(1),x1_kt(2)-r_kt(2),...
        'color','blue','linewidth',lw,'AutoScaleFactor',magquiv/dr1_kt,...
        'MaxHeadSize',0.5);
    dr2_kt = norm(x2_kt-r_kt);
    hgr2 = quiver(r_kt(1),r_kt(2),x2_kt(1)-r_kt(1),x2_kt(2)-r_kt(2),...
        'color','magenta','linewidth',lw,'AutoScaleFactor',magquiv/dr2_kt,...
        'MaxHeadSize',0.5);
    
    hx1 = plot(x1(1,kt),x1(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','blue',...
        'MarkerFaceColor','blue');
    hx2 = plot(x2(1,kt),x2(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','magenta',...
        'MarkerFaceColor','magenta');
    hr = plot(r(1,kt),r(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor','red');
    text(r(1,kt),5+r(2,kt),strcat('$k = ',num2str(kl),'$'),'interpreter','latex','fontsize',ftsz)
end
axis equal
axx = 80;
if sim01
    axx = 30;
end
%xlim([-10+min(x1(1,:)) axx+max(x2(1,:))])
xlim([-50 60])
ylim([-30 60]) %10+max(r(2,:))])
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz)
hplots = [hr hx1 hx2 hgr1 hgr2 hfov1];
%'$\mathbf{p}_{r}(t_{k})+\lambda\mathbf{g}_{FoV1},~ \lambda \geq 0$',...
%'$\mathbf{p}_{r}(t_{k})+\lambda\mathbf{g}_{FoV2},~ \lambda \geq 0$',...
leg = legend(hplots,{...
    '$\mathbf{p}_{r}(t_{k})$',...
    '$\mathbf{p}_{1}(t_{k})$',...
    '$\mathbf{p}_{2}(t_{k})$',...
    '$\mathbf{g}_{r1}(t_{k})$',...
    '$\mathbf{g}_{r2}(t_{k})$',...
    '$\mathbf{h}_{FoVj}(t_{k}), ~j=1,2$',...
    },'Interpreter','latex','location','northeast');
xlabel('$x$ [m]','interpreter','latex','fontsize',ftsz)
ylabel('$y$ [m]','interpreter','latex','fontsize',ftsz)