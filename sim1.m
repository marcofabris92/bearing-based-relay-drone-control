close all
clc
clear all

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
% x10 = r0+20*gbi; % agent 1 is initially located on first border
x10 = r0+20*gfov1; % other simulations on the gfov1 border

v1 = vM*Rz(-pi/2)*gfov1; % velocity v1

Kr_crit = vM/s(gamma); % critical control gain
Kr = 1.5*Kr_crit; % control gain

p.n = 1;
p.v1 = v1;
p.gbi = gbi;
p.Kr = Kr;
p.epsl = 5;
p.vM = vM;

dt = 0.001;
tspan = (0:dt:30);
[t,y] = ode15s(@(t,y) dyn1agent(t,y,p),tspan,[x10; r0]);
if Kr == 1.5*Kr_crit
    p.epsl = 0;
    [tw,yw] = ode45(@(t,y) dyn1agent(t,y,p),tspan,[x10; r0]);
end

x1 = y(:,1:2)';
r = y(:,3:4)';
rw = yw(:,3:4)';

tt = t;

% Fov constraint check + dancing around the bisector
Rz90 = [0 -1; 1 0]';
Pgbi = eye(2)-gbi*gbi';
violation = 0;
strong_violation = 0;
for kt = 1:length(tt)
    gr1 = x1(:,kt)-r(:,kt);
    gr1 = gr1/norm(gr1);
    check1 = gfov1'*Rz90*gr1 >=0 && gfov2'*Rz90*gr1 <=0;
    if check1 == 0 
        if gfov1'*Rz90*gr1 < -10*eps
            gfov1'*Rz90*gr1
            fprintf('FoV constraints have been violated at')
            time = t(kt)
            strong_violation = 1;
        end
        %violation = 1;
    end
end
%violation
strong_violation

mindr1 = +Inf;
mintk = -1;
for tk = 1:length(r(1,:))
    dr1 = norm(x1(:,tk)-r(:,tk));
    if dr1 < mindr1
        mindr1 = dr1;
        mintk = t(tk);
    end
end
mindr1
mintk

figure
grid on
hold on
lw = 2.5;
if Kr == 1.5*Kr_crit
    plot(rw(1,:),rw(2,:),'g','linewidth',2*lw)
end
plot(x1(1,:),x1(2,:),'b','linewidth',lw)
plot(r(1,:),r(2,:),'r','linewidth',lw)
set(gca,'FontSize',ftsz);



L = 5;
DeltaT = (tspan(end)-tspan(1))/L;
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
        'color','blue','linewidth',lw,'AutoScaleFactor',5/dr1_kt,...
        'MaxHeadSize',0.5);
    
    if Kr == 1.5*Kr_crit
        hrw = plot(rw(1,kt),rw(2,kt),'-s','MarkerSize',15,...
            'MarkerEdgeColor','green',...
            'MarkerFaceColor','green');
    end
    hx1 = plot(x1(1,kt),x1(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','blue',...
        'MarkerFaceColor','blue'); 
    hr = plot(r(1,kt),r(2,kt),'-s','MarkerSize',10,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor','red');
    text(r(1,kt),5+r(2,kt),strcat('$k = ',num2str(kl),'$'),'interpreter','latex','fontsize',ftsz)
end
axis equal
xlim([-10+min(x1(1,:)) 50])
ylim([-30 10+max(r(2,:))])
set(gca,'TickLabelInterpreter','latex', 'FontSize',ftsz)
hplots = [hr hx1 hgr1 hfov1];
leg = legend(hplots,{...
    '$\mathbf{p}_{r}(t_{k})$',...
    '$\mathbf{p}_{1}(t_{k})$',...
    '$\mathbf{g}_{r1}(t_{k})$',...
    '$\mathbf{h}_{FoVj}(t_{k}), ~ j=1,2$',...
    },'Interpreter','latex','location','northeast');
if Kr == 1.5*Kr_crit
    hplots = [hrw hr hx1 hgr1 hfov1];
    leg = legend(hplots,{...
        '$\mathbf{p}_{r}(t_{k})$ w/o coll. av.',...
        '$\mathbf{p}_{r}(t_{k})$',...
        '$\mathbf{p}_{1}(t_{k})$',...
        '$\mathbf{g}_{r1}(t_{k})$',...
        '$\mathbf{h}_{FoVj}(t_{k}), ~ j=1,2$',...
        },'Interpreter','latex','location','northeast');
end
xlabel('$x$ [m]','interpreter','latex','fontsize',ftsz)
ylabel('$y$ [m]','interpreter','latex','fontsize',ftsz)