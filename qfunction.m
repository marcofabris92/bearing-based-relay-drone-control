close all
clear all
clc

ftsz = 35;

figure
grid on
hold on

% phi: x must be in [0,gamma]
% gamma: z must be in (0, 90°]
z = [5 25 30 35 40 45 50 55 60 65 70 75 90];
c = @(arg) cos(arg);
s = @(arg) sin(arg);
q = @(xx,zz) c(zz)*s(xx)*abs(c(2*zz-xx))+s(zz)*s(zz-xx)^2+s(zz)^3;
for i = 1:length(z)
    g = z(i);
    x = linspace(0,g,2*10^4);
    y = zeros(1,length(x));
    for j = 1:length(x)
        y(j) = q(x(j)*pi/180,g*pi/180);
    end
    h = plot(x,y,'linewidth',2);
    col = get(h,'Color');
    phistar = 0;
    qstar = 2*s(g*pi/180)^3;
    if g > 30
        phistar = (3*(g*pi/180)/2-pi/4)*180/pi;
        qstar = (3*s(g*pi/180)-1)/2;
    end
    plot([0 phistar],[qstar qstar],'--','color',col)
    plot([phistar phistar],[0 qstar],'--','color',col)
end
for i = 1:length(z)
    g = z(i);
    str = strcat('$\leftarrow \gamma = ',num2str(g),'^{\circ}$');
    text(0,2*s(g*pi/180)^3,str,'fontsize',ftsz,'interpreter','latex')
end

xticks(0:10:90)
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',ftsz);
xlabel('$\phi$ $[^{\circ}]$','fontsize',ftsz,'interpreter','latex')
ylabel('$q_{\gamma}(\phi)$','fontsize',ftsz,'interpreter','latex')
axis([0 90 0 2])




