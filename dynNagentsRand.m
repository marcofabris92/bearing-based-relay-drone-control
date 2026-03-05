function dydt = dynNagentsRand(t,y,p)

n = p.n;
r = y(2*n+2:2*n+3);
x = zeros(2,n);
gr = zeros(2,n);
for i = 0:n-1
    ii = i+1;
    x(:,ii) = y(2*i+1:2*i+2);
    gr(:,ii) = x(:,ii)-r;
    gr(:,ii) = gr(:,ii)/norm(gr(:,ii));
end

v1 = p.v1;
v2 = p.v2;
v3par = p.v3par;
v4par = p.v4par;
v5par = p.v5par;
vM = p.vM;

Kr = p.Kr;
gbi = p.gbi;
gfov1 = p.gfov1;
gfov2 = p.gfov2;

profit1 = -Inf;
k1 = 0;
profit2 = -Inf;
k2 = 0;
for i = 1:n
    profit1_i = gr(:,i)'*gfov1;
    profit2_i = gr(:,i)'*gfov2;
    if profit1_i > profit1
        profit1 = profit1_i;
        k1 = i;
    end
    if profit2_i > profit2
        profit2 = profit2_i;
        k2 = i;
    end
end

k = k1;
if profit2 > profit1
    k = k2;
end

chiN = side(gr,gfov1,gfov2,n);

Pgrk = eye(2)-gr(:,k)*gr(:,k)';
ur = -Kr*Pgrk*gbi;
k
if chiN < 0
    k1
    k2
    Pgr_k1 = eye(2)-gr(:,k1)*gr(:,k1)';
    Pgr_k2 = eye(2)-gr(:,k2)*gr(:,k2)';
    ur = -Kr*(Pgr_k1+Pgr_k2)*gbi;
end

% dynamics of agent 2: piecewise constant
xdot2 = v2;
if t > 26
    xdot2 = t*v2;
end

% dynamics of agent 3: dancing around the bisector
Pgbi = eye(2)-gbi*gbi';
b3 = v3par(1);
c3 = v3par(2);
a3 = sqrt(1-b3^2-c3^2);
omega3 = v3par(3);
xdot3 = -a3*vM*Pgbi*gr(:,3)-...
    b3*vM*sin(omega3*t)*[1 0]'+...
    c3*vM*[0 1]';

% dynamics of agent 4: circle around the bisector
b4 = v4par(1);
c4 = v4par(2);
a4 = sqrt(1-b4^4-c4^2);
omega4 = v4par(3);
xdot4 = -a4*vM*Pgbi*gr(:,4)-...
    b4*vM*cos(omega4*t)*[1 0]'-...
    c4*vM*sin(omega4*t)*[0 1]';

% dynamics of agent 5: Lorenz strange attractor
xl = y(2*n+1);
yl = x(1,5);
zl = x(2,5);
sigmal = v5par(1);
rhol = v5par(2);
betal = v5par(3);
xdot5 = [xl*(rhol-zl)-yl; xl*yl-betal*zl; sigmal*(yl-xl)];
xdot5 = 0.001*xdot5-[0.5*vM*Pgbi*gr(:,5); 0];

dotpr = ur;
if p.epsl > 0
    dotpr = collav(n,p.epsl,p.vM,gbi,x,r,gr,ur);
end

% if p.epsl > 0 
%     mindist = +Inf;
%     greps = zeros(2,1);
%     for ii = 1:n
%         if norm(x(:,ii)-r) < mindist
%             mindist = norm(x(:,ii)-r);
%             greps = greps + gr(:,ii);
%         end
%     end
%     greps = greps / norm(greps);
%     activation = 1+tanh(1*(p.epsl^2-mindist^2));
%     if activation > 1
%         activation = 1;
%     end
%     upsr = (p.vM + dotpr'*greps)*activation;
%     if upsr > 0
%         dotpr = dotpr - upsr*greps;
%     end
% end


dydt = [v1; xdot2; xdot3; xdot4; xdot5; dotpr];
end

