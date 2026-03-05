function dydt = dynNagentsEx(t,y,p)

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

%vM = p.vM;
om = p.omega;
rr = p.rr;
vi1 = p.vi1;
vi2 = p.vi2;
vi3 = p.vi3;

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

% dyn agent 1
xdot1 = [0 0]';
xdot1(1) = -om*rr*cos(om*t);
xdot1(2) = -om*rr*sin(om*t);

% dyn agent 2-4
if t>=0 && t<10
    xdot2 = vi1;
    xdot3 = vi1;
    xdot4 = vi1;
end
if t>=10 && t<20
    xdot2 = vi2;
    xdot3 = vi2;
    xdot4 = vi2;
end
if t>=20 && t<=30
    xdot2 = vi3;
    xdot3 = vi3;
    xdot4 = vi3;
end

% dyn agent 5
xdot5 = [0 0 0]';

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

dydt = [xdot1; xdot2; xdot3; xdot4; xdot5; dotpr];
end

