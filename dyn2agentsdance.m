function dydt = dyn2agentsdance(t,y,p)

x1 = y(1:2);
x2 = y(3:4);
x = [x1 x2];
r = y(5:6);
gr = [(x1-r)/norm(x1-r) (x2-r)/norm(x2-r)];

v1 = p.v1;
v2 = p.v2;
Kr = p.Kr;
gbi = p.gbi;
gfov1 = p.gfov1;
gfov2 = p.gfov2;

gr1 = gr(:,1);
gr2 = gr(:,2);

Pgbi = eye(2)-gbi*gbi';
chi2 = sign(gr1'*Pgbi*gr2);
Pgr1 = eye(2)-gr1*gr1';
Pgr2 = eye(2)-gr2*gr2';

ur = -Kr*(Pgr1+Pgr2)*gbi;
if chi2 >= 0
    gr1gfovj = max(gr1'*gfov1,gr1'*gfov2);
    gr2gfovj = max(gr2'*gfov1,gr2'*gfov2);
    grbar = gr1;
    if gr2gfovj > gr1gfovj
        grbar = gr2;
    end
    Pgrbar = eye(2)-grbar*grbar';
    ur = -Kr*Pgrbar*gbi;
end

omega = 0.5;
b = 0.9;
c = 0.03;
a = sqrt(1-b^2-c^2);
x2dot = -a*norm(v2)*Pgbi*gr2+...
    b*norm(v2)*sin(omega*t)*[1 0]'+...
    c*norm(v2)*[0 1]';

dotpr = ur;
if p.epsl > 0
    dotpr = collav(p.n,p.epsl,p.vM,gbi,x,r,gr,ur);
end

% if p.epsl > 0 
%     mindist = norm(gr1_);
%     greps = gr1;
%     if norm(gr2_) < mindist
%         mindist = norm(gr2_);
%         greps = gr2;
%     elseif norm(gr2_) == mindist
%         greps = greps + gr2;
%         greps = greps / norm(greps);
%     end
%     activation = 1+tanh(1*(p.epsl^2-mindist^2));
%     if activation > 1
%         activation = 1;
%     end
%     upsr = (p.vM + dotpr'*greps)*activation;
%     if upsr > 0
%         dotpr = dotpr - upsr*greps;
%     end
% end

dydt = [v1; x2dot; dotpr];
end



