function dydt = dyn1agent(t,y,p)

x1 = y(1:2);
r = y(3:4);

gr1_ = x1-r;
gr1 = gr1_ / norm(gr1_);

Pgr1 = eye(2)-gr1*gr1';

dotpr = -p.Kr*Pgr1*p.gbi;
if p.epsl > 0
    dotpr = collav(p.n,p.epsl,p.vM,p.gbi,x1,r,gr1,dotpr);
end

% if p.epsl > 0 
%     activation = 1+tanh(1*(p.epsl^2-norm(gr1_)^2));
%     if activation > 1
%         activation = 1;
%     end
%     upsr = (p.vM + dotpr'*gr1)*activation;
%     if upsr > 0
%         dotpr = dotpr - upsr*gr1;
%     end
% end

dydt = [p.v1; dotpr];
end

