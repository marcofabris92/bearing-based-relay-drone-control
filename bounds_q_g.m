clear all
close all
clc

syms g p

assume(in(g,'real') & g>0 & g<=pi/2)
assumeAlso(in(p,'real') & p>=0 & p<=g)

q_g = sin(g)^3+...                      % >= sin(g)^3
    sin(g)*sin(p-g)^2+...               % >= 0 [p=g]
    cos(g)*sqrt(...                     % >= cos(g)*sqrt()
    sin(g)^2+...                        % >= sin(g)^2
    sin(g-p)^2-...                      % >= 0 [p=g]
    2*cos(2*g-p)*sin(g)*sin(g-p)...     % >= -2*cos(2*g)*sin(g)^2 [p=0]
    -(sin(g)^2+sin(g-p)^2)^2);          % >= -4*sin(g)^4 [p=0]

qq = ((sin(g))^2+(sin(g-p))^2-2*cos(2*g-p)*sin(g)*sin(g-p)-((sin(g))^2+(sin(g-p))^2)^2)-...
    ((sin(p)^2)*(cos(2*g-p)^2));

qq = simplify(qq,'IgnoreAnalyticConstraints',false,'Steps',100);
qq = simplify(qq,'IgnoreAnalyticConstraints',true,'Steps',100);
qq = simplify(qq,'Steps',1000)
simplify(subs(qq,[g p],[pi/3 pi/5]),'Steps',1)

q_g = simplify(q_g,'IgnoreAnalyticConstraints',false);
q_g = simplify(q_g,'Steps',10);

% diff(x*cos(x*y), y, 2)
q_g_prime = simplify(diff(q_g,p),'Steps',10);


u = sin(g)^2+sin(g-p)^2-2*cos(2*g-p)*sin(g)*sin(g-p)...     
    -(sin(g)^2+sin(g-p)^2)^2; 
u = simplify(u,'IgnoreAnalyticConstraints',true);
u = simplify(u,'Steps',10);
simplify(subs(u,p,2*g-pi/2),'Steps',10);

qopt = simplify(subs(q_g,p,3/2*g-pi/4),'Steps',100);
simplify(subs(qopt,g,pi/4),'Steps',100);

qqq = -cos(g)*cos(3*g/2 + pi/4)*abs(cos(g/2 + pi/4)) + (sin(g))^3 - (sin(g))^2/2 + (sin(g))/2;
simplify(qqq,'Steps',1000)


simplify(subs(q_g,p,0),'Steps',1000)










