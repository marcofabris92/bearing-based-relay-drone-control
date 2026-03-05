function dotpr = collav(n,epsl,vM,gbi,x,r,gr,ur)

% s = @(angle) sin(angle);
% c = @(angle) cos(angle);
% Rz = @(angle) [c(angle) -s(angle); s(angle) c(angle)];


dr = +Inf;
Vdrstar = [];
for i = 1:n
    dri = norm(x(:,i)-r);
    if dri <= 2*epsl
        if dri < dr
            dr = dri;
            Vdrstar = i;
        elseif dri == dr
            Vdrstar = [Vdrstar i];
        end
    end
end

ar = 0;
eta_eps = 0;
if ~isempty(Vdrstar)
    mingrij = +Inf;
    Vr = [0 0];
    for i = 1:length(Vdrstar)
        for j = 1:length(Vdrstar)
            gri = gr(:,i);
            grj = gr(:,j);
            if gri'*grj < mingrij
                mingrij = gri'*grj;
                Vr = [i j];
            end
        end
    end
    nrbar = gr(:,Vr(1))+gr(:,Vr(2));
    nr = nrbar/norm(nrbar);


    %eta_eps = 1+tanh(1*(epsl^2-dr^2));
    delta = 1/100;
    eta_eps = -dr/(delta*epsl)+(1+delta)/delta;
    if eta_eps > 1
        eta_eps = 1;
    end
    if eta_eps < 0
        eta_eps = 0;
    end

    vMbar = vM/(nr'*gr(:,Vr(1)));
    
    wr = -ur'*nr;
    if wr < vMbar
        ar = (vMbar-wr)/(nr'*gbi);
    end
end

upsilonr = -eta_eps*ar*gbi;
dotpr = ur+upsilonr;


end

