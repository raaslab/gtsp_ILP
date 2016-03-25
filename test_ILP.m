f = [-3;-2;-1];
intcon = 3;
A = [1,1,1];
b = 7;
Aeq = [4,2,1];
beq = 12;
lb = zeros(3,1);
ub = [Inf;Inf;1]; % enforces x(3) is binary
options = optimoptions('intlinprog','Display','off');


problem = struct('f',f,'intcon',intcon,...
    'Aineq',A,'bineq',b,'Aeq',Aeq,'beq',beq,...
    'lb',lb,'ub',ub,'options',options,...
    'solver','intlinprog');



[x,fval,exitflag,output]  = intlinprog(problem);
