s = [1 1 1 2 2 3 3 4 4 5 5 6 7];
t = [2 4 7 3 4 5 6 5 7 8 6 8 8];
weights = [2 1 1 1 1 1 2 1 1 1 1 1 2];
names = {'1' '2' '3' '4' '5' '6' '7' '8'};
G = graph(s,t,weights);

plot(G, 'EdgeLabel', weights);

clus1 = [1 2 4];
clus2 = [3 5 6];
clus3 = [7 8];

f = weights';
intcon = length(t);

A = ;


% A = [1,1,1];
% b = 7;
Aeq = [4,2,1];
beq = 12;
lb = zeros(3,1);
ub = [Inf;Inf;1]; % enforces x(3) is binary