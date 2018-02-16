%% 3 A Naive approach - Random Walk
clear
close all
d = 2;
N = 5000; %Number of particles
walk = zeros(N, 2);
n = [-1,0; 0,1; 1,0; 0,-1];

for i = 3 : N + 1
    newDir = datasample(n,1);
    next = newDir + walk(i-2, :);
    
    if ismember(next, walk, 'rows') == ones(1, 2)
        walk(i-1,:) = next;
        break
    end
    
    walk(i - 1, :) = next;
end
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-1,1),walk(i-1,2),'-*r')
walk(1:i-1,:);

%% Estimating Cn(2)
N = 5000; %Number of particles
walk = zeros(N, 2);
n = [-1,0; 0,1; 1,0; 0,-1];
%c = zeros(10, 1);
nsa = zeros(10,1);

for NbrOfSteps = 1:10
    for NbrOfTries = 1:5000
        for i = 3:NbrOfSteps+2
            newDir = datasample(n,1);
            if i == NbrOfSteps+2
                nsa(NbrOfSteps) = nsa(NbrOfSteps) + 1;
                break
            end
            next = newDir + walk(i-2, :);
            if ismember(next, walk, 'rows') == ones(1, 2)
                walk(i-1,:) = next;
                break
            end
            walk(i - 1, :) = next;
        end
    end
    %c(NbrOfSteps)=(nsa(NbrOfSteps)./N).*(4.^NbrOfSteps);
end
k = 1:10;
k1 = 4.^k;
format long;
c1 = (nsa./N).*k1';

%% 4 Improving - actual self avoiding walk
close all
d = 2;
N = 1e3;
walk = zeros(N, 2);
omega_0 = zeros(N, 1);
g0 = zeros(N,1);
n = [-1,0; 0,1; 1,0; 0,-1];

for i = 3 : N + 1
    % n is the neighbour-vector
    neigh(1,:) = walk(i-2,:) + [-1,0];
    neigh(2,:) = walk(i-2,:) + [0,1];
    neigh(3,:) = walk(i-2,:) + [1,0];
    neigh(4,:) = walk(i-2,:) + [0,-1];
    if ismember(neigh,walk,'rows') == ones(4,1)
        break
    end
    g0(i-2) = 1-(sum(ismember(neigh,walk,'rows'))/4);
    newDir = datasample(n,1);
    while ismember(newDir + walk(i-2,:),walk,'rows') == ones(1,2)
        newDir = datasample(n,1);
    end
    next = newDir + walk(i-2, :);
    walk(i-1, :) = next;
end

%% Estimating Cn(2) using SIS drawing from a gn = SAW
d = 2;
N = 1e3;
nbrOfSteps = 10;
neigh = zeros(4, 2);
walk = zeros(N, 2);
omega = zeros(N, nbrOfSteps);
g = zeros(N,1);
n = [-1,0; 0,1; 1,0; 0,-1];
omega_0 = 1;

for nbrOfSteps = 1:nbrOfSteps
    for nbrOfTries = 1:N
        neigh = zeros(4, 2);
        walk = zeros(N, 2);
        for j = 3:nbrOfSteps+2
            neigh(1,:) = walk(j-2,:) + [-1,0];
            neigh(2,:) = walk(j-2,:) + [0,1];
            neigh(3,:) = walk(j-2,:) + [1,0];
            neigh(4,:) = walk(j-2,:) + [0,-1];
            
            if ismember(neigh,walk,'rows') == ones(4,1)
                break
            end
            
            newDir = datasample(n,1);
            
            while ismember(newDir + walk(j-2,:),walk,'rows') == ones(1,2)
                newDir = datasample(n,1);
            end
            
            g(nbrOfTries) = 4-sum(ismember(neigh, walk, 'rows'));
            next = newDir + walk(j-2, :);
            walk(j-1, :) = next;
            
            if nbrOfSteps == 1
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries)*omega_0;
            else
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries).*omega(nbrOfTries, nbrOfSteps-1);
            end
        end
    end
end
c2 = mean(omega);

%% Estimating Cn(2) using SISR drawing from a gn = SAW
d = 2;
d2 = 2*d;
N = 1e3;
nbrOfSteps = 5;
neigh = zeros(d2, d);
walk = zeros(N, d);
omega = zeros(N, nbrOfSteps);
g = zeros(N,1);
n = [eye(d) ; -1*eye(d)];
omega_0 = 1;

for nbrOfSteps = 1:nbrOfSteps
    for nbrOfTries = 1:N
        neigh = zeros(d2, d);
        walk = zeros(N, d);
        for j = 3:nbrOfSteps+2
            for i = 1:d2
                neigh(i, :) = walk(j-2, :) + n(i, :);
            end
            
            if ismember(neigh,walk,'rows') == ones(d2,1)
                break
            end
            
            newDir = datasample(n,1);
            
            while ismember(newDir + walk(j-2,:),walk,'rows') == ones(1,d)
                newDir = datasample(n,1);
            end
            
            g(nbrOfTries) = d2-sum(ismember(neigh, walk, 'rows'));
            next = newDir + walk(j-2, :);
            walk(j-1, :) = next;
            
            if nbrOfSteps == 1
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries)*omega_0;
            else
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries).*omega(nbrOfTries, nbrOfSteps-1);
            end
        end
    end
    %Selection
    CW=cumsum([0 omega(:, nbrOfSteps)']);
    [~,ind] = histc(rand(1,N),CW/CW(end));
    
    %Mutation 
    walk = walk(ind, :);
end
c2 = mean(omega);

A = [1 1 0; 1 2 log(2); 1 3 log(3)];
B = log(c2);
C = B(1:3)'\A
A2 = exp(C(1))
mu2 = exp(C(2))
gam2 = C(3) + 1

%% Estimating Cn(2) using SISR drawing from a gn = SAW for general d
% Man behöver alla partiklarns senaste position samt alla tidigare positioner
% Iterera genom alla partiklar och kolla deras position för att uppdatera
% de senaste.
d = 2;
d2 = 2*d;
N = 1e3;
nbrOfSteps = 5;
neigh = zeros(d2, d);
walk = zeros(N, d);
omega = zeros(N, nbrOfSteps);
g = zeros(N,1);
n = [eye(d); -1*eye(d)];
omega_0 = 1;

for nbrOfSteps = 1:nbrOfSteps
    for nbrOfTries = 1:N
        neigh = zeros(d2, d);
        walk = zeros(N, d);
        for j = 3:nbrOfSteps+2
            
            for i = 1:d2
            neigh(i, :) = walk(j-2, :) + n(i, :);
            end
            
            if ismember(neigh,walk,'rows') == ones(d2,1)
                break
            end
            
            newDir = datasample(n,1);
            
            while ismember(newDir + walk(j-2,:),walk,'rows') == ones(1,d)
                newDir = datasample(n,1);
            end
            
            g(nbrOfTries) = (d2)-sum(ismember(neigh, walk, 'rows'));
            next = newDir + walk(j-2, :);
            walk(j-1, :) = next;
            
            if nbrOfSteps == 1
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries)*omega_0;
            else
                omega(nbrOfTries, nbrOfSteps) = g(nbrOfTries).*omega(nbrOfTries, nbrOfSteps-1);
            end
        end
    end
    %Selection
    CW=cumsum([0 omega(:, nbrOfSteps)']);
    [~,ind] = histc(rand(1,N),CW/CW(end));
    
    %Mutation 
    walk = walk(ind, :);
    
end
c2 = mean(omega);

%%
figure(1)
plot(walk(1:i-2, 1), walk(1:i-2, 2),'-*')
hold on
plot(walk(i-2, 1), walk(i-2,2), '*r')
figure(2)
plot(1:i,g0(1:i),'*')