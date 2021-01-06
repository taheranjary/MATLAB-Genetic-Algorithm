clear all;
clc;

% INITIALIZE
POPULATION_SIZE = 2000;
MAX_GEN = 1000;

pm = 0.1;   % probability of mutation
pc = 0.9;   % probability of crossover
alpha = 2;  % blend crossover hyperparameter

candidates = 100*rand(POPULATION_SIZE,2);
%candidates = [3 3; 2 2; 9 9; 5 5];

%% BEGIN EVOLUTION
tic;
best = zeros(MAX_GEN,2);
worst = zeros(MAX_GEN,2);

for g = 1:1:MAX_GEN

    fitnesses = fitness(candidates);

    %% Rank candidates according to fitness score
    [out,idx] = sort(fitnesses, 'ascend');
    fitnesses = fitnesses(idx);
    candidates = candidates(idx,:);
    
    % track best and worst candidate in each gen
    best(g,:) = candidates(1,:);
    worst(g,:) = candidates(POPULATION_SIZE,:);

    %% SELECTION - Rank Based Selection
    weights = 1./(fitnesses+1);
    weights = weights./(sum(weights));

    selected_candidates = zeros(POPULATION_SIZE,2);
    ids = 1:1:POPULATION_SIZE;
    for i = 1:1:POPULATION_SIZE
        ix = randsample(ids,1,true,weights);
        selected_candidates(i,:) = candidates(ix,:);
    end

    %% CROSSOVER - Blend Crossover
    offsprings = zeros(POPULATION_SIZE,2);

    for i = 1:2:POPULATION_SIZE
        parent_1 = selected_candidates(i,:);
        parent_2 = selected_candidates(i+1,:);

        if rand() < pc
            % blend crossover
           [child_1, child_2] = blendCrossover(alpha, parent_1, parent_2);

        else
            child_1 = parent_1;
            child_2 = parent_2;
        end
        offsprings(i,:) = child_1;
        offsprings(i+1,:) = child_2;

    end

    %% MUTATION
    for i = 1:1:POPULATION_SIZE
        if rand() < pm
            % mutate
            gi = randi(2);
            u = offsprings(i,gi);
            sd = 0.5;
            gene = sd*randn() + u;
            offsprings(i,gi) = gene; 
            %offsprings(i,:) = mutate(offsprings(i,:));
        end
    end

    %% CYCLE
    candidates = offsprings;

end
toc;
%% ANALYSIS
figure();
hold on;
plot(worst(:,1));
plot(ones(MAX_GEN));
title('gene 1 of worst candidate in each generation');
xlabel('generation');
ylabel('gene value');
%ylabel(['average ' num2str(mean(worst(:,1)))]);
legend('gene 1 of worst candidate','global optimum gene 1');
hold off;
%disp(mean(worst(:,1)));

figure();
hold on;
plot(worst(:,2));
plot(ones(MAX_GEN));
title('gene 2 of worst candidate in each generation');
xlabel('generation');
ylabel('gene value');
%ylabel(['average ' num2str(mean(worst(:,2)))]);
legend('gene 2 of worst candidate','global optimum gene 2');
hold off;
%disp(mean(worst(:,2)));

figure();
hold on;
plot(best(:,1));
plot(ones(MAX_GEN));
title('gene 1 of best candidate in each generation');
xlabel('generation');
ylabel('gene value');
%ylabel(['average ' num2str(mean(best(:,1)))]);
legend('gene 1 of best candidate','global optimum gene 1');
hold off;
%disp(mean(best(:,1)));

figure();
hold on;
plot(best(:,2));
plot(ones(MAX_GEN));
title('gene 2 of best candidate in each generation');
xlabel('generation');
ylabel('gene value');
%ylabel(['average ' num2str(mean(best(:,2)))]);
legend('gene 2 of best candidate','global optimum gene 2');
%disp(mean(best(:,2)));

disp(['last best candidate solution: (' num2str(best(end,1)) ', ' num2str(best(end,2)) ')']);
disp(['average worst candidate: (' num2str(mean(worst(:,1))) ', ' num2str(mean(worst(:,2))) ')']);
disp(['average best candidate: (' num2str(mean(best(:,1))) ', ' num2str(mean(best(:,2))) ')']);


%% Function Defs
function y = fitness(x)
    y = arrayfun(@rosenbrock,x(:,1),x(:,2));
end

function y = rosenbrock(x1,x2)
    y = 100*((x2-(x1^2))^2) + (1-x1)^2;
end

% Blend Crossover
function [child_1, child_2] = blendCrossover(alpha, parent_1, parent_2)
    child_1 = zeros(1,2);
    child_2 = zeros(1,2);

    % blending first gene only
    edge_1 = parent_1(1) - alpha*( parent_2(1) - parent_1(1) );
    edge_2 = parent_2(1) + alpha*( parent_2(1) - parent_1(1) );
    a = min( edge_1, edge_2 );
    b = max( edge_1, edge_2 );
    
    g1 = (b-a)*rand() + a;
    g2 = (b-a)*rand() + a;
    
    child_1(1) = g1;
    child_2(1) = g2;
    
    % blending second gene only
    edge_1 = parent_1(2) - alpha*( parent_2(2) - parent_1(2) );
    edge_2 = parent_2(2) + alpha*( parent_2(2) - parent_1(2) );
    a = min( edge_1, edge_2 );
    b = max( edge_1, edge_2 );
    
    g1 = (b-a)*rand() + a;
    g2 = (b-a)*rand() + a;
    
    child_1(2) = g1;
    child_2(2) = g2;   
    
end