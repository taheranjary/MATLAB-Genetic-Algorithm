clear all;
clc;

rng default

rosenbrock = @(x)( 100*(x(:,2) - x(:,1).^2).^2 + (1 - x(:,1)).^2 );

%% Genetic Algorithm
fprintf('\nGA \n ');
tic;
[x,fval,exitflag,output,population,scores] = ga(rosenbrock,2)
toc;
fprintf('\nGA options \n ');
[user,sys] = memory%place in the begging of your program
tic;
options = optimoptions('ga','MaxGenerations',1000,'MaxStallGenerations',Inf,'InitialPopulationMatrix',100*rand(2000,2),'PopulationSize',2000,'PlotFcn',@gaplotbestf);
[x,fval,exitflag,output,population,scores] = ga(rosenbrock,2,options)
toc;
[user2,sys2] = memory%place in the end of your program
memory_used_in_bytes=user2.MemAvailableAllArrays-user.MemAvailableAllArrays

%% Simulated Annealing
fprintf('\nSimulated Annealing \n ');
x0 = [10 5];
tic;
[x,fval,exitflag,output] = simulannealbnd(rosenbrock,x0)
toc;
%% Particle Swarm
fprintf('\nParticle Swarm \n ');
tic;
[x,fval,exitflag,output] = particleswarm(rosenbrock,2)
toc;
%% Pattern Search
fprintf('\nPattern Search \n ');
x0 = [5, 5];
tic;
[x,fval,exitflag,output] = patternsearch(rosenbrock,x0)
toc;
%% Function Defs
function memAnalysis()
    [user,sys] = memory;
    
    disp('user: ');
    disp(['MaxPossibleArrayBytes: ' num2str(user.MaxPossibleArrayBytes)]);
    disp(['MemAvailableAllArrays: ' num2str(user.MemAvailableAllArrays)]);
    fprintf(['MemUsedMATLAB: ' num2str(user.MemUsedMATLAB) '\n']);
    
    disp('sys: ');
    disp(['VirtualAddressSpace Available: ' num2str(sys.VirtualAddressSpace.Available)]);
    disp(['VirtualAddressSpace Total: ' num2str(sys.VirtualAddressSpace.Total)]);
    disp(['SystemMemory Available: ' num2str(sys.SystemMemory.Available)]);
    disp(['PhysicalMemory Available: ' num2str(sys.PhysicalMemory.Available)]);
    disp(['PhysicalMemory Total: ' num2str(sys.PhysicalMemory.Total)]);
end
% function y = rosenbrock(x)
%     x1 = x(1);
%     x2 = x(2);
%     y = 100*((x2-(x1^2))^2) + (1-x1)^2;
% end

% var = 1:4;
% W = [0.5, 1/6, 1/6, 1/6];
% X = randsample(var,1,true,W)