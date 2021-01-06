clear all;
clc;

rng default;
rosenbrock = @(x)( 100*(x(:,2) - x(:,1).^2).^2 + (1 - x(:,1)).^2 );

%profile -memory on

%memories = zeros(6,8,'double');
%[user,sys] = memory;
%memories(1,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(1,:) = memAnalysis();

%% Genetic Algorithm
%fprintf('\nGA \n ');
x = ga(rosenbrock,2);
%disp(x);

%[user,sys] = memory;
%memories(2,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(2,:) = memAnalysis();

%fprintf('\nGA options \n ');
options = optimoptions('ga','MaxGenerations',1000,'MaxStallGenerations',Inf,'InitialPopulationMatrix',100*rand(2000,2),'PopulationSize',2000,'PlotFcn',@gaplotbestf);
x = ga(rosenbrock,2,options);
%disp(x);

%[user,sys] = memory;
%memories(3,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(3,:) = memAnalysis();

%% Simulated Annealing
%fprintf('\nSimulated Annealing \n ');
x0 = [10 5];
x = simulannealbnd(rosenbrock,x0);
%disp(x);

%[user,sys] = memory;
%memories(4,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(4,:) = memAnalysis();

%% Particle Swarm
%fprintf('\nParticle Swarm \n ');
%profile clear;
x = particleswarm(rosenbrock,2);
%profile report;
%disp(x);

%[user,sys] = memory;
%memories(5,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(5,:) = memAnalysis();

%% Pattern Search
%fprintf('\nPattern Search \n ');
x0 = [5, 5];
x = patternsearch(rosenbrock,x0);
%disp(x);

%[user,sys] = memory;
%memories(6,:) = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
%memories(6,:) = memAnalysis();

%%
%disp(memories);

%% Function Defs
% function m = memAnalysis()
%     [user,sys] = memory;
%     
%     disp('user: ');
%     disp(['MaxPossibleArrayBytes: ' num2str(user.MaxPossibleArrayBytes)]);
%     disp(['MemAvailableAllArrays: ' num2str(user.MemAvailableAllArrays)]);
%     disp(['MemUsedMATLAB: ' num2str(user.MemUsedMATLAB)]);
%     
%     fprintf('\nsys: \n');
%     disp(['VirtualAddressSpace Available: ' num2str(sys.VirtualAddressSpace.Available)]);
%     disp(['VirtualAddressSpace Total: ' num2str(sys.VirtualAddressSpace.Total)]);
%     disp(['SystemMemory Available: ' num2str(sys.SystemMemory.Available)]);
%     disp(['PhysicalMemory Available: ' num2str(sys.PhysicalMemory.Available)]);
%     disp(['PhysicalMemory Total: ' num2str(sys.PhysicalMemory.Total)]);
%     
%     m = [user.MaxPossibleArrayBytes, user.MemAvailableAllArrays, user.MemUsedMATLAB, sys.VirtualAddressSpace.Available, sys.VirtualAddressSpace.Total, sys.SystemMemory.Available, sys.PhysicalMemory.Available, sys.PhysicalMemory.Total];
% end