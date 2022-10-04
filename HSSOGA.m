% V2.0

% A Hybrid Genetic Algorithm and Sperm Swarm Optimization (HGASSO)

%Shehadeh, H. A., Mustafa, H. M., & Tubishat, M. (2022). A Hybrid Genetic Algorithm and Sperm Swarm Optimization (HGASSO) for Multimodal Functions. International Journal of Applied Metaheuristic Computing (IJAMC), 13(1), 1-33.?

  
clc;clear;close all
%% Problem Definition
gen =1;
for inde = 1:gen
CostFunction=@(x) Sphere(x);     % Cost Function

nVar=30;            % Number of Decision Variables

VarSize=[1 nVar];   % Decision Variables Matrix Size

VarMin=-600;         % Lower Bound of Variables
VarMax=600;         % Upper Bound of Variables

MaxIt=1000;         % Maximum Number of Iterations

nPop=100;           % Population Size (Swarm Size)

VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;


% GA Parameters
%--------------------

pc=0.7;                 % Crossover Percentage
nc=2*round(pc*nPop/2);  % Number of Offsprings (also Parnets)
gamma=0.4;              % Extra Range Factor for Crossover

pm=0.3;                 % Mutation Percentage
nm=round(pm*nPop);      % Number of Mutants
mu=0.1;                 % Mutation Rate
beta=8;

empty_individual.Position=[];
empty_individual.Cost=[];

pop=repmat(empty_individual,nPop,1);
GlobalBest.Cost=inf;

empty_sperm.Position=[];
empty_sperm.Cost=[];
empty_sperm.Velocity=[];
empty_sperm.Best.Position=[];
empty_sperm.Best.Cost=[];

%% Initialization
for i=1:nPop
    
    % Initialize Position
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    
    % Evaluation
    pop(i).Cost=CostFunction(pop(i).Position);
    
    sperm(i).Position=pop(i).Position;
    
    % Initialize Velocity
    sperm(i).Velocity=zeros(VarSize);
    
    % Evaluation
    sperm(i).Cost=pop(i).Cost;
    
    % Update Personal Best
    sperm(i).Best.Position=sperm(i).Position;
    sperm(i).Best.Cost=sperm(i).Cost;
    
    % Update Global Best
    if sperm(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=sperm(i).Best;
        
    end
end

Costs=[pop.Cost];
[Costs, SortOrder]=sort(Costs);
pop=pop(SortOrder);

% Store Best Solution
BestSol=pop(1);

% Array to Hold Best Cost Values
BestCost=zeros(MaxIt,1);

% Store Cost
WorstCost=pop(end).Cost;
%% Main Loop
for it=1:MaxIt
    
    P=exp(-beta*Costs/WorstCost);
    P=P/sum(P);
    
    popc=repmat(empty_individual,nc/2,2);
    for k=1:nc/2
        
        % Select Parents Indices (Roulette Wheel Selection)
        i1=RouletteWheelSelection(P);
        i2=RouletteWheelSelection(P);
        p1=pop(i1);
        p2=pop(i2);
        
        % Apply Crossover
        [popc(k,1).Position, popc(k,2).Position]=Crossover(p1.Position,p2.Position,gamma,VarMin,VarMax);
        
        % Evaluate Offsprings
        popc(k,1).Cost=CostFunction(popc(k,1).Position);
        popc(k,2).Cost=CostFunction(popc(k,2).Position);
        
    end
    popc=popc(:);
    
    % Mutation
    popm=repmat(empty_individual,nm,1);
    for k=1:nm
        
        % Select Parent
        i=randi([1 nPop]);
        p=pop(i);
        
        % Apply Mutation
        popm(k).Position=Mutate(p.Position,mu,VarMin,VarMax);
        
        % Evaluate Mutant
        popm(k).Cost=CostFunction(popm(k).Position);
        
    end
    
    % Create Merged Population
    pop=[pop
        popc
        popm]; %#ok
    
    % Sort Population
    Costs=[pop.Cost];
    [Costs, SortOrder]=sort(Costs);
    pop=pop(SortOrder);
    
    % Update Worst Cost
    WorstCost=max(WorstCost,pop(end).Cost);
    
    % Truncation
    pop=pop(1:nPop);
    Costs=Costs(1:nPop);
    
    % Store Best Solution Ever Found
    BestSol=pop(1);
    
    
    
    % Store Best Cost Ever Found
    
    for i=1:nPop
        if pop(i).Cost<=sperm(i).Cost
            sperm(i).Position = pop(i).Position;
            sperm(i).Cost = pop(i).Cost;
        end
        Cx(i) = sperm(i).Cost;
    end
    [BestCost(it),r]=min(Cx);
    GlobalBest.Cost=sperm(r).Cost;
    GlobalBest.Position=sperm(r).Position;
    
    BstCostGA(it)=BestCost(it);
   
   %% SSO Velocity and Position Update 
    for i=1:nPop

         sperm(i).Velocity=rand()*(log10((7-14)*rand(1,1)+7))*sperm(i).Velocity...
             +(log10((7-14)*rand(1,1)+7))* (log10((35.5-38.5)*rand(1,1)+35.5))*(sperm(i).Best.Position-sperm(i).Velocity)...
             +(log10((7-14)*rand(1,1)+7))* (log10((35.5-38.5)*rand(1,1)+35.5))*(GlobalBest.Position-sperm(i).Position);
         sperm(i).Velocity
        % Apply Velocity Limits
        sperm(i).Velocity = max(sperm(i).Velocity,VelMin);
        sperm(i).Velocity = min(sperm(i).Velocity,VelMax);
        
        % Update Position
        sperm(i).Position = sperm(i).Position + sperm(i).Velocity;
        
        % Velocity Mirror Effect
        IsOutside=(sperm(i).Position<VarMin | sperm(i).Position>VarMax);
        sperm(i).Velocity(IsOutside)=-sperm(i).Velocity(IsOutside);
        
        % Apply Position Limits
        sperm(i).Position = max(sperm(i).Position,VarMin);
        sperm(i).Position = min(sperm(i).Position,VarMax);
  
        % Evaluation
        sperm(i).Cost = CostFunction(sperm(i).Position);
        
        % Update Personal Best
        if sperm(i).Cost<sperm(i).Best.Cost
            
            sperm(i).Best.Position=sperm(i).Position;
            sperm(i).Best.Cost=sperm(i).Cost;
            
            % Update Global Best
            if sperm(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=sperm(i).Best;
                
            end
            
        end
        
    end
    
    for i=1:nPop
        if sperm(i).Cost<=pop(i).Cost
            pop(i).Position = sperm(i).Position;
            pop(i).Cost = sperm(i).Cost;
            
        end
        Cx(i)=pop(i).Cost;
    end
    
    [BestCost(it),r]=min(Cx);
    GlobalBest.Cost=pop(r).Cost;
    GlobalBest.Position=pop(r).Position;
    
    BestCostHSSOGA(it,inde)=BestCost(it);
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
   
end

BestCostHSSOGA;
BestSol = GlobalBest;
GlobalBest;
GlobalBestHSSOGA (inde,1) = GlobalBest.Cost;
% BestSol.Position
for i=1:nVar
    results(i,inde)=BestSol.Position(i);
end
results;

inde= inde+1;
end

% BestCostHSSOGA
% GlobalBestHSSOGA
% results

gbcHSSOGA = 0;
for i = 1:gen
    gbcHSSOGA = gbcHSSOGA+GlobalBestHSSOGA(i);
    
end
gbcHSSOGAave = gbcHSSOGA/gen;
%gbcHSSOGAave
allTable(1,3) = gbcHSSOGAave;
    
for i = 1:1000
    bcHSSOGA = 0;
    for j =1:gen
        bcHSSOGA = bcHSSOGA + BestCostHSSOGA(i,j);
    end
    bcHSSOGAave(i,1) = bcHSSOGA/gen;
    allTable(i,1) = bcHSSOGA/gen;
end
% bcHSSOGAave
    
for i = 1:30
    gbHSSOGA = 0;
    for j =1:gen
        gbHSSOGA = gbHSSOGA + results(i,j);
    end
    gbHSSOGAave(i,1) = gbHSSOGA/gen;
    allTable(i,5) = gbHSSOGA/gen;
end
%gbHSSOGAave

figure;
%plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;