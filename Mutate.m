%Mutate

function y = Mutate(x, mu, VarMin, VarMax)

    nVar = numel(x);
    
    nmu = ceil(mu*nVar);
    
    for i=1:nmu
        j = randsample(nVar, 1);
    
        sigma = 0.1*(VarMax-VarMin);
    
        y = x;
    
        test = randn(size(j));
        y(j) = x(j)+sigma*test;

        y = max(y, VarMin);
        y = min(y, VarMax);
    end
end