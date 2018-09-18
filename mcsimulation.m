function [mu,sigma,MU] = mcsimulation(model,inputU,seed)
% Simulation of the dynamic GP model, where the output variance is
% propagated using the Monte Carlo method
%
%% Syntax
%  [mu, sigma] = mcsimulation(model,inputU,seed)
% 
%% Description
% Idea: at every time step the output of GP model is approximated with
% Nsamples samples, which are used as the future inputs of the GP model.
% Samples are re-used if necessary (ie. y(k-1) for y(k-2) if lag=2 etc.) 
% Uses routines gpx and gmx_sample. 
% 
% Input:
% * hyp      ... the structure of optimized hyperparameters 
% * inf      ... the function specifying the inference method 
% * meanfunc ... the prior mean function
% * cov      ... the specified covariance function, see help covFun for more info 
% * lik      ... the likelihood function
% * input    ... the input part of the training data,  NxD matrix
% * target   ... the output part of the training data (ie. target), Nx1 vector 
% * test     ... the input matrix for simulation, kxD vector, see
%                construct.m for more info  
% * lag      ... the order of the model (number of used lagged outputs) 
% * Nsamples ... the number of samples used in algorithm (ie. runs of simulation) 
% 
% Output:
% * mu    ... the mean predicted output 
% * s2    ... the associated variances (with noise variances)
% * MU    ... the matrix of all predicted means, kxNsamples
% * SIG2  ... the associated predicted variances 
% 
% 
%% 
% * Written by Wilmer Ariza November 2016
% * Based in the code Written by J. Prikryl, November 2010
% * Based on the work of C.E. Rasmussen, A. Girard, K. Azman. 
%

% meanfunc ... 'mean' is used as a matlab core function in this file
Nsamples=50;lag=1;
% Preallocate mu and sigma
[m,n]=size(seed);
mu = zeros ( m, n );
sigma = zeros ( m, n );

% 1st step - input is a point
AUX=seed;
X_prime{1}=1;
for j=2:n+1
    X_prime{j}=[AUX(1,j-1) inputU(1,end-1) inputU(1,end)];
end
[y_p, s2] = multigpPosteriorMeanVar(model, X_prime,'false');
for j=2:n+1
    A(1,j-1)=cell2mat(y_p(j).');
    B(1,j-1)=cell2mat(s2(j).');
end
mu(1,:)=A;
sigma(1,:)=B;
%% First Point montecarlos exploration

for i=1:n
  MU(1,:,i) = mu(1,i)*ones(1,Nsamples); 
  SIG2(1,:,i) = sigma(1,i)*ones(1,Nsamples);
end



for i=1:n
  pdf = gmx_sample(mu(1,i),sigma(1,i),Nsamples);
  PDF(:,lag,i) = pdf;
end

test_ = repmat( inputU(2,:), Nsamples, 1 );
X_prime{1}=1;
for j=2:n+1
    X_prime{j}=[PDF(:,lag,j-1) test_(:,end-1) test_(:,end)];
end
[MUX, SIG2X] = multigpPosteriorMeanVar(model, X_prime,'false');
A=zeros(Nsamples,n);B=zeros(Nsamples,n);


for j=2:n+1
    A(:,j-1)=cell2mat(MUX(j).');
    B(:,j-1)=cell2mat(SIG2X(j).');
end
MUX=A;
SIG2X=B;



mu(2,:)=mean(MUX);
for i=1:n
    sigma(2,i) = mean(SIG2X(:,i)) + mean((MUX(:,i)-mu(2,i)).^2);
end

for i=1:n
  MU(2,:,i) = MUX(:,i);
  SIG2(2,:,i) = SIG2X(:,i);
end

% steps from 3 on ...
for k=3:m

    for i=1:n
        pdf = gmx_sample(MU(k-1,:,i),SIG2(k-1,:,i),Nsamples);
        PDF(:,lag,i) = pdf;
    end
    test_ = repmat( inputU(k,:), Nsamples, 1 );
    X_prime{1}=1;
    for j=2:n+1
        X_prime{j}=[PDF(:,lag,j-1) test_(:,end-1) test_(:,end)];
    end
    [MUX, SIG2X] = multigpPosteriorMeanVar(model, X_prime,'false');
    A=zeros(Nsamples,n);B=zeros(Nsamples,n);
    for j=2:n+1
        A(:,j-1)=cell2mat(MUX(j).');
        B(:,j-1)=cell2mat(SIG2X(j).');
    end
    MUX=A;
    SIG2X=B;

    mu(k,:)=mean(MUX);
    for i=1:n
        sigma(k,i) = mean(SIG2X(:,i)) + mean((MUX(:,i)-mu(k,i)).^2);
    end
    
    for i=1:n
        MU(k,:,i) = MUX(:,i);
        SIG2(k,:,i) = SIG2X(:,i);
    end

end

return