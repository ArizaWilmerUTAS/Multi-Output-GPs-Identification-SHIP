% Process data obtained from the simulation of a ship to obtain a GP
% Dynamic System
%
%
%  Written by Wilmer Ariza based in the work of Jus Kocijan ,Dr Alvarez
%  and,Dr Lawrence
%% Load data
clc
clear
% load inputYk_1                                       % load test data
load outputY                                      % load test data  outputs
load inputU                                     % load test data control signals
figure_c=1;
%% Training

% Build training data (delayed outputs y first)
% Output=Output(:,1:100);
% inputU=inputU(:,1:100);
[m,n]=size(Output);
Output=Output';
inputY=[zeros(1,m);Output];
inputY(end,:)=[];
inputYT=inputY(1:round(n/2),:);% input Y training data
inputYV=inputY(round(n/2)+1:end,:);%input Y validation data
OutputT=Output(1:round(n/2),:);%output training data
OutputV=Output(round(n/2)+1:end,:);%output validation data
inputU=inputU';
[m,n]=size(inputU');
inputU=[zeros(1,m);inputU];
inputU(end,:)=[];
inputUT=inputU(1:round(n/2),:);%input control training data
inputUV=inputU(round(n/2)+1:end,:);%input control validation data
% inputUT=inputU;
% OutputT=Output;
% inputYT=inputY;
[m,n]=size(inputUT);
Tx=Output(end,1);% TIme data

input = [inputYT(:,2:end) inputUT(:,2:3)] ; %assemble data for normzalization
target = OutputT(:,2:end);


% normalize training data
[nInput,inputMin,inputMax, nTarget, targetMin, targetMax] = preNorm(input,target);

dataSetName = 'WILMERreal';
%

experimentNo = 1;
[~,n]=size(nInput);
for j=1:n-2
    if j==1
        XTemp{j}=[ nInput(:,j) nInput(:,end-1) nInput(:,end) ];% Input Y Input URPM
    else
        XTemp{j}=[ nInput(:,j) nInput(:,end-1) nInput(:,end) ];
    end
    yTemp{j}=[nTarget(:,j)];
end


options = multigpOptions('ftc');
options.kernType = 'ggwhite';
options.optimiser = 'optimiMinimize';%'optimiMinimize','scg'
options.nlf = 1;

X = cell(size(yTemp, 2)+options.nlf,1);
y = cell(size(yTemp, 2)+options.nlf,1);

for j=1:options.nlf
    y{j} = [];
    X{j} = [1 1 1 ];
end
for i = 1:size(yTemp, 2)
    y{i+options.nlf} = yTemp{i};
    X{i+options.nlf} = XTemp{i};
end

q = 3;
d = size(yTemp, 2) + options.nlf;

%%  Creates the model
model = multigpCreate(q, d, X, y, options);

params = modelExtractParam(model);
index = paramNameRegularExpressionLookup(model, 'multi .* inverse .*');
params(index) = log(10 +  5*rand(1,length(index)));
index = paramNameRegularExpressionLookup(model, 'multi .* inverse width output .*');
params(index) = log(10 +  2*rand(1,length(index)));
index = paramNameRegularExpressionLookup(model, 'multi .* variance');
params(index) = log(10 +  2*rand(1,length(index)));
index = paramNameRegularExpressionLookup(model, 'multi .* sensitivity');
params(index) = 1+ rand(1, length(index))
%params=[-0.705718051131773,0,6.33045758968866,1.47365006311238,3.59324131462178,1.44638764331041,3.76420250788585,1.44912791163510,3.61995164811514,1.44408347908599,-12.7681380068096,-6.96131676228187,-6.85149344665989,-6.39042235784673];
model = modelExpandParam(model, params);

display = 1;
iters = 300;
%% Trains the model
init_time = cputime;
model = multigpOptimise(model, display, iters);
elapsed_time = cputime - init_time;
params = modelExtractParam(model)

%% Save the results.
capName = dataSetName;
capName(1) = upper(capName(1));
save(['dem' capName num2str(experimentNo) '.mat'], 'model');
%% Predict values
Xt{1}=[1 1 1 ];
for j=2:n-1
    if j==2
        Xt{j}=[nInput(:,j-1) nInput(:,end-1) nInput(:,end)];
    else
        Xt{j}=[nInput(:,j-1) nInput(:,end-1) nInput(:,end)];
    end
end

[mu, s2naive] = multigpPosteriorMeanVar(model, Xt,'false');
for j=2:n-1
    A(:,j-1)=cell2mat(mu(j).');
    B(:,j-1)=cell2mat(s2naive(j).');
end
mu=A;
s2naive=B;
%y = postNorm(mu, targetMin, targetMax);

%%
[m,n]=size(inputUT);
T=linspace(0,Tx/2,m)';
names=cell(3,1);
names={'Udot m2/s Training';'V m/s Training';'r rad/s Training';'p rad/s Training';'p rad/s Training'};
% figure (1)
Yprime = postNorm(mu, targetMin, targetMax);
Ys = postNorm(nTarget, targetMin, targetMax);
figure_c=PLOTGPSHIP(T,Ys,Yprime,s2naive,names,figure_c);

%% Validatio of model
inputv = [inputYV(:,2:end) inputUV(:,2:3)] ;
[m,n]=size(inputv);
% normalize training data
[nInputv] = preNorm(inputv,inputMin,inputMax);
A=zeros(m,n-2);B=zeros(m,n-2);
Xt{1}=[1 1 1 ];
for j=2:n-1
    if j==2
        Xt{j}=[nInputv(:,j-1) nInputv(:,end-1) nInputv(:,end)];
    else
        Xt{j}=[nInputv(:,j-1) nInputv(:,end-1) nInputv(:,end)];
    end
end
%Xt ={1;nInput(:,1);nInput(:,2);nInput(:,3);nInput(:,4);nInput(:,5)};
%Xt = [-1.25 1.27];
[muv, s2v] = multigpPosteriorMeanVar(model, Xt,'false');
for j=2:n-1
    A(:,j-1)=cell2mat(muv(j).');
    B(:,j-1)=cell2mat(s2v(j).');
end
muv=A;
s2v=B;
%y = postNorm(mu, targetMin, targetMax);

[m,n]=size(inputUV);
T=linspace(0,Tx/2,m)';
names=cell(3,1);
names={'Udot m2/s Validation';'V m/s Validation';'r rad/s Validation';'p rad/s Validation';'p rad/s Training'};
% figure (1)
targetv = OutputV(:,2:end);
[nTargetv] = preNorm(targetv,targetMin, targetMax);

Yprime = postNorm(muv, targetMin, targetMax);
Ys = postNorm(nTargetv, targetMin, targetMax);
figure_c=PLOTGPSHIP(T,Ys,Yprime,s2v,names,figure_c);


%%

y_start=[nInput;nInputv];

[m,n]=size(y_start);

AUX=[y_start(1,1:end-2)];
A=0;
B=0;
for i=1:m
    
    X_prime{1}=[1 1 1];
    for j=2:n-1
        if j==2
            X_prime{j}=[AUX(1,j-1) y_start(i,end-1) y_start(i,end)];
        else
            X_prime{j}=[AUX(1,j-1) y_start(i,end-1) y_start(i,end)];
            
        end
        
    end
    
    
    
    [y_p, s2] = multigpPosteriorMeanVar(model, X_prime,'false');
    for j=2:n-1
        A(1,j-1)=cell2mat(y_p(j).');
        B(1,j-1)=cell2mat(s2(j).');
    end
    y_prime(i,:)=A;
    s2_prime(i,:)=B;
    if i>=m
    else
        
        AUX=A;
    end
    
end
T=linspace(0,Tx,m)';
% figure (3)
names={'Udot m2/s Simulation';'V m/s Simulation';'r rad/s Simulation';'p rad/s Simulation',;'p rad/s Simulation'};
Yprime = postNorm(y_prime, targetMin, targetMax);
Ys = postNorm([nTarget;nTargetv], targetMin, targetMax);
%%
load inputU;
inputU=inputU';
[m,n]=size(inputU);
T=linspace(0,Tx,m)';
f1=figure(20);
set(gcf,'Color','w');
%
subplot(2,1,1);
plot(T,inputU(:,end-1),'LineWidth',1);

title('U [RPM]');
set(gca,'fontsize',12);
set(gca,'Linewidth', 1);
%
subplot(2,1,2);
plot(T,inputU(:,end),'LineWidth',1);

title('RUDDER ANGLE [rad]');
set(gca,'fontsize',12);
set(gca,'Linewidth', 1);
saveas(f1,'input','epsc');
% %% Neural Network 1
% X1=inputUT(:,2:3)';
% T1=OutputT(:,2:5)';
% X=num2cell(X1,1);
% T1=num2cell(T1,1);
% net = layrecnet(1:2,40);
% 
% %view(net)
% [Xs,Xi,Ai,Ts] = preparets(net,X,T1);
% net.trainParam.epochs=10000;
% net.layers{1}.transferFcn = 'tansig';
% net.trainFcn = 'trainlm';
% net = train(net,Xs,Ts,Xi,Ai);
% Xs=inputU(:,2:3)';
% Xs=num2cell(Xs,1);
% Y = net(Xs,Xi,Ai);
% Ts=Output(:,2:5)';
% Ts=num2cell(Ts,1);
% perf = perform(net,Y,Ts);
% Y=cell2mat(Y);
% ynn1=Y';
% %ynn1=[cell2mat(ydata(7))',cell2mat(ydata(5))',cell2mat(ydata(3))',cell2mat(ydata(1))'];
% %% Neural Network 2
% X1=inputUT(:,2:3)';
% T1=OutputT(:,2:5)';
% X=num2cell(X1,1);
% T1=num2cell(T1,1);
% net = layrecnet(1:3,40);
% 
% %view(net)
% [Xs,Xi,Ai,Ts] = preparets(net,X,T1);
% net.trainParam.epochs=10000;
% net.layers{1}.transferFcn = 'tansig';
% net.trainFcn = 'trainlm';
% net = train(net,Xs,Ts,Xi,Ai);
% Xs=inputU(:,2:3)';
% Xs=num2cell(Xs,1);
% Y = net(Xs,Xi,Ai);
% Ts=Output(:,2:5)';
% Ts=num2cell(Ts,1);
% perf = perform(net,Y,Ts);
% Y=cell2mat(Y);
% ynn2=Y';
% %ynn2=[cell2mat(ydata(7))',cell2mat(ydata(5))',cell2mat(ydata(3))',cell2mat(ydata(1))'];
%%
load ynn1.mat
ynn1=resample(ynn1',1,2);
load ynn2.mat
ynn2=resample(ynn2',1,2);
figure_c=PLOTGPNN(T,Ys,Yprime,ynn1,ynn2,s2_prime,names,figure_c);

%%

Ytotal=Ys;
YGps=Yprime;
Ynn1=ynn1;
Ynn2=ynn2;

pressgp=mean(sum((Ytotal-YGps).^2));
pressnn1=mean(sum((Ytotal-Ynn1).^2));
pressnn2=mean(sum((Ytotal-Ynn2).^2));
[m,n]=size(Ytotal);
% Root Mean Squared Error
RMSEGps=mean(sqrt(mean((Ytotal - YGps).^2)));
RMSEnn1=mean(sqrt(mean((Ytotal - Ynn1).^2)));
RMSEnn2=mean(sqrt(mean((Ytotal - Ynn2).^2)));

fprintf(1, 'Root mean squre error\n')
fprintf(1, '             GPs                 NN1                  NN2\n')    % Column Titles
fprintf(1, '\t\t%2d\t\t%2d\t\t%2d\n', [RMSEGps,RMSEnn1,RMSEnn2])      % Write Rows

fprintf(1, 'Predicted residual error sum of squares\n')
fprintf(1, '             GPs                 NN1                  NN2\n')    % Column Titles
fprintf(1, '\t\t%2d\t\t%2d\t\t%2d\n', [pressgp,pressnn1,pressnn2])      % Write Rows


