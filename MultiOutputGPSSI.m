% Process data obtained from the simulation of a ship to obtain a GP
% Dynamic System
%  
% 
%  Written by Wilmer Ariza
%% Load data
clc
clear 

% load inputYk_1                                       % load test data                                         
load outputY                                       % load test data                                         
load inputU                                       % load test data
figure_c=1;
%% Training

% Build training data (delayed outputs y first) 
% inputY=inputY';
% [m,n]=size(inputY);
% inputYT=inputY(1:round(m/2),:);
% inputYV=inputY(round(m/2)+1:end,:);
[m,n]=size(Output);
Output=Output';
inputY=[zeros(1,m);Output];
inputY(end,:)=[];
inputYT=inputY(1:round(n/2),:);
inputYV=inputY(round(n/2)+1:end,:);
OutputT=Output(1:round(n/2),:);
OutputV=Output(round(n/2)+1:end,:);
inputU=inputU';
inputUT=inputU(1:round(n/2),:);
inputUV=inputU(round(n/2)+1:end,:);
% inputUT=inputU;
% OutputT=Output;
% inputYT=inputY;
[m,n]=size(inputUT);
Tx=Output(end,1);

input = [inputYT(:,2:end) inputUT(:,2:3)] ; 
target = OutputT(:,2:end); 


% normalize training data
[nInput,inputMin,inputMax, nTarget, targetMin, targetMax] = preNorm(input,target);     

dataSetName = 'WILMERreal';
%

experimentNo = 1;
[~,n]=size(nInput);
for j=1:n-2
    if j==1
        XTemp{j}=[nInput(:,j) nInput(:,end-1) nInput(:,end)];
    else
    XTemp{j}=[nInput(:,j) nInput(:,end-1) nInput(:,end)];
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
   X{j} = 1;
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
params(index) = log(100);
%params=[120   0.00001     0     1     0     0     1     0     0     2    -10   -10    -10];
model = modelExpandParam(model, params);

display = 1;
iters = 3000;
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
Xt{1}=1;
for j=2:n-1
    if j==2
    Xt{j}=[nInput(:,j-1) nInput(:,end-1) nInput(:,end)];
    else
        Xt{j}=[nInput(:,j-1) nInput(:,end-1) nInput(:,end)];
    end
end
%Xt ={1;nInput(:,1);nInput(:,2);nInput(:,3);nInput(:,4);nInput(:,5)};
%Xt = [-1.25 1.27];
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
figure_c=PLOTGP(T,Ys,Yprime,s2naive,names,figure_c);

%% Validatio of model
inputv = [inputYV(:,2:end) inputUV(:,2:3)] ; 
[m,n]=size(inputv);
% normalize training data
[nInputv] = preNorm(inputv,inputMin,inputMax);    
A=zeros(m,n-2);B=zeros(m,n-2);
Xt{1}=1;
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
figure_c=PLOTGP(T,Ys,Yprime,s2v,names,figure_c);


%% Simulatio Naive DATA
% load outputYS                                       % load test data                                         
% load inputUS                                       % load test data
% 
% [m,n]=size(Output);
% Output=Output';
% inputYS=[zeros(1,m);Output];
% inputYS(end,:)=[];
% inputYTS=inputYS(1:round(n/2),:);
% inputYVS=inputYS(round(n/2)+1:end,:);
% OutputTS=Output(1:round(n/2),:);
% OutputVS=Output(round(n/2)+1:end,:);
% inputU=inputU';
% inputUTS=inputU(1:round(n/2),:);
% inputUVS=inputU(round(n/2)+1:end,:);
% % inputUT=inputU;
% % OutputT=Output;
% % inputYT=inputY;
% [m,n]=size(inputUTS);
% 
% input = [inputYTS(:,2:end) inputUTS(:,2:3)] ; 
% target = OutputTS(:,2:end); 
% 
% 
% % normalize training data
% [nInput] = preNorm(input,inputMin,inputMax);
% [nTarget] = preNorm(target,targetMin, targetMax);



%% 

y_start=[nInput;nInputv];
%[nInput] = preNorm(y_start,inputMin,inputMax);  
%inputU(11,:)=[];
[m,n]=size(y_start);

AUX=[y_start(1,1:end-2)];
A=0;
B=0;
for i=1:m
%     [inputU(i,2:3)] = preNorm(inputU(i,2:3),inputMin,inputMax); 
%     [y_start] = preNorm(y_start,inputMin,inputMax);
    
%     [AUX] = preNorm(AUX,inputMin,inputMax);
    X_prime{1}=1;
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
names={'Udot m2/s Simulation';'V m/s Simulation';'r rad/s Simulation';'p rad/s Simulation',;'p rad/s Training'};
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
%% Neural Network
% load inputU;
% load outputY;
Ts = 3 ;% Sample time is 0.2 sec
z = iddata(Output(:,2:end),inputU(:,2:3),Ts);
ze = iddata(OutputT(:,2:end),inputUT(:,2:3),Ts);

zv = z;
na=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 2];
nb=[1 1;1 1;1 1;1 1];
nk=[0 0;0 0;0 0;0 0];
 
net = feedforwardnet(10);
NL=neuralnet(net);
Orders = [na,nb,nk];
sys = nlarx(ze,Orders,NL);
figure(21)
compare(zv, sys)
axh = findobj( gcf, 'Type', 'Line' );
xdata = get(axh, 'XData');  %data from low-level grahics objects
ydata = get(axh, 'YData');
xdata=xdata(9:end);
ydata=ydata(9:end);
%ydata=cell2mat(ydata);
ynn1=[cell2mat(ydata(7))',cell2mat(ydata(5))',cell2mat(ydata(3))',cell2mat(ydata(1))'];
%%
na=[4 0 0 0;0 4 0 0;0 0 4 0;0 0 0 4];
nb=[2 2;2 2;2 2;2 2];
nk=[0 0;0 0;0 0;0 0];
net = feedforwardnet(10);
NL=neuralnet(net);
Orders = [na,nb,nk];
sys = nlarx(ze,Orders,NL);
figure(21)
compare(zv, sys);
axh = findobj( gcf, 'Type', 'Line' );
xdata = get(axh, 'XData');  %data from low-level grahics objects
ydata = get(axh, 'YData');
xdata=xdata(9:end);
ydata=ydata(9:end);
%ydata=cell2mat(ydata);
ynn2=[cell2mat(ydata(7))',cell2mat(ydata(5))',cell2mat(ydata(3))',cell2mat(ydata(1))'];
%%
figure_c=PLOTGPNN(T,Ys,Yprime,ynn1,ynn2,s2_prime,names,figure_c);

%%





%% Montecarlo Simulation
% 
% [mu,sigma,MU]=mcsimulation(model,nInput(:,end-1:end),nInput(:,1:end-2));
% names={'Udot m/s Simulation MC';'V m/s Simulation MC';'r rad/s Simulation MC';'p rad/s Simulation MC'};
% [m,n]=size(nTarget);
% T=linspace(0,500,m)';
% PLOTGP(T,nTarget,mu,sigma,names,figure_c);

