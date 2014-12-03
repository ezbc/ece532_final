clc;
close;
clear all;
format compact;

%% Generate Freq Matrix
% {
% Number of data points
N = 100;
tic

% load Activity data
data = load('data/activity_train.mat');
train = data.activity_train;
data = load('data/activity_test.mat');
test = data.activity_test;

descTrain = train(1:N, 3);
descTest = test(1:N,3);

salaryTrain = train.SalaryNormalized(1:N);
salaryTest = test.SalaryNormalized(1:N);

% Grab every word as a keyword
all_words = true;
if all_words
    words = {'RemoveTheInitialWordFromwords'};
    for i = 1:N
        text = strsplit(descTrain.FullDescription{i}, ' ');
        for j = 1:length(text)
            %Remove all non English letters characters and numbers
            text{j} = regexprep(text{j},'[^a-zA-Z0-9]','');
            % turn to all letters to lower case
            text{j} = lower(text{j});
            %Add text to words
            if length(text{j}) > 2
                words = [words, text{j}]; %#ok<AGROW>
            end
        end
    end
    words = words(2:end); %remove the initial
    keywords = unique(words);
end
keywords = sort(keywords);
toc

tic
% Get frequencies of keywords, or A matrix, of the train Set
%{
freq_matrixTrain = zeros(N, length(keywords));
for i = 1:N
    for j = 1:length(keywords)
        freq = length(strfind(descTrain.FullDescription{i}, keywords{j})) * length(keywords{j});
        sentence_length = length(descTrain.FullDescription{i});
        freq_matrixTrain(i, j) = freq / sentence_length;
    end
end
toc
%}

%this is functionally equivalent but twice the speed...
tic
nKeys = length(keywords);
freq_matrixTrain = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTrain.FullDescription,keywords{ikeys});
    for idesc = 1:N;
        freq_matrixTrain(idesc, ikeys) = length(a{idesc}) * length(keywords{ikeys}) / length(descTrain.FullDescription{idesc});
    end
end
toc
    
    
% Get frequencies of keywords, or A matrix, of the test Set
freq_matrixTest = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTest.FullDescription,keywords{ikeys});
    for idesc = 1:N;
        freq_matrixTest(idesc, ikeys) = length(a{idesc}) * length(keywords{ikeys}) / length(descTest.FullDescription{idesc});
    end
end

%% Lasso It
tic
lambda = 10; 
maxIter = 100;
eps = 10^-6;
% xhat = Lasso( freq_matrixTrain, lambda,salaryTrain,maxIter,eps );
% xhat = Lasso2( freq_matrixTrain, lambda,salaryTrain,maxIter,eps ); 
A = freq_matrixTrain; b = salaryTrain;
xhat = zeros(size(A,2),1);
% size(lambda*(A'*A) )
% alpha = 1./(lambda*(A'*A)); %Something seems wrong here.
[u,s,v] = svd(A'*A);
alpha = 1/s(1,1);
delta = 10;
iter = 0;

while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat)
    xhat = xnext;
end
if iter >= maxIter;
    warning('MATLAB:LassoWarmup','Lasso exited at max iterations (%d) with delta %3.4f > eps %3.4f', maxIter, delta, eps);
end;

%Find and display key words by checking xhat values greater than eps
[~,ind] = sort(xhat,'descend');
for i = 1:30;
    fprintf('%d @ %3.4f = [%s]\n', i, xhat(ind(i)), keywords{ind(i)} );
end

%Predict the salary based only on the words 
% predSalaryTest = freq_matrixTest * xhat;

%get error of predicted salary
% error = norm(predSalaryTest - salaryTest)

toc