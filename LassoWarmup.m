clc;
close;
clear all;
format compact;

%% Generate Freq Matrix
% {
% Number of data points
N = 2;
tic

% load warmup data
data = load('data/warmup_train.mat');
train = data.warmup_train;
data = load('data/warmup_test.mat');
test = data.warmup_test;

descTrain = train(1:2, 3);
descTest = test(1:1,3);

salaryTrain = table2array (train(:,11) );
salaryTest = table2array (test(:,11) );

% Grab every word as a keyword
all_words = true;
if all_words
    words = {'RemoveTheInitialWordFromwords'};
    for i = 1:2
        desc = descTrain.FullDescription{i};
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
ignore = {'be' 'at' 'you' 'we' 'the' 'and' 'it' 'them' 'a' 'these' ...
          'those' 'with' 'can' 'for' 'an' 'is' 'or' 'of' 'are' 'has' 'have' ...
          'in' 'or' 'to' 'they' 'he' 'she' 'him' 'her' 'also'...
          '', 'able','all','as','but','by','cv','every','from','get','had','if','its',...
          'not','on','only','our','put','per','so','that','this','what','will','year','years','your'};
keywords = setdiff(keywords, ignore);
keywords = sort(keywords);

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
freq_matrixTest = zeros(1,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTest.FullDescription,keywords{ikeys});
    for idesc = 1:1;
        freq_matrixTest(idesc, ikeys) = length(a{idesc}) * length(keywords{ikeys}) / length(descTest.FullDescription{idesc});
    end
end

%% Lasso It
tic
lambda = 10; % gives 2 words
% lambda = 5; % gives 5 words
% lambda = 2; % gives 7 words
% lambda = 1; % gives 10 words



maxIter = 10000;
eps = 10^-3;

A = freq_matrixTrain; b = salaryTrain;
xhat = zeros(size(A,2),1);
[u,s,v] = svd(A'*A);
alpha = 1/s(1,1);
delta = 10;
iter = 0;

while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat);
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

% %Predict the salary based only on the words 
predSalaryTest = freq_matrixTest * xhat

% get error of predicted salary
error = norm(predSalaryTest - salaryTest)


toc