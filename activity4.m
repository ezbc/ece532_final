%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141105
clc; clear all; close all; format compact;

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
toc

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

% Perform least squares fit
xhat = pinv(freq_matrixTrain)*salaryTrain;

salaryTest_hat = freq_matrixTest * xhat;

disp('L2 norm of predicted vs. test data set')
disp(norm(salaryTest - salaryTest_hat)) 

disp('Best predictor of salary')
disp(keywords(xhat == max(xhat)))

% Should be 'surrounding' with a weight of 5e5




