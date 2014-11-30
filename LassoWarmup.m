clc;
close;

% % Number of data points
% % N = 100;
% 
% % load Activity data
% data = load('data/activity_train.mat');
% train = data.activity_train;
% data = load('data/activity_test.mat');
% test = data.activity_test;
% % 
descTrain = train(:, 3);
descTest = test(:,3);

salaryTrain = table2array(train(:, 11));
salaryTest =table2array(test(:,11));
% 
% % Grab every word as a keyword
% all_words = true;
% if all_words
%     words = {' '};
%     for i = 1:N
%         text = strsplit(descTrain{i,1}{1}, ' ');
%         for j = 1:length(text)
%             %Remove all non English letters characters and numbers
%             text(j) = regexprep(text(j),'[^a-zA-Z0-9]','');   
%             % turn to all letters to  lower case
%             text(j) = lower(text(j));
%             %Add text to words
%             words = [words, text(j)];
%         end
%     end
%     keywords = unique(words);
% end
% 
% 
% % 
% % Get frequencies of keywords, or A matrix, of the train Set
% freq_matrixTrain = zeros(N, length(keywords));
% for i = 1:N
%     for j = 1:length(keywords)
%         freq = length(strfind(descTrain{i,1}{1}, keywords{j}));
%         sentence_length = length(descTrain{i,1}{1});
%         freq_matrixTrain(i, j) = freq /sentence_length;
%         
%         
%     end
% end
% 
% % Get frequencies of keywords, or A matrix, of the test Set
% freq_matriTest = zeros(N, length(keywords));
% for i = 1:N
%     for j = 1:length(keywords)
%         freq = length(strfind(descTest{i,1}{1}, keywords{j}));
%         sentence_length = length(descTest{i,1}{1});
%         freq_matrixTest(i, j) = freq /sentence_length;
%     end
% end
% 
lambda = 10^-2;
maxIter = 10;
eps = 10^-6;
xhat = Lasso( freq_matrixTrain, lambda,salaryTrain,maxIter,eps );
% xhat = lasso(freq_matrixTrain,salaryTrain);
% predSalaryTest = freq_matrixTest * xhat;
% 
