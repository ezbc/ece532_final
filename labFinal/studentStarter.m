%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141204
clc; clear all; close all; format compact;

%% Load Data
load('labData.mat');

%% Activity 1
% ------------------------------------------------------------------------------
fprintf('----Activity 1----\n');
desc1.id = 1; %from the activity1 training set
desc2.id = 2; %from the activity1 training set
desc3.id = 1; %from the test set

%count unique words
fprintf('desc1:\n');
desc1.words = strsplit(data.activity1.train.desc{desc1.id});
desc1.words = strrep(desc1.words, ',','');
desc1.words = strrep(desc1.words, '.','');
desc1.words = lower(desc1.words);
desc1.words = sort(desc1.words); 
desc1.uwords = unique(desc1.words);
for i = 1:length(desc1.uwords);
    desc1.a(i) = sum( strcmp(desc1.words, desc1.uwords{i}) );
    fprintf('%d = %s\n', desc1.a(i), desc1.uwords{i});
end

fprintf('\ndesc2:\n');
desc2.words = strsplit(data.activity1.train.desc{desc2.id});
desc2.words = strrep(desc2.words, ',','');
desc2.words = strrep(desc2.words, '.','');
desc2.words = lower(desc2.words);
desc2.words = sort(desc2.words); 
desc2.uwords = unique(desc2.words);
for i = 1:length(desc2.uwords);
    desc2.a(i) = sum( strcmp(desc2.words, desc2.uwords{i}) );
    fprintf('%d = %s\n', desc2.a(i), desc2.uwords{i});
end

fprintf('\ndesc3:\n');
desc3.words = strsplit(data.activity1.test.desc{desc3.id});
desc3.words = strrep(desc3.words, ',','');
desc3.words = strrep(desc3.words, '.','');
desc3.words = lower(desc3.words);
desc3.words = sort(desc3.words); 
desc3.uwords = unique(desc3.words);
for i = 1:length(desc3.uwords);
    desc3.a(i) = sum( strcmp(desc3.words, desc3.uwords{i}) );
    fprintf('%d = %s\n', desc3.a(i), desc3.uwords{i});
end

%combine word lists
comb = [desc1.uwords, desc2.uwords];
comb = unique(comb);
comb = sort(comb);

i1 = 1; %per-description counters
i2 = 1;
A = zeros(2,length(comb));
for i = 1:length(comb);
    if i1 <= length(desc1.uwords) && strcmp( comb{i}, desc1.uwords{i1} )
        A(1,i) = desc1.a(i1);
        i1 = i1+1;
    end
    if i2 <= length(desc2.uwords) && strcmp( comb{i}, desc2.uwords{i2} )
        A(2,i) = desc2.a(i2);
        i2 = i2+1;
    end
end
[m,n] = size(A);

fprintf('\nA = \n'); nA = 10;
for i = 1:nA-1; fprintf('%d & ',A(1,i)); end; fprintf('%d\\cdots\\\\\n', A(1,nA));
for i = 1:nA-1; fprintf('%d & ',A(2,i)); end; fprintf('%d\\cdots\n', A(2,nA));

% Display the frequency table
fprintf('\nFreq table 1\n')
for i = 1:11;
    disp([desc1.uwords{i} ' & ' num2str(A(1,i)) ' \\'])
end
fprintf('\nFreq table 2\n')
for i = 1:11;
    disp([desc2.uwords{i} ' & ' num2str(A(2,i)) ' \\'])
end



%% Activity 2
% ------------------------------------------------------------------------------
% Now predict salary of a new description
%...
fprintf('\n\n----Activity 2----\n');

%% Activity 3
% ------------------------------------------------------------------------------

%% Activity 4
% ------------------------------------------------------------------------------
fprintf('\n\n----Activity 4----\n');

%% Activity 5
% ------------------------------------------------------------------------------
% Number of data points
N = 100;

% Grab every word in the description
descTrain = data.activity4.train.desc(1:N);
descTest = data.activity4.test.desc(1:N);

%get the salary
salaryTrain = data.activity4.train.salary(1:N);
salaryTest = data.activity4.test.salary(1:N);

tic
all_words = true;
if all_words
    words = {'RemoveTheInitialWordFromwords'};
    for i = 1:N
        text = strsplit(descTrain{i}, ' ');
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
%remove the ignore words that probably dont contribute anything to the data
%set.
ignore = {'be' 'at' 'you' 'we' 'the' 'and' 'it' 'them' 'a' 'these' ...
          'those' 'with' 'can' 'for' 'an' 'is' 'or' 'of' 'are' 'has' 'have' ...
          'in' 'or' 'to' 'they' 'he' 'she' 'him' 'her' 'also'...
          '', 'able','all','as','but','by','cv','every','from','get','had','if','its',...
          'not','on','only','our','put','per','so','that','this','what','will','year','years','your'};
keywords = setdiff(keywords, ignore);
keywords = sort(keywords);

%calculate the frequency matrix
nKeys = length(keywords);
freqMatrixTrain = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTrain,keywords{ikeys});
    for idesc = 1:N;
        freqMatrixTrain(idesc, ikeys) = length(a{idesc}) / length(keywords{ikeys}) / length(descTrain{idesc});
    end
end

% Get frequencies of keywords, or A matrix, of the train Set
tic % time computation from here to toc
nKeys = length(keywords);
freqMatrixTrain = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTrain,keywords{ikeys});
    for idesc = 1:N;
        freqMatrixTrain(idesc, ikeys) = length(a{idesc}) / length(keywords{ikeys}) / length(descTrain{idesc});
    end
end

% Get frequencies of keywords, or A matrix, of the test Set
freqMatrixTest = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTest,keywords{ikeys});
    for idesc = 1:N;
        freqMatrixTest(idesc, ikeys) = length(a{idesc}) / length(keywords{ikeys}) / length(descTest{idesc});
    end
end
toc

%% Activity 5a) - Lasso Implementation
% ------------------------------------------------------------------------------
lambda = .1; 
maxIter = 1e3;
eps = 10^-5;

% Initilize the two weight vectors xhat and nextXhat to all zeros

% Take the SVD of the frequency matrix A

% Initilize alpha to 1/(largest eigen value of A)

% Set delta, usually the norm of xhat and nextXhat, to a large initial value

% Complete while loop
% while( ) % add conditions of loop
   
    %increase iteration number
    
    % set xhat to nextXhat;
    
    % compute y
    
    % compute nextXhat
   
    % compute delta as the 2-norm of xhat and nextXhat

    % delta = norm(xnext - xhat);
    
% end

%% Activity 5b) - using lasso with different lambda values and finding 
%                 words and erros
% ------------------------------------------------------------------------------
fprintf('\n\n----Activity 5b----\n');


% lambda = .001; 
% lambda = .00001; maxIter = 1e5;
% lambda = .1; maxIter = 1e4;
%N=300:
% lambda = .01; maxIter = 1e4;
% lambda = .001; maxIter = 1e5;
%N=500:
% lambda = .1; maxIter = 1e3;

%% Activity 5c) - Experimentation of choice of lambda
% ------------------------------------------------------------------------------
fprintf('\n\n----Activity 5c----\n');

%% Activity 6a) - Predict salaries of a Moder Major General and an EE Grad
% ------------------------------------------------------------------------------
% Major general
fprintf('\n\n----Activity 6a----\n');
Amjr = zeros(1,nKeys);
for i = 1:nKeys;
    Amjr(i) = length( strfind(data.activity6.major, keywords{i}) ) / length(keywords{i}) / length(data.activity6.major);
end

% EE grad
Aeeg = zeros(1,nKeys);
for i = 1:nKeys;
    Aeeg(i) = length( strfind(data.activity6.eengr, keywords{i}) ) / length(keywords{i}) / length(data.activity6.eengr);
end

%% Activity 6b)
% ------------------------------------------------------------------------------
fprintf('\n\n----Activity 6a----\n');


