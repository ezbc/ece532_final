%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141204
clc; clear all; close all; format compact;

%% Load Data
data = load('../data/warmup_train.mat');
train = data.warmup_train;
data = load('../data/warmup_test.mat');
test = data.warmup_test;
data = load('../data/warmup_valid.mat');
valid = data.warmup_valid;
clear data;

%% Activity 1
desc1.id = 1; %random pedagogical choice 1
desc2.id = 2;
desc3.id = 1;

%count unique words
fprintf('desc1:\n');
desc1.words = strsplit(train.FullDescription{desc1.id});
desc1.words = strrep(desc1.words, ',','');
desc1.words = strrep(desc1.words, '.','');
desc1.words = lower(desc1.words);
desc1.words = sort(desc1.words); 
desc1.uwords = unique(desc1.words);
for i = 1:length(desc1.uwords);
    desc1.a(i) = sum( strcmp(desc1.words, desc1.uwords{i}) );
    fprintf('%d = %s\n', desc1.a(i), desc1.uwords{i});
end

fprintf('desc2:\n');
%...

fprintf('desc3:\n');
%...

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

% Display the frequency table
fprintf('\n\nFreq table 1\n')
for i = 1:11;
    %...
end
fprintf('\n\nFreq table 2\n')
for i = 1:11;
    %...
end

%% Activity 2
% Now predict salary of a new description
%...

disp(['Norm of b - bhat = ' num2str(norm(b3-bhat))])

%% Activity 3

%% Activity 4

%% Activity 5

%loading data needed for this activity
data = load('../activity_data/activity4.mat');
train = data.train;
test = data.test;

% Grab every word in the description
descTrain = train.FullDescription(1:N);
descTest = test.FullDescription(1:N);

%get the salary
salaryTrain = train.SalaryNormalized(1:N);
salaryTest = test.SalaryNormalized(1:N);

% Calculating the frequency matrix (A) of the data
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
%remove the ignore words that probably dont contribute anything to the dataset.
ignore = {'be' 'at' 'you' 'we' 'the' 'and' 'it' 'them' 'a' 'these' ...
          'those' 'with' 'can' 'for' 'an' 'is' 'or' 'of' 'are' 'has' 'have' ...
          'in' 'or' 'to' 'they' 'he' 'she' 'him' 'her' 'also'...
          '', 'able','all','as','but','by','cv','every','from','get','had','if','its',...
          'not','on','only','our','put','per','so','that','this','what','will','year','years','your'};
keywords = setdiff(keywords, ignore);
keywords = sort(keywords);

%calculate the frequency matrix
nKeys = length(keywords);
freq_matrixTrain = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTrain,keywords{ikeys});
    for idesc = 1:N;
        freq_matrixTrain(idesc, ikeys) = length(a{idesc}) / length(keywords{ikeys}) / length(descTrain{idesc});
    end
end

%% Activity 5a) - Lasso Implementation
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
    
    %compute nextXhat
   
    % compute delta as the 2-norm of xhat and nextXhat
    delta = norm(xnext - xhat);
    
% end

%% Activity 5b) - using lasso with different lambda values and finding 
%                 words and erros


% lambda = .001; 
% lambda = .00001; maxIter = 1e5;
% lambda = .1; maxIter = 1e4;
%N=300:
% lambda = .01; maxIter = 1e4;
% lambda = .001; maxIter = 1e5;
%N=500:
% lambda = .1; maxIter = 1e3;

%% Activity 5c) - Experimentation of choice of lamda



%% Activity 6a)

% Predict salaries of a Moder Major General and an EE Grad
fid = fopen('majorGeneral.txt');
major = ' ';
for i = 1:37
    major = strcat(major, fgetl(fid));
end
fclose(fid);

% ...and an EE Grad


predMajor = Amjr * xhat
predEEGrad = Aeeg * xhat


%% Activity 6b)

