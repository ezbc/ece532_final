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

%% WarmUp
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
desc2.words = strsplit(train.FullDescription{desc2.id});
desc2.words = strrep(desc2.words, ',','');
desc2.words = strrep(desc2.words, '.','');
desc2.words = lower(desc2.words);
desc2.words = sort(desc2.words); 
desc2.uwords = unique(desc2.words);
for i = 1:length(desc2.uwords);
    desc2.a(i) = sum( strcmp(desc2.words, desc2.uwords{i}) );
    fprintf('%d = %s\n', desc2.a(i), desc2.uwords{i});
end

fprintf('desc3:\n');
desc3.words = strsplit(test.FullDescription{desc3.id});
desc3.words = strrep(desc3.words, ',','');
desc3.words = strrep(desc3.words, '.','');
desc3.words = lower(desc3.words);
desc3.words = sort(desc3.words); 
desc3.uwords = unique(desc3.words);
for i = 1:length(desc3.uwords);
    desc3.a(i) = sum( strcmp(desc3.words, desc3.uwords{i}) );
    fprintf('%d = %s\n', desc3.a(i), desc3.uwords{i});
end

% Define common words that should be ignored in frequency matrix
ignore = {'be' 'at' 'you' 'we' 'the' 'and' 'it' 'them' 'a' 'these' ...
          'those' 'with' 'can' 'for' 'an' 'is' 'or' 'of' 'are' 'has' 'have' ...
          'in' 'or' 'to' 'they' 'he' 'she' 'him' 'her' 'also'};

% remove the ignored words
%ignore_indices = find(ismember(desc1.uwords, ignore));
%desc1.uwords(ignore_indices) = [];
%ignore_indices = find(ismember(desc2.uwords, ignore));
%desc2.uwords(ignore_indices) = [];

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

% for i = 1:n; fprintf('%d %d %s\n',A(1,i),A(2,i), comb{i}); end;
fprintf('A = \n'); nA = 10;
for i = 1:nA-1; fprintf('%d & ',A(1,i)); end; fprintf('%d\\cdots\\\\\n', A(1,nA));
for i = 1:nA-1; fprintf('%d & ',A(2,i)); end; fprintf('%d\\cdots\n', A(2,nA));
fprintf('\n\n');




% Display the frequency table
fprintf('\n\nFreq table 1\n')
for i = 1:11;
    disp([desc1.uwords{i} ' & ' num2str(A(1,i)) ' \\'])
end
fprintf('\n\nFreq table 2\n')
for i = 1:11;
    disp([desc2.uwords{i} ' & ' num2str(A(2,i)) ' \\'])
end

fprintf('\n\nb matrix\n')
b = [train.SalaryNormalized(desc1.id); train.SalaryNormalized(desc2.id)];
xhat = pinv(A)*b;
[xsrt,isrt] = sort(xhat,'descend');
for i = 1:11;
    if i < 11
        fprintf('%3.4f & ``%s" \\text{from} [%d, %d] \\\\ \n', xsrt(i), comb{isrt(i)}, A(1,isrt(i)), A(2,isrt(i)));
    end
end

% Now predict salary of a new description
i1 = 1; %per-description counters
A3 = zeros(1,length(comb));
for i = 1:length(comb);
    if i1 <= length(desc3.uwords) && strcmp( comb{i}, desc3.uwords{i1} )
        A3(1,i) = desc3.a(i1);
        i1 = i1+1;
    end
end

b3 = test.SalaryNormalized(desc3.id);

bhat = A3*xhat;

disp(['Norm of b - bhat = ' num2str(norm(b3-bhat))])


%% LSE figure
%{
x = 0:.01:1;
data = randn(1,length(x)) + .5;
fit = polyfit(x,data,1);
y = polyval(fit,x);

msz = 20; lwd = 3;
figure(); hold on;
plot(x, data, 'b.','MarkerSize',msz);
plot(x, y, 'r-','LineWidth',lwd);
xlim([0,1]); ylim([0,1]);
fs = 20;
ylabel('Normalized Salary','FontSize',fs); xlabel('Normalized Occurrences of `analyse`','FontSize',fs);
set(gca,'FontSize',fs);
%}

%% Regularized SVD

%% Lasso Regularization

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

%% Calculating the frequency matrix (A) of the data
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
freq_matrixTrain = zeros(N,nKeys);
for ikeys = 1:nKeys;
    a = strfind(descTrain,keywords{ikeys});
    for idesc = 1:N;
        freq_matrixTrain(idesc, ikeys) = length(a{idesc}) / length(keywords{ikeys}) / length(descTrain{idesc});
    end
end
%% Activity 5a) - Lasso Implementation

% Initilize the two weight vectors xhat and xhat' to all zeros

% Take the SVD of the frequency matrix A

% Initilize alpha to 1/(largest eigen value of A)

% Set delta, usually the norm of xhat and xhat', to a large initial value

% Complete while loop
while(delta == inf ) %% add conditions of loop
   
    %increase iteration number
    
    % set xhat to xhat';
    
    % compute y
    
    %compute xhat'
   
    % compute delta as the 2-norm of xhat and xhat'
    delta = norm(xnext - xhat);
    
end

%% Activity 5b) - using lasso with different lambda values and finding 
%                 words and erros


% lambda = .001; 
% lambda = .00001; maxIter = 1e5;
lambda = .1; maxIter = 1e4;
%N=300:
% lambda = .01; maxIter = 1e4;
% lambda = .001; maxIter = 1e5;
%N=500:
% lambda = .1; maxIter = 1e3;

%% Activity 5c) - Experimentation of choice of lamda


%% Ramifications
Pb = A*A';
Ps = A'*A;
