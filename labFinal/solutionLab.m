%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141204
clc; clear all; close all; format compact;
set(0,'defaultaxescolor',.7*[1,1,1]);

%% Load Data
load('labData.mat');

%% Generate figure 1 - not required for students
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

%% Activity 1) - pinv() Solution
% {
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

fprintf('\nb matrix\n')
b = data.activity1.train.salary;
xhat = pinv(A)*b;
[xsrt,isrt] = sort(xhat,'descend');
for i = 1:11;
    if i < 11
        fprintf('%3.4f & ``%s" \\text{from} [%d, %d] \\\\ \n', xsrt(i), comb{isrt(i)}, A(1,isrt(i)), A(2,isrt(i)));
    end
end

%% Activity 2) - Estimate Salary of Third Description
fprintf('\n\n----Activity 2----\n');
% Now predict salary of a new description
i1 = 1; %per-description counters
A3 = zeros(1,length(comb));
for i = 1:length(comb);
    if i1 <= length(desc3.uwords) && strcmp( comb{i}, desc3.uwords{i1} )
        A3(1,i) = desc3.a(i1);
        i1 = i1+1;
    end
end

b3 = data.activity1.test.salary(desc3.id);

bhat = A3*xhat;
fprintf('Salary: Actual %3.2f vs Predicted %3.2f  difference = %3.2f\n', b3, bhat, norm(b3-bhat));
%}

%% Activity 3
% Following the same method as in the warm-up, describe how the
% frequency matrix will change with added descriptions. If we include
% every word from each of the descriptions in our frequency matrix,
% is it likely that the number of salaries will be greater than the
% number of words in the frequency matrix? Describe how the
% computation time scales with the size of the frequency matrix. Is
% it practical to include every word from each description in the
% frequency matrix?

%% Activity 4
% With the data in data.activity4, construct a frequency matrix of
% words, excluding common words you are confident will not be
% important to a description's salary. Derive a linear least squares
% solution for $\hat{\bm{x}}$ in Equation~\ref{eq:linsolve} with the
% first 10 descriptions in the training set. Predict the salaries of
% the first 10 descriptions in the test set. Report the $L_2$ norm of
% the difference between the true salaries and predicted salaries for
% the test set. Report the word which is best able to predict the
% salary of a description. State whether this problem is over
% determined or under determined.

% {
fprintf('\n\n----Activity 4----\n');
% Number of data points
N = 10;

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
tic
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

% Predict
xhat = pinv(freqMatrixTrain)*salaryTrain;
bhat = freqMatrixTest*xhat;


% Most useful words:
[xsrt,isrt] = sort(xhat,'descend');
for i = 1:11;
    if i < 11
        fprintf('%3.4f & ``%s"\n', xsrt(i), keywords{isrt(i)});
    end
end

fprintf('\nSalary: Actual vs. Predicted\n');
fprintf(' %3.2f vs %3.2f\n',salaryTest, bhat);
fprintf('L2 norm = %3.4f\n', norm(salaryTest-bhat));
%}

%% Activity 5 - Construct A
% Number of data points
N = 500;

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
tic
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
fprintf('\n\n----Activity 5a----\n');
A = freqMatrixTrain; b = salaryTrain;
xhat = zeros(size(A,2),1);
[u,s,v] = svd(A'*A);
alpha = 1/s(1,1);

lambda = 100; 
maxIter = 1e3;
eps = 10^-5;
delta = 10;
iter = 0;

tic
while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat);
    xhat = xnext;
end
toc
if iter >= maxIter;
    warning('MATLAB:LassoWarmup','Lasso exited at max iterations (%d) with delta %3.4f > eps %3.4f', maxIter, delta, eps);
end
xhat100 = xhat;

% Number of words
fprintf('For lambda = %d, there are %d words\n', lambda, numel(find(xhat ~= 0)));

% Display top key words
[~,ind] = sort(xhat,'descend');
for i = 1:(numel(find(xhat ~= 0))+1);
    fprintf('%d (%d) @ %3.4f = [%s]\n', i, ind(i), xhat(ind(i)), keywords{ind(i)} );
end

% Histogram
bhat = freqMatrixTest * xhat;
figure('Name','5a - Salary Residual Histogram'); hold on;
hist( abs(salaryTest - bhat), N );


%% Activity 5b) - using lasso with different lambda values
fprintf('\n\n----Activity 5b----\n');
%N=100
% lambda = .001; maxIter = 1e4; %max its, delta 68, error 1e5
% lambda = .00001; maxIter = 1e5; %max its, delat .7 error 1e5, 4.4min
% lambda = .1; maxIter = 1e4;
%N=300:
% lambda = .01; maxIter = 1e4; % max its, delta 1513, error 2e5, 5.2 min
% lambda = .001; maxIter = 1e5; %max its, delat 20, error 2e5, 32 min
%N=500:
% lambda = .1; maxIter = 1e3; %max its, delta 32940, error 3e5, 5 min

tic
lambda = 500;
xhat = zeros(size(A,2),1); iter = 0; delta = 10;
while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat);
    xhat = xnext;
end
xhat500 = xhat;

lambda = 10;
xhat = zeros(size(A,2),1); iter = 0; delta = 10;
while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat);
    xhat = xnext;
end
xhat010 = xhat;

lambda = 1;
xhat = zeros(size(A,2),1); iter = 0; delta = 10;
while( iter < maxIter && delta > eps)
    y = xhat + alpha*A'*(b-A*xhat);
    xnext = sign(y) .* max([abs(y) - alpha*lambda,zeros(size(A,2),1)],[],2);
    iter = iter+ 1;
    delta = norm(xnext - xhat);
    xhat = xnext;
end
xhat001 = xhat;
toc

r500 = salaryTest - freqMatrixTest*xhat500;
r100 = salaryTest - freqMatrixTest*xhat100;
r010 = salaryTest - freqMatrixTest*xhat010;
r001 = salaryTest - freqMatrixTest*xhat001;
[~,ind] = sort(r001);

figure('Name','5b - Salary Residuals'); hold on;
plot( 1:N, r500(ind), 'ko-');
plot( 1:N, r100(ind), 'bo-');
plot( 1:N, r010(ind), 'ro-');
plot( 1:N, r001(ind), 'co-');
xlabel('Sorted by residual'); ylabel('Salary Test - Predict');
legend('500','100','10','1');

%% Activity 5c) - 
% Experiment with the Lasso method using smaller values of $\lambda$.
% What $\lambda$ value do you think keeps the number of words
% reasonable while still including enough words in the weight matrix
% to accurately predict the salary? This can be an art, there is not

%% Activity 6a) - Predict salaries of a Moder Major General and an EE Grad
fprintf('\n\n----Activity 6a----\n');
Amjr = zeros(1,nKeys);
for i = 1:nKeys;
    Amjr(i) = length( strfind(data.activity6.major, keywords{i}) ) / length(keywords{i}) / length(data.activity6.major);
end

Aeeg = zeros(1,nKeys);
for i = 1:nKeys;
    Aeeg(i) = length( strfind(data.activity6.eengr, keywords{i}) ) / length(keywords{i}) / length(data.activity6.eengr);
end

[~,ind] = min([norm(r500),norm(r100),norm(r010),norm(r001)]);
xhat = [xhat500,xhat100,xhat010,xhat001]; xhat = xhat(:,ind);
predMajor = Amjr * xhat;
predEEGrad = Aeeg * xhat;

fprintf('Predicted Major General %3.4f\n', predMajor);
fprintf('Predicted Electrical Grad %3.4f\n', predEEGrad);

%% Activity 6b) - 
% The salaries fall substantially below the poverty line in any country;
% why is this the case?  Compare the most important keywords against
% those contained in \texttt{majorGen} and \texttt{elecGrad}; it may help
% to consider how the predicted salary is computed.

figure('Name','6b - Occurrences of Keywords'); hold on;
col = gray(N);
[~,inds] = sort(sum(freqMatrixTrain,2),'ascend');
for i = N:-1:1;
    ii = inds(i);
    ind = find(freqMatrixTrain(ii,:) > 0);
    plot( ind, freqMatrixTrain(ii,ind), ':', 'Color',col(i,:) );
end
ind = find(Amjr > 0); h1 = plot( ind, Amjr(ind), 'r*' ); set(h1,'DisplayName','Major');
ind = find(Aeeg > 0); h2 = plot( ind, Aeeg(ind), 'b*' ); set(h2,'DisplayName','EEGrad');
xlabel('Keyword Index'); ylabel('Counts');
legend([h1,h2])




