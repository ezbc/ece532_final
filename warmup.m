%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141105
clc; clear all; close all; format compact;

%% Load Data
data = load('data/warmup_train.mat');
train = data.warmup_train;
data = load('data/warmup_test.mat');
test = data.warmup_test;
data = load('data/warmup_valid.mat');
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

% % % Define common words that should be ignored in frequency matrix
% ignore = {'be' 'at' 'you' 'we' 'the' 'and' 'it' 'them' 'a' 'these' ...
%           'those' 'with' 'can' 'for' 'an' 'is' 'or' 'of' 'are' 'has' 'have' ...
%           'in' 'or' 'to' 'they' 'he' 'she' 'him' 'her' 'also'};
% 
% % remove the ignored words
% ignore_indices = find(ismember(desc1.uwords, ignore));
% desc1.uwords(ignore_indices) = [];
% ignore_indices = find(ismember(desc2.uwords, ignore));
% desc2.uwords(ignore_indices) = [];

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
% 
% % Display the frequency table
% fprintf('\n\nFreq table 1\n')
% for i = 1:11;
%     disp([desc1.uwords{i} ' & ' num2str(A(1,i)) ' \\'])
% end
% 
% fprintf('\n\nFreq table 2\n')
% for i = 1:11;
%     disp([desc2.uwords{i} ' & ' num2str(A(2,i)) ' \\'])
% end
% 
fprintf('\n\nb matrix\n')
b = [train.SalaryNormalized(desc1.id); train.SalaryNormalized(desc2.id)];
% xhat = pinv(A)*b;
xhat = (A'*A)\A'*b;

% [xsrt,isrt] = sort(xhat,'descend');
% for i = 1:length(comb);
%     
% %     foundIn = '';
% %     if A(1,i) ~= 0; foundIn = [foundIn,'d1 ']; end;
% %     if A(2,i) ~= 0; foundIn = [foundIn,'d2']; end;
% %     fprintf('%d: %3.4f = "%s" from [%s]\n', i, xsrt(i), comb{i}, foundIn);
% 
% %     fprintf('%d: %3.4f = "%s" from [%d %d]\n', i, xsrt(i), comb{i}, A(1,i), A(2,i));
%     if i < 11
%         fprintf('%3.4f & "%s" from [%d, %d] \\\\ \n', xsrt(i), comb{isrt(i)}, A(1,isrt(i)), A(2,isrt(i)));
%     end
% end
% 
% % ---------------------------------------
% % Now predict salary of a new description
% % ---------------------------------------

i1 = 1; %per-description counters
A3 = zeros(1,length(comb));
for i = 1:length(comb);
    if i1 <= length(desc3.uwords) && strcmp( comb{i}, desc3.uwords{i1} )
        A3(1,i) = desc3.a(i1);
        i1 = i1+1;
    end
end

bhat = A3*xhat

disp(['Norm of b - bhat = ' num2str(norm(test.SalaryNormalized(desc3.id)-bhat))])




%% LSE figure
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