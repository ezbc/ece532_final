close;
clear all;
clc;

%train = load('data/train.m');
%test = load('data/test.m');
%valid = load('data/valid.m');

%data = containers.Map;
%sub_data = containers.Map;

% Train
%sub_data('desc') = train(:, 3);
%sub_data('salary_raw') = train(:, 10);
%sub_data('salary_norm') = train(:, 11);
%data('train') = sub_data;

% % Test
% sub_data('desc') = test(:, 3);
% sub_data('salary_raw') = test(:, 10);
% sub_data('salary_norm') = test(:, 11);
% data('test') = sub_data;
%
% % valid
% sub_data('desc') = valid(:, 3);
% sub_data('salary_raw') = valid(:, 10);
% sub_data('salary_norm') = valid(:, 11);
% data('valid') = sub_data;

N = 100;

% get subset of data
%full_train = train;
%train = full_train(1:N, :);

data = load('data/train.mat');
train = data.sub_train;
data = load('data/test.mat');
test = data.sub_test;
data = load('data/valid.mat');
valid = data.sub_valid;

desc = train(:, 3);


% Grab every word as a keyword
all_words = false;
if all_words
    words = {' '};
    for i = 1:N
        text = strsplit(desc{i,1}{1}, ' ');
        for j = 1:length(text)
            % Remove periods
            text(j) = strrep(text(j), '.', '');
            % Add text to words
            words = [words, text(j)];
        end
    end
    keywords = unique(words);
end

% List keywords to count
keywords = {'math' 'mathematical' 'systems' 'analyst' 'and'};

% Get frequencies of keywords, or A matrix
freq_matrix = zeros(N, length(keywords));
for i = 1:N
    for j = 1:length(keywords)
        freq = length(strfind(desc{i,1}{1}, keywords{j}));
        sentence_length = length(desc{i,1}{1});
        freq_matrix(i, j) = freq /sentence_length;
    end
end

% Get salaries, or b matrix
salary = table2array(train(:, 11));

[u, s, v] = svd(freq_matrix, 'econ');

x_hat = pinv(freq_matrix) * salary;

% Get freq matrix of valid data
freq_matrix_valid = zeros(N, length(keywords));
for i = 1:N
    for j = 1:length(keywords)
        freq = length(strfind(desc{i,1}{1}, keywords{j}));
        sentence_length = length(desc{i,1}{1});
        freq_matrix_valid(i, j) = freq /sentence_length;
    end
end

% 









