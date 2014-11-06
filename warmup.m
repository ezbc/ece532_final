%% Elijah Bernstein-Cooper, Ahmed Saif, Ben Conrad - ECE532 Project - 141105
clc; clear all; close all; format compact;

%% Load Data
data = load('data/train.mat');
train = data.sub_train;
data = load('data/test.mat');
test = data.sub_test;
data = load('data/valid.mat');
valid = data.sub_valid;
clear data;

%% WarmUp
desc1.id = 1; %random pedagogical choice 1
desc2.id = 10;

%count unique words
fprintf('desc1:\n');
desc1.words = strsplit(test.FullDescription{desc1.id});
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
desc2.words = strsplit(test.FullDescription{desc2.id});
desc2.words = strrep(desc2.words, ',','');
desc2.words = strrep(desc2.words, '.','');
desc2.words = lower(desc2.words);
desc2.words = sort(desc2.words); 
desc2.uwords = unique(desc2.words);
for i = 1:length(desc2.uwords);
    desc2.a(i) = sum( strcmp(desc2.words, desc2.uwords{i}) );
    fprintf('%d = %s\n', desc2.a(i), desc2.uwords{i});
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

b = [test.SalaryNormalized(desc1.id); test.SalaryNormalized(desc2.id)];
xhat = pinv(A)*b;
[xsrt,isrt] = sort(xhat,'descend');
for i = 1:length(comb);
    
%     foundIn = '';
%     if A(1,i) ~= 0; foundIn = [foundIn,'d1 ']; end;
%     if A(2,i) ~= 0; foundIn = [foundIn,'d2']; end;
%     fprintf('%d: %3.4f = "%s" from [%s]\n', i, xsrt(i), comb{i}, foundIn);

%     fprintf('%d: %3.4f = "%s" from [%d %d]\n', i, xsrt(i), comb{i}, A(1,i), A(2,i));
    fprintf('%d: %3.4f = "%s" from [%d %d]\n', i, xsrt(i), comb{isrt(i)}, A(1,isrt(i)), A(2,isrt(i)));
end







