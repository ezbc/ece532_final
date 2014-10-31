close;
clear all;
clc;

N = 1000;

test = readtable('data/Test_rev1.xlsx');
train = readtable('data/Train_rev1.xlsx');
valid = readtable('data/Valid_rev1.xlsx');

sub_test = test(1:N, :);
sub_train = train(1:N, :);
sub_valid = valid(1:N, :);

save('data/test.mat', 'sub_test');
save('data/train.mat', 'sub_train');
save('data/valid.mat', 'sub_valid');

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






