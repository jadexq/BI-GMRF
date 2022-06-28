clear; clc;

%% load tensor_toolbox "MATLAB Tensor Toolbox, Version 3.1"
addpath('./tensor_toolbox-v3.1');
savepath;

% load gray scale images
% load alpha true shape
I1 = imread('./input/alphaTrue.png');
disp(size(I1)); % figure; imagesc(I1);
% load beta true shape
I2 = imread('./input/betaTrue.png');
disp(size(I2)); % figure; imagesc(I2);

%%
% Downsize using interpolation
% Set target dimension
targetdim = round([64 64]);
disp(targetdim);
% alpha true
I1_down = array_resize(I1, targetdim); % default method is interpolate
disp(size(I1_down)); % figure; imagesc(I1_down);
% beta true
I2_down = array_resize(I2, targetdim);% default method is interpolate
disp(size(I2_down)); % figure; imagesc(I2_down);


%% Generate data for 2D simulation
% m_i (v) = alpha(v) * x_i + lambda1(v) * c_1i + epsilon_i(v)
% y_i = sum_v beta(v) m_i(v) + gamma * x_i + delta_i

% constants
nRep = 1; % number of replications
n = 500; % sample size
v1 = 64; % dimension of 2D images
v2 = 64;
vv = v1*v2;

% true values
xx1 = -0.2; % alpha true value scale
xx2 = 0.2; % beta true value scale
sigE = 1; % sd of epsilon_i(v)
sigD = 0.1; % sd of delta_i
gamma = -0.3; % gamma true value
% alpha and beta true images
alpha = 1 - (I1_down == 255);
alpha = xx1 * alpha;
figure; imagesc(alpha);
beta = 1 - (I2_down == 255);
beta = xx2* beta;
figure; imagesc(beta);
% lambda1 true value
lambda1 = normrnd(0.0, 0.03, [v1, v2]);
% vectorize by col
alphaVec = reshape(alpha, 1, []);
betaVec = reshape(beta, 1, []);
lambda1Vec = reshape(lambda1, 1, []);

%% generate data
x = zeros(n*nRep, 1);
c1 = zeros(n*nRep, 1);
y = zeros(n*nRep, 1);
m = zeros(n*nRep, vv);

for r = 1:nRep
    % x and z
    x_r = round(rand(1,n)); % binary p=0.5
    c1_r = normrnd(0,1,[1,n]); % std normal

    % m_i (v) = alpha(v) * x_i + lambda1(v) * c_1i + epsilon_i(v)
    m_r = zeros(n, vv);
    for i = 1:n
        E = normrnd(0, sigE, [1, vv]);
        m_r(i,:) = alphaVec * x_r(i) + lambda1Vec * c1_r(i) + E;
    end

    % y_i = sum_v beta(v) m_i(v) + gamma * x_i + delta_i
    y_r = zeros(n,1);
    for i = 1:n
        sumB = dot(betaVec, m_r(i,:));
        y_r(i) = sumB + gamma * x_r(i) + normrnd(0, sigD);
    end
    % assign
    x((n*(r-1)+1):(n*r), 1) = x_r';
    c1((n*(r-1)+1):(n*r), 1) = c1_r';
    y((n*(r-1)+1):(n*r), 1) = y_r';
    m((n*(r-1)+1):(n*r), :) = m_r;
end

%% write .txt
% x, c1, y, m
dlmwrite('./output/x_gen.txt', x);
dlmwrite('./output/c1_gen.txt', c1);
dlmwrite('./output/y_gen.txt', y);
dlmwrite('./output/m_gen.txt', m ,' ');
% alpha, lambda1, beta true values
dlmwrite('./output/alpha_true.txt', alpha ,' '); % as matrix
dlmwrite('./output/alphaVec_true.txt', alphaVec' ,' '); % as vector
dlmwrite('./output/beta_true.txt', beta ,' ');
dlmwrite('./output/betaVec_true.txt', betaVec' ,' ');
dlmwrite('./output/lambda1_true.txt', lambda1 ,' ');
dlmwrite('./output/lambda1Vec_true.txt', lambda1Vec' ,' ');
dlmwrite('./output/gamma_true.txt', gamma);




