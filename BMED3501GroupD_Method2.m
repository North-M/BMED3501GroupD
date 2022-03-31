close all;


%define the basic parameter
theta = (0:0.4:179.6);
N = 256;
%normalize the image to calculate cost function
normal = @(X) (X - min(min(X))) / (max(max(X)) - min(min(X)));
%define radon transform
A = @(x) radon(x,theta);
A_T = @(x) imresize(iradon(x,theta,"pchip","none"),[N,N]);
%define cost function
CostFun = @(x,y) sum(sum(power((x-y),2)));
%generate a phantom image and calculate the radon transform
x = phantom(N);
p = A(x);

%define the parameter for iteration
epoch = 1000;
lambda = 0.1;

L = 5e-4;
I = zeros(N,N);
t = 1;
Ip = I;
y = I;
costp = 99999;
%iteration
for i = 1:epoch
    d = y - L * A_T(A(y) - p);
    I = max((abs(d) - lambda * L),0) .* sign(d);
    tp = (1+sqrt(1+4*t^2))/2;
    y = I + (t - 1)/tp * (I - Ip);

    cost = CostFun(y,x);
    if mod(i,100) == 0
        figure(1);hold;
        subplot(2,5,i/100);
        imshow(I);
        title(sprintf("%d iteration cost=%.2f",i,CostFun(y,x)));
    end
    if cost > costp
        figure(2);imshow(Ip);title(sprintf("pchip: %d iteration cost=%.2f",i - 1,costp));
        break;
    end
    Ip = I;
    t = tp;
    costp = cost;
end
