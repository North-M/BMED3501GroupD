close all

%normalize the image to calculate cost function
normal = @(X) (X - min(min(X))) / (max(max(X)) - min(min(X)));

%calculate the cost function by summing the square difference
CostFun = @(x,y) sum(sum(power((x-y),2)));

%create a phaontom image in 256x256 matrix
I = phantom(256); 
NI = normal(I);


%conduct radon transform
theta = (0:0.4:179.6);
R = radon(I,theta);

%add white Gaussian noise to signal
R_noise = awgn(R,20,'measured');

%filter the signal by hamming window
win = zeros(3,3);
h = hamming(7);
for i = (1:3)
    for j = (1:3)
        win(i,j) = h(4 - round(sqrt((i-2) * (i-2) + (j-2) * (j-2))));
    end
end


FilterR = filter2(win,R_noise);

%differentiate R with x(theta) and y(s)
[DX, DY] = gradient(R);
DY = DY / (pi);
DYY = diff(R);
DYY = DYY / 2 / pi;

[DX_noise, DY_noise] = gradient(R_noise);

%conduct hilbert tranform with y(s), the image part of the output H is the
%result of hilbert transform
H = hilbert(DY);
H_noise = hilbert(DY_noise);
HI = imag(H);
HI_noise = imag(H_noise);


%backprojection
F = iradon(HI, theta, 'linear', 'none');
F_noise = iradon(HI_noise,theta,'linear','none');

%normalize the output image
NF = normal(imresize(F,[256,256]));
NF_noise = normal(imresize(F_noise,[256,256]));

%calculate the cost function for image with noise and filtered image
CostFun(NI,NF_noise)
CostFun(I,imresize(F,[256,256]))

figure(1);imagesc(I);colormap(gray(64));colorbar;title("original image");
figure(2);imagesc(R);colormap(gray(64));colorbar;title("Rf");
figure(3);imagesc(DY);colormap(gray(64));colorbar;title("D(g(s,θ))");
figure(4);imagesc(HI);colormap(gray(64));colorbar;title("H(D(Rf))");
figure(5);imagesc(NF);colormap(gray(64));colorbar;title("reconstructed image");
figure(6);imagesc(NF_noise);colormap(gray(64));colorbar;title("With noise");