function [h,q] = Double_layer_model(K1, K2, Fs,h1, h2)
% Domain information.
nx = 101;
nx1 = 70;
delta = 100;
h(1) = h1;
h(nx) = h2;
B(1) = -2*delta^2*Fs-K1*h1^2;
B(nx-2) = -2*delta^2*Fs-K2*h2^2;
B(2:nx-3) =  -2*delta^2*Fs;
x1(1:nx1-2) = K1;
x1(nx1-1:nx-3) = K2;
x2(1:nx1-2) = -2*K1;
x2(nx1:nx-2) = -2*K2;
x2(nx1-1) = -(K1 + K2);
x3(1:nx1-2) = K1;
x3(nx1-1:nx-3) = K2;
A = diag(x1,1) + diag(x2) + diag(x3, -1);
ht = B/A;
htt = ht.^(1/2);
h(2:nx-1) = htt;
for i = 1:nx-2
    if i < (nx-2)/2
    q(i) = -K1*(h(i+2)-h(i))/(2*delta);
    else
    q(i) = -K2*(h(i+2)-h(i))/(2*delta);
    end
end