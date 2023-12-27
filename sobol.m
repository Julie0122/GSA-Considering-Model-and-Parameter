function LY=sobollon(X0)
Y1 = icdf('Normal', X0, 0, 1);%%将均匀分布转换为正态分布
LY=exp(Y1);%%正态分布转换为对数正态分布
end

function Y = soboln(X0)
Y = icdf('Normal', X0, 0, 1);%%将均匀分布转换为正态分布
end

function X0 = sobolr(N, P)%%参数个数
S = 1000;
L= 100;
p = sobolset(P, 'Skip', S, 'Leap', L);
p = scramble(p, 'MatousekAffineOwen');
X0 = net(p, N);
end