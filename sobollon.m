function LY=sobollon(X0)
Y1 = icdf('Normal', X0, 0, 1);%%将均匀分布转换为正态分布
LY=exp(Y1);%%正态分布转换为对数正态分布
end