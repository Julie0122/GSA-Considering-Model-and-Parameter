function Y = soboln(X0)
Y = icdf('Normal', X0, 0, 1);%%将均匀分布转换为正态分布
end