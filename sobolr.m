function X0 = sobolr(N, P)
S = 1000;
L= 100;
p = sobolset(P, 'Skip', S, 'Leap', L);
p = scramble(p, 'MatousekAffineOwen');
X0 = net(p, N);
end