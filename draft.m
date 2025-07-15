clc;clear
% 用来测试cluster 的CV哪些项是零的一段代码
tr = 3;
I = eye(2);
Z = [1,0;0,-1];

cosh(tr)
sinh(tr)

% % CM for TMSV states
V = [kron(I,cosh(tr)), kron(Z,sinh(tr));
     kron(Z,sinh(tr)), kron(I,cosh(tr))];

% V = Sigma;

Zk = kron(I, Z);
J = kron(I, [1,1i;1,-1i])/2;

Vt = Zk*J*V*J'*Zk;
X = kron(I, [0,1;1,0]);
Hin = -0.5*X*Vt;