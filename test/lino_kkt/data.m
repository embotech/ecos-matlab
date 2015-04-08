%Generate random data for KKT-System

%Maximal dimensions
maxEq = 100;
maxStates = 1000;
maxConeDim = 100;

%random dimensions
n = randi(maxEq);
m = randi(maxStates)+(n-1);
k = randi(maxConeDim);

%random A and G
A = rand(n,m)*200-100; %random A with values [-100,100]
G = rand(k,m)*200-100; %random G with values [-100,100]

%random s and z with s(0), z(0) [0,100]
s = zeros(k,1);
z = zeros(k,1);

s(1) = rand*100;
z(1) = rand*100;

for i=2:k
    s(i) = 2*rand*s(1)/(k-1)-(s(1)/(k-1));
    z(i) = 2*rand*z(1)/(k-1)-(z(1)/(k-1));
end
assert(s(1)-norm(s(2:k)) >= 0,'s not in SOCP');
assert(z(1)-norm(z(2:k)) >= 0,'z not in SOCP');

%RHS with random values [-100,100]
bx = rand(m,1)*200-100;
by = rand(n,1)*200-100;
bz = rand(k,1)*200-100;


