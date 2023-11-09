FileData = load("shape.mat");
M = [x;y;z];
M = transpose(M);
csvwrite(['NX_test.csv'],M)