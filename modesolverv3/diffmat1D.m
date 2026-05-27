function D = diffmat1D(coord)

N = length(coord);
d = diff(coord); % differential estimation between coordinates

rows = 1:N-1;
cols1 = 1:N-1;
cols2 = 2:N;

vals1 = -1 ./ d;
vals2 =  1 ./ d;

D = sparse([rows rows], [cols1 cols2], [vals1 vals2], N, N); 

end
