infile = 'benchmarks/toy/toy1';
outfile = 'outfile.txt';

fd = fopen(infile);

chipx = fscanf(fd,'%f',1);
chipy = fscanf(fd,'%f',1);
fundunit = fscanf(fd,'%f',1);

ngates = fscanf(fd,'%f',1);
nnets = fscanf(fd,'%f',1);

b = zeros(ngates,nnets);

for i = 1:ngates
    gnum = fscanf(fd,'%f',1);
    assert(i==gnum);
    nconn = fscanf(fd,'%f',1);
    c_net = fscanf(fd,'%f', nconn); 
    b(gnum,c_net) = 1;
end

npins = fscanf(fd,'%f',1);

pinc = zeros(npins,nnets);
pinp = zeros(npins,2);

for i = 1:npins
    pnum = fscanf(fd,'%f',1);
    assert(i==pnum);
    net = fscanf(fd,'%f',1);
    pinc(pnum,net) = 1;
    pos = fscanf(fd,'%f',2);
    pinp(pnum,:) = pos;
end

fclose(fd);

% Assemble connectivity matrix

c = zeros(ngates,ngates);

% TODO: Verify
for i = 1:nnets
    c(b(:,i)==1,b(:,i)==1) += 1;
end

% Gates don't self-connect
c(logical(eye(size(c)))) = 0;

cdiag = sum(c,1)';
for i = 1:nnets
    cdiag(b(:,i)==1) += sum(pinc(:,i));
end


A = diag(cdiag) - c;

bxy = zeros(ngates,2);
for i = 1:nnets
    if nnz(pinc(:,i)) ~= 0
        % This will break if more than one net connects a gate and a pin
        bxy(b(:,i)==1,:) = repmat(-pinp(pinc(:,i) == 1,:),nnz(b(:,i)),1);
    end
end


x_soln = A \ -bxy(:,1);
y_soln = A \ -bxy(:,2);

fd = fopen(outfile,'w');
fprintf(fd,'%d %f %f\n',[ 1:size(x_soln,1) ; x_soln'; y_soln'])
fclose(fd);

