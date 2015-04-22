function [TotalCost] = func(X,state)

alpha = 10; % Half Perimeter Wirelength MAGIC parameter % alpha = gridlength*radius
gridlength = 5; % Second MAGIC parameter % alpha = gridlength*radius
r = 2; % Another MAGIC parameter % alpha = gridlength*radius

ngates = state.ngates;
chipx = state.chipx;
chipy = state.chipy;
nnets = state.nnets;
pinc = state.pinc;
pinp = state.pinp;
nconn = state.nconn;
fundunit = state.fundunit;
b = state.b;

cellx = X(1:ngates);
celly = X((ngates+1):(2*ngates));


gridx(:,1) = 0:gridlength:chipx;
gridy(:,1) = 0:gridlength:chipy;
ngrids = (gridlength+1)^2; % No. of grids

temp1 = zeros(nnets,1);
temp2 = zeros(nnets,1);
temp3 = zeros(nnets,1);
temp4 = zeros(nnets,1);
halfperi = zeros(nnets,1);
for i = 1:nnets
    n = length(find(b(:,i)==1)) + length(find(pinc(:,i)==1));
    xi(1:n,1) = [cellx(b(:,i)==1);pinp((pinc(:,i)==1),1)]; % X coord of cells&pads connected to ith net
    yi(1:n,1) = [celly(b(:,i)==1);pinp((pinc(:,i)==1),2)];
    for j = 1:n
        temp1(i) = temp1(i) + exp(xi(j)/alpha);
        temp2(i) = temp2(i) + exp(-xi(j)/alpha);
        temp3(i) = temp3(i) + exp(yi(j)/alpha);
        temp4(i) = temp4(i) + exp(-yi(j)/alpha);
    end    
    halfperi(i,1) = alpha*(log(temp1(i)) + log(temp2(i)) + log(temp3(i)) + log(temp4(i)));
end

cost1 = sum(halfperi); % First part of cost function - Half-perimeter

area = zeros(ngates,1);
for i = 1:ngates
       area(i) = nconn(i)*fundunit;
end
capacity = sum(area)/ngrids; % TODO Verify
K = zeros(ngates,length(gridx),length(gridy));

[gx,gy,i] = ndgrid(gridx,gridy,1:ngates);


K = (area(i(:))/(r^2).*potential2(cellx(i(:)),celly(i(:)),gx(:),gy(:),r) - capacity).^2;


cost2 = sum(K); % Second part of cost function - Overlap Minimization
        
left = 0;
right = chipx;
top = chipy;
bottom = 0;
PenaltyLeft = zeros(ngates,1);
PenaltyRight = zeros(ngates,1);
PenaltyTop = zeros(ngates,1);
PenaltyBottom = zeros(ngates,1);
for i = 1:ngates
    if (cellx(i) < left)
        PenaltyLeft(i) = ((cellx(i) - left)/alpha)^2;
    elseif (cellx(i) > right)
        PenaltyRight(i) = ((cellx(i) - right)/alpha)^2;
    elseif (celly(i) < bottom)
        PenaltyBottom(i) = ((celly(i) - bottom)/alpha)^2;
    elseif (celly(i) > top)
        PenaltyTop(i) = ((celly(i) - top)/alpha)^2;
    end
end

cost3 = sum(sum(PenaltyLeft) + sum(PenaltyRight) + sum(PenaltyTop) + sum(PenaltyBottom));
TotalCost = cost1+cost2+cost3;
% Add weights

end
