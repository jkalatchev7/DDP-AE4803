function [Cost] = calcCost_ec(x_nom,u_nom,target,Q,Q_f,R,dt)
[numOfStates,horizon] = size(x_nom);
 Cost = 0;
 
 for j =1:(horizon-1)
     
    Cost = Cost + 0.5 * u_nom(:,j)' * R * u_nom(:,j) * dt + 0.5 * (x_nom(:,j)-target(:,j))' * Q * (x_nom(:,j)-target(:,j)) * dt;
     
 end
 
 TerminalCost= (x_nom(:,horizon) - target(:,horizon))'*Q_f * (x_nom(:,horizon) - target(:,horizon));
 
 Cost = Cost + TerminalCost;
end