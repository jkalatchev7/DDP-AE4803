function [Cost] = calcCost(x_nom,u_nom,target,Q,Q_f,R,dt)
[numOfStates,horizon] = size(x_nom);
 Cost = 0;
 
 for j =1:(horizon-1)
     
    Cost = Cost + 0.5 * u_nom(:,j)' * R * u_nom(:,j) * dt + 0.5 * (x_nom(:,j)-target)' * Q * (x_nom(:,j)-target) * dt;
     
 end
 
 TerminalCost= (x_nom(:,horizon) - target)'*Q_f * (x_nom(:,horizon) - target);
 
 Cost = Cost + TerminalCost;
end