function [polys] = getPolys(states,controls,nx,dt,params)

    polys = zeros(size(states,1)-1,4,nx);

    for ii=1:(size(states,1)-1)     
        % Extract points
        x_i   = states(ii,:);
        x_ip1 = states(ii+1,:);
        f_i   = f(states(ii,:),controls(ii,:),params)';
        f_ip1 = f(states(ii+1,:),controls(ii+1,:),params)';
        
        
        % Create the polynomials (in descending order of power)
        polys(ii,4,:) = x_i;
        polys(ii,3,:) = f_i;
        polys(ii,2,:) = 3/dt^2*(x_ip1-x_i) - 1/dt*(f_ip1 + 2*f_i);
        polys(ii,1,:) = -2/dt^3*(x_ip1-x_i) + 1/dt^2*(f_i + f_ip1);
    end
    
    
end