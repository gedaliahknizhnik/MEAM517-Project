function [w_star,u_star] = getDes(t,t_stars,w_stars,u_stars,polys)
% GETDES(t,t_stars,w_stars,u_stars,polys) finds the the optimal 
%   trajectory point corresponding to the current time t.
% ---------------------------------------------------------------
% INPUTS:
%   t - the current time
%   t_stars - a vector of times corresponding to the SNOPT outputs w_stars
%       and u_stars
%   w_stars - a matrix of SNOPT state outputs
%   u_stars - a vector of SNOPT control outputs
%   polys - the polynomial trajectory output by SNOPT (currently unused)
% OUTPUTS:
%   w_star - the current optimal trajectory point.
%   u_star - the current optimal control point.

    % Interpolate between states using the cubic splines used by the
    % optimal trajectory generator.
    for jj=2:size(t_stars)
        if t <= t_stars(jj) 
            ind = jj-1;
            break
        end
    end

    dt = t - t_stars(ind);
    dts = [dt^3;dt^2;dt^1;1];
    w_star = squeeze(polys(ind,:,:))'*dts;

    % Linearly interpolate in control. 
    u_star = interp1(t_stars,u_stars,t);
end