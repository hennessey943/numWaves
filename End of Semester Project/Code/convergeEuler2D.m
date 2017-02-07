function convergeEuler2D(Case,iPlot)
if Case==1
    %smooth initial data
    
    %Richardson Extrapolation
elseif Case==2
    %Riemann initial data (shock tube)
    %Exact Error comparison
end

if iPlot==1
    %Plot Solution
    if Case==2
        %Plot error and solution
        surf(x,y,u(1,:,:))
        return
    end
    
end