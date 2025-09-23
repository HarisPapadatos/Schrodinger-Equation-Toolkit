classdef Calculus
    methods (Static)
        function L = laplacian1D(n, width, periodic)
            if periodic
                h = width / n;
            else
                h = width / (n+1);
            end
            e = ones(n,1);
            L = spdiags([e -2*e e], -1:1, n, n) / h^2;
            if periodic
                L(end,1) = 1/h^2;
                L(1,end) = 1/h^2;
            end
        end

    end
end