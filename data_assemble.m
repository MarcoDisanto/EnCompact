function X = data_assemble(x, proc_grid)
% function that assembles 3D arrays of variables following the order
% imposed by the cartesian topology. It needs a cell array of 3D arrays and
% the processes grid as input.

    dim = size(proc_grid);
    
    if (size(x)~=dim(1)*dim(2)*dim(3))
        print('error');
    else

        % creating an intermediate cell array to store data assembled along x
        X1 = cell(dim(2), dim(3));

        for k = 1 : dim(3)
            for j = 1 : dim(2)

                i = 1;
                X1{j, k} = x{proc_grid(i, j, k)};

                for i = 2 : dim(1)
                    X1{j, k} = cat(1, X1{j, k}, x{proc_grid(i, j, k)});
                end

            end
        end

        clear x; % no need for this anymore

        % creating an intermediate cell array to store data assembled along y
        X2 = cell(dim(3), 1);

        for k = 1 : dim(3)

            j = 1;
            X2{k} = X1{j, k};

            for j = 2 : dim(2)
                X2{k} = cat(2, X2{k}, X1{j, k});
            end

        end

        clear X1; % no need for this anymore

        k = 1;
        X = X2{k};

        for k = 2 : dim(3)
            X = cat(3, X, X2{k});
        end

        clear X2; % no need for this anymore
    end


end

