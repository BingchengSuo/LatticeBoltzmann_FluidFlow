using Plots
function main()
    function distance(x1, y1, x2, y2)
        return ((x2-x1)^2 + (y2-y1)^2)^0.5;
    end;


    # define grid and time constants
    Nx = 400;   # x dimension 
    Ny = 100;   # y dimension 
    Nt = 50;  # iteration time 
    τ  = 0.53;  # collision time scale / kinematic viscosity

    # defind lattice constants 
    NL = 9; # 9 different velocities for adjacent grids
    cxs = ([0, 0, 1, 1, 1, 0, -1, -1, -1]); # x coordinates of the nine adjacent grids
    cys = ([0, 1, 1, 0, -1, -1, -1, 0, 1]); # y coordinates of the nine adjacent grids
    weights = ([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36]);

    # define initial conditions
    F = ones(Ny, Nx, NL) + 0.01 * randn((Ny, Nx, NL)); 
    F[:,:,3] .= 2.3; # set a right velocity by assigning additional velocities to the third node 

    # obstacle
    grids = Array{Int64}(undef, Ny, Nx);
    for y = 1:Ny
        for x = 1:Nx
            if distance(Nx//4, Ny//2, x, y) < 13
                grids[y,x] = 1;
            end;
        end;
    end;

    cylinder = fill(false, (Ny, Nx));


    for y=1:Ny
        for x=1:Nx
            if distance(Nx//4, Ny//2, x, y) < 13.0
                cylinder[y,x] = true;
            end;
        end;
    end;
    

    # main loop
    for it = 1:Nt

        # streaming
        for i = 1:NL
            for (cx, cy) in zip(cxs, cys)
                F[:,:,i] .= circshift(F[:,:,i], cx);
                F[:,:,i] .= circshift(F[:,:,i], cy);
            end;
        end;

        # boundary
        bndryF = F[cylinder, :];
        bndryF = bndryF[:, [1, 6, 7, 8, 9, 2, 3, 4, 5]];

        # fluid variables
        ρ = sum(F, dims=3);
        temp_F_x = zeros(Ny, Nx, NL);
        temp_F_y = zeros(Ny, Nx, NL);

        for y = 1:Ny # I cannot find convenient way to do 3_d product 
            for x = 1:Nx
                temp_F_x[y,x,:] = F[y,x,:] .* cxs;
                temp_F_y[y,x,:] = F[y,x,:] .* cys;
            end;
        end;
        ux = sum(temp_F_x, dims=3) ./ ρ;
        uy = sum(temp_F_y, dims=3) ./ ρ;

        F[cylinder, :] = bndryF;
        ux[cylinder,:] .= 0;
        uy[cylinder,:] .= 0;


        # collision
        F_eq = zeros(size(F)); # equilibrium
        for i = 1:NL
            for (cx, cy, w) in zip(cxs, cys, weights)
                F_eq = ρ .* w .* (1 .+ 3 .* (cx.*ux .+ cy.*uy) .+ 9 * (cx.*ux .+ cy.*uy).^2 ./ 2 .- 3 .* (ux.^2 .+ uy.^2)./2);
            end;
        end;
        F = F .+ (-1/τ) .* (F .- F_eq);
        display(heatmap((ux.^2 .+ uy.^2)[:,:,1], c=:diverging_bkr_55_10_c35_n256, aspect_ratio=1));
    end;
end;

