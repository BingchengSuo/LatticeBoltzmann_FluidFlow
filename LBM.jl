using Plots
function main()

    function distance(x1, y1, x2, y2)
        return sqrt((x2-x1)^2 + (y2-y1)^2);
    end;

    # define grid and time constants
    Nx = 400;   # x dimension 
    Ny = 100;   # y dimension 
    Nt = 3000;  # iteration time 
    τ  = 0.53;  # collision time scale / kinematic viscosity

    # defind lattice constants 
    NL = 9; # 9 different velocities for adjacent grids
    cxs = ([0, 0, 1, 1, 1, 0, -1, -1, -1]); # x coordinates of the nine adjacent grids
    cys = ([0, 1, 1, 0, -1, -1, -1, 0, 1]); # y coordinates of the nine adjacent grids
    weights = ([4/9, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36, 1/9, 1/36]);

    # define initial conditions
    F = ones(Ny, Nx, NL) + 0.01 * randn((Ny, Nx, NL)); 
    F[:,:,4] .= 2.3; # set a right velocity by assigning additional velocities to the third node 

    # obstacle
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
        for (i, cx, cy) in zip(1:NL, cxs, cys)
            F[:,:,i] = circshift(F[:,:,i], (0, cx));
            F[:,:,i] = circshift(F[:,:,i], (cy, 0));
        end;


        # boundary
        bndryF = F[cylinder, :];
        bndryF = bndryF[:, [1, 6, 7, 8, 9, 2, 3, 4, 5]];
        

        # fluid variables
        ρ = sum(F, dims=3)[:,:,1];
        tempx = zeros(Ny, Nx);
        tempy = zeros(Ny, Nx);
        for i=1:NL
            tempx += (F[:,:,i] * cxs[i]);
            tempy += (F[:,:,i] * cys[i]);
        end
        ux = tempx ./ ρ;
        uy = tempy ./ ρ;
    
        F[cylinder, :] = bndryF;
        ux[cylinder] .= 0;
        uy[cylinder] .= 0;

        # collision
        F_eq = zeros(size(F)); # equilibrium
        for (i, cx, cy, w) in zip(1:NL, cxs, cys, weights)
            F_eq[:,:,i] = ρ .* w .* (1 .+ (3 .* (cx.*ux .+ cy.*uy)) .+ ((9/2) .* (cx.*ux .+ cy.*uy).^2) .- ((3/2) .* (ux.^2 .+ uy.^2)));
        end;
        F += (-1/τ) .* (F .- F_eq);
        if it%100 == 0
            display(heatmap((uy.^2 + ux.^2).^0.5, aspect_ratio=1));
        end;
    end;
end;
