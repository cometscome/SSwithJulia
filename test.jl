@everywhere include("SS_p.jl")
@everywhere import SS



@everywhere function make_1Dtb(N,periodic=true) #1D tight binding model
    mat_H = spzeros(N,N)
    t = 1.0
    
    for i in 1:N
        for dx in -1:1
            j = i + dx
            if periodic
                if j > N
                    j = j -N
                elseif j < 1
                    j = j +N
                end
            end
            
            if 1 <= j <= N
                if dx == 1
                    mat_H[i,j] = -t
                elseif dx == -1
                    mat_H[i,j] = -t
                elseif dx ==0
                    mat_H[i,j] = 1.5
                end
                    
            end
            
        end
    end
    
    return mat_H
end


@everywhere function make_2Dtb(Nx,Ny,periodic=true) #2D tight binding model
    mat_H = spzeros(Nx*Ny,Nx*Ny)
    t = 1.0
    μ = -1.5
    
    for ix in 1:Nx
        for iy in 1:Nx
            for dx in -1:1
                for dy in -1:1
                    jx = ix + dx
                    jy = iy + dy
                    v = 0.0
                    if periodic
                        if jx > Nx
                            jx = jx - Nx
                        elseif jx < 1
                            jx = jx + Nx
                        end
                        if jy > Ny
                            jy = jy - Ny
                        elseif jy < 1
                            jy = jy + Ny
                        end                        
                    end
                    if 1 <= jx <= Nx && 1 <= jy <= Ny
                        if dx == 0 
                            if dy == 0
                                v = -μ
                            elseif dy == 1
                                v = -t
                            elseif dy == -1
                                v = -t
                            end
                        elseif dx == 1                            
                            if dy == 0
                                v = -t
                            end
                        elseif dx == -1                            
                            if dy == 0
                                v = -t
                            end
                        end
                        i = (iy-1)*Nx+ix
                        j = (jy-1)*Ny+jx
                        mat_H[i,j] = v
                    end
                end
            end
        end
    end

    
    
    return mat_H
end


@everywhere function make_2Dtbsc(Nx,Ny,periodic=true) #2D tight binding model with superconducting SNS π junction
    mat_H = spzeros(Nx*Ny*2,Nx*Ny*2)
    N = Nx*Ny
    t = 1.0
    μ = -1.5
    Δ = 1.0
    
    for ix in 1:Nx
        for iy in 1:Nx
            for dx in -1:1
                for dy in -1:1
                    jx = ix + dx
                    jy = iy + dy
                    v = 0.0
                    if periodic
                        if jx > Nx
                            jx = jx - Nx
                        elseif jx < 1
                            jx = jx + Nx
                        end
                        if jy > Ny
                            jy = jy - Ny
                        elseif jy < 1
                            jy = jy + Ny
                        end                        
                    end
                    i = (iy-1)*Nx+ix
                    j = (jy-1)*Ny+jx
                    if 1 <= jx <= Nx && 1 <= jy <= Ny
                        if dx == 0 
                            if dy == 0
                                v = -μ                                
                                vd = ifelse(jx < Nx/2,Δ,-Δ)
                                
                                if Nx/2-Nx/6 <= jx <= Nx/2 + Nx/6
                                  vd = 0.0  
                                end
                                mat_H[i,j+N] = vd
                                mat_H[i+N,j] = vd                                   
                                    
                            elseif dy == 1
                                v = -t
                            elseif dy == -1
                                v = -t
                            end
                        elseif dx == 1                            
                            if dy == 0
                                v = -t
                            end
                        elseif dx == -1                            
                            if dy == 0
                                v = -t
                            end
                        end

                        mat_H[i,j] = v
                        mat_H[i+N,j+N] = -v
                    end
                end
            end
        end
    end

    
    
    return mat_H
end




@everywhere function main()
    Nx = 12
    Ny = 12

    
    ρ = 0.2
    γ = 0.0
    ε =0
#    mat_H = make_1Dtb(N) #1D system
    N = Nx*Ny
#    mat_H = make_2Dtb(Nx,Ny,false) #2D system
    N = Nx*Ny*2
    mat_H = make_2Dtbsc(Nx,Ny,false) #2D system with superconducting SNS π-junction
    
    println("Dimension of the matrix: ",N)
    rε = Float64[]
    if N <= 8192
        println("Doing the full diagonalization...")
        @time ε,vec_w = eig(full(mat_H))
        println("Done")
        

        integers = Int64[]
        for i in 1:N
            push!(integers,i)
        end      
        
        
        rε = Float64[]
        is = 0
        for i in 1:N
            if γ-ρ <= ε[i] <= γ+ρ
                is += 1
                push!(rε,ε[i])
            end
        end         
    end
   

   

    #println(ε)
    if N < 2000
      
        for i in 1:is
            println(i, "       ",rε[i])
        end
       
        
    end
    

    
    #println(typeof(mat_H[1,1]))
    println("Doing the SS method...")
    @time eigenvalues,residuals,eigenvectors,num = SS.eigensystem(mat_H,N,ρ,γ)
    println("done.")
    
    if N <= 8192
        println("number    eigenvalue    residual    original     difference")
        for i in 1:num
            println(i,"     ",eigenvalues[i],"    ",residuals[i],"     ",rε[i], "    ",eigenvalues[i]-rε[i])
        end
    else
        println("number    eigenvalue    residual    ")
        for i in 1:num
            println(i,"     ",eigenvalues[i],"    ",residuals[i])
        end        
    end
    println("-----------------------------------------") 
#plot(integers[:],ε[:],label="Original") 
#    plot(integers[1:num],eigenvalues[1:num],label="SS method") 

#    plot(integers[1:N],ε[1:N],label="Original")      
    
    return rε,eigenvalues,residuals,eigenvectors,num

    
end    




@time rε,eigenvalues,residuals,eigenvectors,num = main()



println("Eigenvalues are calculated!")


