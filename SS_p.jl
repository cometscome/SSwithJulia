#Parallel version
@everywhere include("shiftedCG.jl")

module SS
    import shiftedcg
    function eigensystem(mat_H,N,ρ,γ,α=0.1,L0=10,Nq=64,M=10,κ=2,δ=1e-14,numr = 0,cg=true) #for hermitian matrix
        hatV = zeros(Complex128,N,L0)
        hatV = rand([-1.0,1.0],N,L0) #each element should be -1 or 1
        println("Estimating number of eigenvalues...")
        println("------------------------------------------------------------")
    
        γshift = γ
        for i in 1:N
            γ += -γshift
            mat_H[i,i] += γshift
        end        
       
        vec_Sk,conv = calc_SkinEq22(hatV,mat_H,N,ρ,γ,α,L0,Nq,M,cg)
    
        if conv == false
            hatV = rand([-1.0,1.0],N,L0) #each element should be -1 or 1
            eshift = 0.5
            γ += -eshift
            for i in 1:N
                mat_H[i,i] += eshift
            end
        
            vec_Sk,conv = calc_SkinEq22(hatV,mat_H,N,ρ,γ,α,L0,Nq,M,cg)
        
            γ += eshift
            for i in 1:N
                mat_H[i,i] += -eshift
            end        
        end
       
    
        vec_s0 = zeros(Complex128,N,L0)
        vec_s0[1:N,1:L0] = vec_Sk[1:N,1:L0]
        mst = calc_msinEq26(vec_s0,hatV,L0)
    
        println("done.")
        println("------------------------------------------------------------")
        println("Estimated number of eigenvalues: ",mst)
        L = Int64(round(κ*mst/M,RoundNearestTiesAway))
        println("L is ",L)
    
        println("Calculating eigenvalues and eigenvectors...")
        println("------------------------------------------------------------")    

        hatV = zeros(Complex128,N,L)
        hatV = rand(N,L)*2.0-1.0 #each element should be in (-1,1)     
    
        ε,vec_x,ms = calc_SS(hatV,mat_H,N,ρ,γ,α,L,Nq,M,δ,cg)
    

    
        if numr >=1
            for r in 1:numr
                mat_c = rand(ms,L)*2.0-1.0
                hatV = vec_x*mat_c
                ε,vec_x,ms = calc_SS(hatV,mat_H,N,ρ,γ,α,L,Nq,M,δ,cg)            
            end
        end

        ε += -γshift    
        for i in 1:N
            γ += γshift
            mat_H[i,i] += -γshift
        end        
    
    
    
        is = 0
        eigenvalues = Float64[]
        residuals = Float64[]
        eigenvectors = []

        for i in 1:ms
            Hx = mat_H*vec_x[:,i]
            rup = Hx-ε[i]*vec_x[:,i]
            resiup = sqrt(rup'*rup)
            residown1 = sqrt(Hx'*Hx)
            residown2 = sqrt(vec_x[:,i]'*vec_x[:,i])
            resi = real(resiup/(residown1+abs(ε[i])*residown2))
            if resi < 0.1 && γ-ρ <= ε[i] <= ρ+γ
                is += 1
                push!(eigenvalues,ε[i])
                push!(residuals,resi)
                push!(eigenvectors,vec_x[:,i])
            end
        end    
    
    


        #println("-----------------------------------------")   
        println("done.")
        println("------------------------------------------------------------")
        return eigenvalues,residuals,eigenvectors,is
 
    


    
    end

    function calc_SS(hatV,mat_H,N,ρ,γ,α,L,Nq,M,δ,cg)
        vec_Sk,conv = calc_SkinEq22(hatV,mat_H,N,ρ,γ,α,L,Nq,M,cg)
       
        (U,Σ,W ) = svd(vec_Sk) 
        j = 0
        ratio = 1
        while ratio > δ
            j += 1
            ratio = Σ[j]/Σ[1]
            #println(ratio)
            if j == length(Σ)
                ratio = 0.0
            end
        end
        ms = j
        mat_Q = zeros(Complex128,N,ms)
        mat_Q[1:N,1:ms] = U[1:N,1:ms]
    
    
        mat_Ht = mat_Q'*mat_H*mat_Q
        mat_Ht = (mat_Ht'+mat_Ht)/2

        ε,vec_w = eig(mat_Ht)
        vec_x = mat_Q*vec_w    
    
        return ε,vec_x,ms
    end


    function calc_SkinEq22(hatV,mat_H,N,ρ,γ,α,L,Nq,M,cg)
        vec_z = zeros(Complex128,Nq)
        vec_w = zeros(Complex128,Nq)
        vec_Sk = zeros(Complex128,N,L*M)
        conv = true
        for j in 1:Nq
            θj = (2π/Nq)*(j-1/2)
            vec_z[j] = γ + ρ*(cos(θj)+im*α*sin(θj))
            vec_w[j] = α*cos(θj)+im*sin(θj)
        end
    
        
        yarray = pmap(i -> calc_Ski(i,hatV,mat_H,N,ρ,γ,α,L,Nq,M,vec_z,vec_w,cg),1:L)
        #=
        yarray = Array{Any}(L)
        for i in 1:L
            yarray[i] = calc_Ski(i,hatV,mat_H,N,ρ,γ,α,L,Nq,M,vec_z,vec_w,cg)
        end
            =#
    
        for i in 1:L
            vec_y = yarray[i][1]
            conv = yarray[i][2]
            if conv == false
                return vec_Sk,conv
            end
            for k in 0:M-1
                ik = k*L+i
                for j in 1:Nq
                    vec_Sk[:,ik] += ρ*vec_w[j]*vec_z[j]^k*vec_y[:,j]/Nq
                end
            end
        end
    
        return vec_Sk,conv
    end

    function calc_Ski(i,hatV,mat_H,N,ρ,γ,α,L,Nq,M,vec_z,vec_w,cg)
        println("i = ",i,"/",L)
        if cg==false
            vec_y = directsolver(-mat_H,N,hatV[:,i],vec_z,Nq)
            conv = true
        else
            vec_y,y,conv = shiftedcg.solve(-mat_H,N,hatV[:,i],vec_z,Nq) #(zj I - H) Yj = V
        end
        return vec_y,conv

    end

    function directsolver(mat_A,N,b,vec_z,Nq)      
        vec_y = zeros(Complex128,N,Nq)
        mat_X = spzeros(N,N)
        mat_I = speye(N,N)
        for j in 1:Nq
            mat_X[:,:] = mat_A[:,:]
            zj = vec_z[j]
            mat_X = mat_X + zj*mat_I        
            y = mat_X \b
            vec_y[:,j] = y[:]
        
        end
        return vec_y
        
    
    end

    function calc_msinEq26(vec_s0,hatV,L0)
        m = 0.0
        for i in 1:L0
            m += hatV[:,i]'*vec_s0[:,i]
        end
        m = real(m)/L0
        #println(m)
        ms = Int64(round(m))

        return ms
    end

end

