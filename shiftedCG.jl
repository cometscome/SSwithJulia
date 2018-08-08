module Shifted_cg
    export solve
    if v"0.6.9" <= VERSION #The version of Julia is higher than 0.6.9
        using Distributed
        using SparseArrays
    end
    function solve(mat_A,n,b,vec_σ,N,eps=1e-12)
        println("--------------------------------------------------------")
        println("Solving linear equations with the shifted CG method...")
        tb = typeof(b[1]) 
        x = zeros(tb,n)
        r = zeros(tb,n)
        p = zeros(tb,n)
        Ap = zeros(tb,n)
        r[:] = b[:]
        p[:] = b[:]
        αm = 1.0
        βm = 0.0
        tσ = typeof(vec_σ[1])
        vec_x = zeros(tσ,n,N)
        vec_r = zeros(tσ,n,N)
        vec_p = zeros(tσ,n,N)
        ρm = ones(tσ,N)
        ρ0 = ones(tσ,N)
        ρp = ones(tσ,N)
        for j in 1:N
            vec_r[:,j] = b[:] 
            vec_p[:,j] = b[:] 
        end

        k = 0
        hi = 1.0
#        println(eps)
#        println(hi)
        while abs(hi) > eps
            A_mul_B!(Ap,mat_A,p)
            #Ap = mat_A*p
            pAp = p'*Ap
            rr = r'*r
            αk = rr/pAp
            x += αk*p
            r += -αk*Ap
            βk = r'*r/rr
            p = r + βk*p
            for j in 1:N
                ρkj = ρ0[j]
                
                if abs(ρkj) < eps
                    continue
                end
                #ρp[j] = ifelse(abs(ρkj) > eps,ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_σ[j])+αk*βm*(ρkmj-ρkj)),ρkj)
                #αkj = ifelse(abs(ρkj) > eps,(ρp[j]/ρkj)*αk,0.0)
                ρkmj =ρm[j]
                ρp[j] = ρkj*ρkmj*αm/(ρkmj*αm*(1.0+αk*vec_σ[j])+αk*βm*(ρkmj-ρkj))
                αkj = (ρp[j]/ρkj)*αk
                vec_x[:,j] = vec_x[:,j]+αkj*vec_p[:,j]
                βkj = (ρp[j]/ρkj)^2*βk
                vec_p[:,j] = ρp[j]*r+βkj*vec_p[:,j]
                #vec_p[:,j] = ifelse(abs(ρkj) > eps,ρp[j]*r+βkj*vec_p[:,j],vec_p[:,j])

            end

            ρm[:] = ρ0[:]
            ρ0[:] = ρp[:]
            αm = αk
            βm = βk
    
            ρMAX = maximum(abs.(ρp))^2
            hi = abs(rr*ρMAX)
            #println(hi," ",k)
            k +=1
            if k > 10*n
                println(hi," ",k)
                println("not converged! there might be zero eigenvalue")
                conv = false
                return vec_x,x,conv
            end
        end
        println("done.")
        println("error = ",hi)
        println("Num. of iterations = ", k)
        conv = true

        println("--------------------------------------------------------")
        

        return vec_x,x,conv
    end



end