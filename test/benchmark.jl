

reload("EMST.jl")
using EMST

function f_benchmark( d , n , reps ; nmin=64)
    ti = Vector{Float64}()
    for zi=1:reps
        tic()
        x_sol = EMST.compute_emst(rand(d,n);nmin=nmin)
        time_i = toc()
        push!(ti,time_i)
    end
    return ( mean(ti) , ti )
end


function test_optimal_nmin( d::Vector{Int64} , n::Vector{Int64} , nmins::Vector{Int64} , reps )
    # results[i,j][k] contains results for d[i] , n[j], and min=nmins[k]
    results      = Array{Vector{Float64},2}(length(d),length(n))
    results_full = Array{Vector{Vector{Float64}},2}(length(d),length(n))

    for zi in 1:length(d)
        for zj in 1:length(n)
            results_ij = Vector{Float64}(length(nmins))
            full_results_ij = Vector{Vector{Float64}}(length(nmins))
            for zk in 1:length(nmins)
                (t_m,t_all) = f_benchmark( d[zi] , n[zj] , reps ; nmin=nmins[zk])
                results_ij[zk] = t_m
                full_results_ij[zk] = t_all
            end
            results[zi,zj]      = results_ij
            results_full[zi,zj] = full_results_ij
        end
    end
    return ( results , results_full )
end




(results , results_full) = test_optimal_nmin( [2,4,20] , [8000] , [48,64,80,96,128,256,512] , 4)


(results , results_full) = test_optimal_nmin( [2,4,20] , [1000,2000,40000,8000,16000] , [48,64,80,96,128,256,512] , 4)



results[2,2]
