# https://sudonull.com/post/7909-Julia-Scripts-and-parsing-command-line-arguments

using DelimitedFiles, SharedArrays, Random, LinearAlgebra, Clustering, Base
using DataFrames, CSV
using Statistics, ArgParse

nmats = 100 # iteraciones para samplear los eigenvalues de las correlaciones
n_kmedoids = 100 # iteraciones para calcular el minimo de kmedoids
n_kmedoids_out = 10 # Output de de las iteraciones del minimo de los kmedoids



function max_eigv(gene) #returns largest eigenvalue of a randomized matrix generated from data.
    gran = zeros(size(gene))
    for i = 1:size(gene,2)
        gran[:,i] = gene[shuffle(1:end),i] #Shuffle over genes
    end
 #    @time mr = corr_mat(gran)
    # println(mr[[end-1],end])
    mr = cor(transpose(gran))
    # println(mr[[end-1],end])
    return maximum(eigvals(mr))
end



#############
## MAIN
#############

## JULIA_NUM_THREADS=$(nproc)
# julia Julia/02-clusters-chrom.jl -c Data/clean/Healthy-clean.tsv

s = ArgParseSettings()
@add_arg_table! s begin
    # "-c"
    #     help = "Flag to calculate number of clusters"
    #     action = :store_true
    "expfile"
        help = "Expression Gene Matrix (read a clean format exp matrix)"
        arg_type = String
        required = true
end
parsed_args = parse_args(s)
# parsed_args = Dict([("c", true), ("expfile", "Data/clean/Basal-clean.tsv")])
println(parsed_args)
fname = basename(parsed_args["expfile"])
subtype = split(fname,'-')[1]
println("Working on " * subtype)
flush(stdout)

Random.seed!(34568)


##============ 0.5 - Input
#! gname s1 s2 ... gstart chromname
df = CSV.read(parsed_args["expfile"], delim='\t')

#! Selection per chromosome
for g in groupby(df, :chromname)
    cname = unique(g.chromname)[1]
    println("Working on chromosome: ", cname)    
    #! Extract expression by chromosome
    gene = convert(Array{Float64,2},g[:, 2:end-2]) 
    println(size(gene))

    ##============ 01 - First correlation
    println("First correlation (m)")
    @time  m = cor(transpose(gene))
    
    ##============ 02 - Eigenvalues
    println("Get Eigenvalues from ",nmats," shuffles")
    l_evmax = SharedArray{Float64, 1}((nmats))
    @time for i in 1:nmats #doing the surrogate matrices and saving max eigv.
        # @time l_evmax[i] = max_eigv(gene)
        l_evmax[i] = max_eigv(gene)
        # println("Iteration: ",i,", max eigenvalue: ",l_evmax[i])
        flush(stdout)
    end
    #! Get eigenvalue of gene
    m_ev = filter(x->x>0.0000001, eigvals(m)) # Drop Zeros
    #! Todos los eigenvalues de la correlacion que sean mayores al mayor de las 
    #! correlaciones muestreadas; el numero de estos eigenvalues es el numero de cluesters
    #! nc = 13, for Basal-All.txt with nmats=100
    nc = length(filter(x->x>maximum(l_evmax),m_ev))
    println("Number of clusters (nc): ",nc)

    
    ##============ 03 - KMedoids
    println("Starting kmedoids sampling")
    #! Make m in [0,1] insted [-1,1]
    dmat = map(x-> 1-abs(x),m)
    print(size(dmat))
    ens_c = []
    for i = 1:n_kmedoids
        push!(ens_c, kmedoids(dmat,nc))
        if (i % n_kmedoids_out) == 0
          println(i)
        end
    end

    #! De todas los clustereos cual es el que tiene el costo total minimo
    inx_c = findmin(map(x -> x.totalcost, ens_c))[2]
    c = ens_c[inx_c]
    println("Cluster counts: ",c.counts)


    ##============ 04 - Saving
    df = DataFrame([g[:,1],c.assignments,c.costs,g[:,end-1],g[:,end]])
    header = ["gname","clusterid","assignmentcost","gstart","chromosome"]
    rename!(df, Symbol.(header))
    odir = "Data/Clustered/" * subtype * "/"
    mkpath(odir)
    fname = odir*subtype*"-chr"*cname*"-clusters.tsv"
    CSV.write(fname,df)
end



