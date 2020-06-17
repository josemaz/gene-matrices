# https://sudonull.com/post/7909-Julia-Scripts-and-parsing-command-line-arguments

using DelimitedFiles, SharedArrays, Random, LinearAlgebra, Clustering, Base
using DataFrames, CSV
using Statistics, ArgParse
using ProgressMeter

nmats = 100 # iteraciones para samplear los eigenvalues de las correlaciones
n_kmedoids = 100 # iteraciones para calcular el minimo de kmedoids
n_kmedoids_out = 10 # Output de de las iteraciones del minimo de los kmedoids


############# FUNCTIONS
function max_eigv(gene) 
    #returns largest eigenvalue of a randomized matrix generated from data.
    gran = zeros(size(gene))
    for i = 1:size(gene,2)
        gran[:,i] = gene[shuffle(1:end),i]
    end
 #    @time mr = corr_mat(gran)
    # println(mr[[end-1],end])
    mr = cor(transpose(gran))
    # println(mr[[end-1],end])
    return maximum(eigvals(mr))
end
#------------------------


############### MAIN
## JULIA_NUM_THREADS=$(nproc)
s = ArgParseSettings()
@add_arg_table! s begin
    "-c"
        help = "Flag to calculate number of clusters"
        action = :store_true
    "expfile"
        help = "Matriz de expresion (read a clean format exp matrix)"
        arg_type = String
        required = true
end
parsed_args = parse_args(s)
# parsed_args = Dict([("c", true), ("expfile", "Data/clean/Basal-clean.tsv")])
println(parsed_args)
fname = basename(parsed_args["expfile"])
subtype = split(fname,'-')[1]
println("Working on " * subtype)

println("Kamedoids clustering method on Cancer Expression Matrix")
flush(stdout)

Random.seed!(34568)

### INPUT FILES
#! gname s1 s2 ... gstart chromname
expr = readdlm(parsed_args["expfile"],skipstart=1)
#! cut cols: gname, gstart, chromname
gene = convert(Array{Float64,2},expr[:, 2:end-2])


###################
println("First correlation (m)")
@time  m = cor(transpose(gene))
gnames = expr[:,1]
gstarts = expr[:,end-1]
chroms= expr[:,end]
df = DataFrame(src = [], dst = [], src_gstart = Int64[],
    dst_gstart = Int64[], src_chrom = [], dst_chrom = [],
    pearson = Float64[])
@showprogress for i in 1:length(gnames)
    for j in i+1:length(gnames)
        # println(i," ",j," ",gnames[i]," ",gnames[j]," ",
        #     gstarts[i]," ",gstarts[j]," ",chroms[i]," ",
        #     chroms[j]," ", m[i,j])
        row = [gnames[i] gnames[j] gstarts[i] gstarts[j] chroms[i] chroms[j] m[i,j]]
        push!(df, row)
    end
end
odir = "Data/Pearson/"
mkpath(odir)
oname = odir * subtype * "-allchroms.tsv"
CSV.write(oname,df)


###################
println("Automatic number of clusters detection for all chromosomes")
if parsed_args["c"] == true
    println("Get Eigenvalues from shuffles")
    l_evmax = SharedArray{Float64, 1}((nmats))
    for i in 1:nmats #doing the surrogate matrices and saving max eigv.
        @time l_evmax[i] = max_eigv(gene)
        println("Iteration: ",i,", max eigenvalue: ",l_evmax[i])
        flush(stdout)
    end
    #! get eigenvalue of gene
    m_ev = filter(x->x>0.0000001, eigvals(m))
    #! Todos los eigenvalues de la correlacion que sean mayores al mayor de las 
    #! correlaciones muestreadas; el numero de estos eigenvalues es el numero de cluesters
    #! nc = 13, for Basal-All.txt with nmats=100
    nc = length(filter(x->x>maximum(l_evmax),m_ev))
    println("Number of clusters (nc): ",nc)
end


###################
println("Starting Kmedoids sampling for all chromosomes")
# make m in [0,1] insted [-1,1]
# correalacion solo en valores positivos
dmat = map(x-> 1-abs(x),m)

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


###################
#! Saving
# writedlm("out.tsv", [expr[:,1] c.assignments c.acosts expr[:,end-1] expr[:,end]])
# header = ["gname" "clusterid" "assignmentcost" "gstart" "chromosome"]
# data = [expr[:,1] c.assignments c.costs expr[:,end-1] expr[:,end]]
# writedlm("out.tsv", [header ; data])
df = DataFrame([expr[:,1],c.assignments,c.costs,expr[:,end-1],expr[:,end]])
header = ["gname","clusterid","assignmentcost","gstart","chromosome"]
rename!(df, Symbol.(header)) # Force set header names
odir = "Data/Clustered/"
mkpath(odir)
oname = odir * subtype * "-all-clusters.tsv"
CSV.write(oname,df)


