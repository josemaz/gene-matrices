# https://sudonull.com/post/7909-Julia-Scripts-and-parsing-command-line-arguments
# import Pkg; Pkg.add("PyPlot")
using DelimitedFiles, SharedArrays, Random, LinearAlgebra, Clustering, Base
using DataFrames, CSV
using Statistics, ArgParse

nmats = 100

function max_eigv(gene) #returns largest eigenvalue of a randomized matrix generated from data.
    gran = zeros(size(gene))
    for i = 1:size(gene,2)
        gran[:,i] = gene[shuffle(1:end),i]
    end
 #    @time mr = corr_mat(gran)
    # println(mr[[end-1],end])
    mr = cor(transpose(gran))
    # println(mr[[end-1],end])
    return eigvals(mr)
end

parsed_args = Dict([("c", true), ("expfile", "Data/clean/Basal-clean.tsv")])

##============ 0.5 - Input
#! gname s1 s2 ... gstart chromname
df = CSV.read(parsed_args["expfile"], delim='\t')

df = df[df.chromname .== "17", :]
gene = convert(Array{Float64,2},df[:, 2:end-2])

evs = SharedArray{Float64, 2}((nmats,size(gene)[1]))
@time for i in 1:nmats 
    @time evs[i,:] = max_eigv(gene) # return array of ngenes size
end
# All randomized eigenvalues
new_evs = vec(evs)
new_evs = filter(x->x>0.0000001, new_evs) #Drop Zeros

# Original eigenvalues
mr = cor(transpose(gene))
mr_e = eigvals(mr)
mr_e = filter(x->x>0.0000001, mr_e) #Drop Zeros

using PyPlot

#Figura 2a
plt.clf()
# color = '0.75'
# color = '#eeefff'
plt.hist(new_evs,density=true,bins=100,label="Randomized",color=[:C1])
# grid("on")
xlabel("Eigenvalues")
ylabel("Probability density")
legend()
ax = gca()
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
plot(size=(1200,800))
savefig("Plots/Basal-ch17-eigendist.png")


#Figura 2b
plt.clf()
plt.hist(mr_e,density=true,bins=100,label="Empirical data")
plt.hist(new_evs,density=true,bins=100,alpha=0.7,label="Randomized")
# grid("on")
xlabel("Eigenvalues")
ylabel("Probability density")
legend()
ax = gca()
ax.spines["top"].set_visible(false)
ax.spines["right"].set_visible(false)
plot(size=(1200,800))
savefig("Plots/Basal-ch17-eigenRandEmpirical.png")
