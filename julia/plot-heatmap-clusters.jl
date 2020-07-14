using CSV, DataFrames, Plots

clusfile = "Data/Clustered/Healthy/Healthy-chr1-clusters.tsv"
print(clusfile)
df = CSV.read(clusfile, delim=',', copycols=true)
print(first(df,5))

nclusters = maximum(df[!,"clusterid"])
sort!(df, [:chromosome, :clusterid])

df = df[1:5,:]
M = Array{Int64}(zeros(nrow(df), nrow(df)))
genes = df[!,:gname]
cl = df[:,:clusterid]
for i in 1:length(genes)
	clidsrc = cl[i]
	M[i,i] = cl[i]
	for j in i+1:length(genes)
		# println("$i $j")
		cliddst = cl[j]
		clidsrc == cliddst ? M[i,j] = clidsrc : M[i,j] = 0
	end
end

# print(M)
# myplot = spy(M)
# draw(PNG("Plots/myplot.png", 3inch, 3inch), myplot)
heatmap(M, title="Threshold")
plot(size=(1200,1200))

savefig("Plots/clustering.png")