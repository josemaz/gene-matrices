function get_nnd(s)
    nnd = []
    clusters = unique(s)
    for c in clusters
        g_c = findall(x-> x==c, s)
        forwd = diff(g_c)
        backwd = map(x->abs(x),diff(reverse(g_c)))
        nnd_s = map((x,y)-> minimum([x,y]),forwd,backwd)
        push!(nnd,nnd_s) 
    end
    return vcat(nnd...)
end
function get_H(s::Array{Int64,1})
    T = Dict{Any,Float64}()
    n = length(s)
    for i = 1:length(s)
        T[s[i]] = get(T, s[i],0) + 1/n
    end
    p = collect(values(T))
    return mapreduce(x-> -x*log(x),+,p)
end
function get_freq(s)
    T = Dict{Int64,Int64}()
    n = length(s)
    for i = 1:n
        T[s[i]] = get(T, s[i],0) + 1
    end
    return T
end
function get_nnd_distro(s,dmax::Int64)
    P = zeros(dmax)
    n = length(s)
    for i in 1:dmax
        ndata = findall(x-> x==i, s)
        if !isempty(ndata)
            P[i] = (length(ndata) + 1/dmax)/(n + 1/dmax)
        else
            P[i] = (1/dmax) / (n+1/dmax)
        end
    end
    return P
end
KS_d(x,y) = maximum(abs.(cumsum(x) - cumsum(y)))
function get_gene_cor(gene_exp)
    pcor = Float64[]
    g_di = Float64[]
    lg = size(gene_exp)[1]
    for i in 1:lg-1
        for j in i+1:lg
            #try
            push!(pcor, cor(gene_exp[i,2:end-1], gene_exp[j,2:end-1]))
            push!(g_di, abs(gene_exp[i,end]- gene_exp[j,end]))
        end
    end
    g_ord = sortperm(g_di)
    return [g_di[g_ord] pcor[g_ord]]
end
function get_gene_null(gene_exp)
    pcor = Float64[]
    g_di = Float64[]
    lg = size(gene_exp)[1]
    for i in 1:lg-1
        for j in i+1:lg
            #try
            push!(pcor, cor(shuffle(gene_exp[i,2:end-1]), gene_exp[j,2:end-1]))
            push!(g_di, abs(gene_exp[i,end]- gene_exp[j,end]))
        end
    end
    g_ord = sortperm(g_di)
    return [g_di[g_ord] pcor[g_ord]]
end
function get_smooth_cor(gene_cor; num_w=1000, q_1=0.25, q_2=0.75)
    w_s = round((gene_cor[:,1][end] - gene_cor[:,1][1] )/ num_w)
    med_dis = Float64[]
    med_cor = Float64[]
    q_cor = []
    for nw in 0:num_w-1
        ini = nw*w_s
        fin = nw*w_s + w_s
        ix_dw = findall(x-> ini<x<fin, gene_cor[:,1])
        if isempty(ix_dw); continue; end
        push!(med_dis, median(gene_cor[ix_dw,1]))
        m_cor = median(gene_cor[ix_dw,2])
        push!(med_cor, m_cor)
        qs = quantile(gene_cor[ix_dw,2],[q_1,q_2])
        push!(q_cor,[m_cor-qs[1] qs[2]-m_cor])
    end
    q_cor =vcat(q_cor...)
    return med_dis, med_cor, q_cor
end
moving_average(vs,n) = [sum(@view vs[i:(i+n-1)])/n for i in 1:(length(vs)-(n-1))]
function get_plot_cor(gene_cor,gene_null)
    dis_cor, val_cor, q_cor = get_smooth_cor(gene_cor)
    dis_null, val_null, q_null = get_smooth_cor(gene_null)
    mv_w = 11
    dis_cor_mv = moving_average(dis_cor,mv_w)
    val_cor_mv = moving_average(val_cor,mv_w)
    q_cor_mv = [moving_average(q_cor[:,1],mv_w) moving_average(q_cor[:,2],mv_w)]
    dis_null_mv = moving_average(dis_null,mv_w)
    val_null_mv = moving_average(val_null,mv_w)
    q_null_mv = [moving_average(q_null[:,1],mv_w) moving_average(q_null[:,2],mv_w)]
    pt = plot(leg=false, #title="$gtitle",
        size=(1100,1000)
    )
    scatter!(pt,gene_cor[:,1]/10^8,gene_cor[:,2],
        alpha=0.08,
        ylims=(-0.3,1.0),
        color="#36983B"
    )
    scatter!(pt,gene_null[:,1]/10^8,gene_null[:,2], 
        alpha=0.08,
        color="#EA5108"
    )
    plot!(pt,dis_cor_mv/10^8, val_cor_mv,ribbon=q_cor_mv, color="#36983B",linewidth=2)
    plot!(pt,dis_null_mv/10^8, val_null_mv,ribbon=q_null_mv, color="#EA5108",linewidth=2)
    return pt
end
function get_genes_ch(gene_tags, genes_data)
    g_ch = []
    g_st = []
    ix_x = []
    genes= []
    for g in gene_tags
        ix = findfirst(x-> x==g,genes_data[:,3])
        if isnothing(ix)
            #println("gene $g not found")
            continue
        end
        push!(ix_x, ix)
        push!(genes, g)
        push!(g_ch, genes_data[ix,5])
        push!(g_st, genes_data[ix,4])
    end
    return [genes g_st g_ch]
end
