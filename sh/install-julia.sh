cd $HOME
wget "https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz"
tar xf julia-1.4.2-linux-x86_64.tar.gz
echo "export PATH=$HOME/julia-1.4.2/bin:$PATH" >> .bashrc
source .bashrc
julia -e 'import Pkg; Pkg.add("ProgressMeter")'
julia -e 'import Pkg; Pkg.add("Clustering")'
julia -e 'import Pkg; Pkg.add("DataFrames")'
julia -e 'import Pkg; Pkg.add("CSV")'
julia -e 'import Pkg; Pkg.add("ArgParse")'
julia -e 'import Pkg; Pkg.add("Plots")'
julia -e 'import Pkg; Pkg.add("Gadfly")'
julia -e 'import Pkg; Pkg.add("GRUtils")'
julia -e 'import Pkg; Pkg.add("PyPlot")'
julia -e 'import Pkg; Pkg.add("GR")'

rm -rf julia-1.4.2-linux-x86_64.tar.gz*
echo "To remove julia: rm -rf $HOME/julia-1.4.2"
echo "JULIA INSTALLED SUCCESSFULLY"
echo "Quit your session please..."
