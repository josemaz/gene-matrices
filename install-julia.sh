cd $HOME
wget "https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz"
tar xf julia-1.4.2-linux-x86_64.tar.gz
export PATH=$HOME/julia-1.4.2/bin:$PATH
echo "JULIA INSTALLED SUCCESSFULLY"
echo "Quit your session please..."
rm -rf julia-1.4.2-linux-x86_64.tar.gz*