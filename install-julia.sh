cd $HOME
wget "https://julialang-s3.julialang.org/bin/linux/x64/1.4/julia-1.4.2-linux-x86_64.tar.gz"
tar xf julia-1.4.2-linux-x86_64.tar.gz
echo "export PATH=$HOME/julia-1.4.2/bin:$PATH" >> .bashrc
rm -rf julia-1.4.2-linux-x86_64.tar.gz*
echo "To remove julia: rm -rf $HOME/julia-1.4.2"
echo "JULIA INSTALLED SUCCESSFULLY"
echo "Quit your session please..."
