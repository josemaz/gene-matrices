# split -b 30m gene-expression.tgz data.

cat data.* > gene-expression.tgz
tar xzvf gene-expression.tgz
cd Data/Expression/
md5sum -c checksums.md5
echo "DATA INSTALLED SUCCESSFULLY"

# clean
# rm -rf Data gene-expression.tgz
