# split -b 30m gene-expression.tgz data.

cat data.* > gene-expression.tgz
tar xzvf gene-expression.tgz
cd Data/Expression/
md5sum -c checksums
[[ $? == 0 ]] \
   && echo "DATA INSTALLED SUCCESSFULLY" \
   && cd ../.. && rm -rf gene-expression.tgz && exit 0
echo "Something was wrong"
exit 15


# clean
# rm -rf Data gene-expression.tgz
