#!/bin/bash

for i in Clustered Pearson
do 
	echo "Downloading in Data/${i}"
	mkdir -p Data/${i}
	cd Data/${i}
	url="http://espinal3.inmegen.gob.mx:5270/genecorr/${i}"
	wget -c -o down.log -r -nH --cut-dirs=2 --no-parent --reject="index.html*" ${url}
	[[ $? -ne 0 ]] && echo "Download error: check down.log" && exit 15
	rm -rf ${i}
	cd ../..
done

