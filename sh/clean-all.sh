
cd Data
read -p "WARNING: You want delete Data (y/n)?" choice
case "$choice" in 
  y|Y ) echo "Delete"		
		rm -rf Expression 
		rm -rf Clustered
		rm -rf Pearson
		;;
  n|N ) echo "no"
		;;
  * ) echo "invalid"
		;;
esac


cd Plots
read -p "WARNING: You want delete Plots (y/n)?" choice
case "$choice" in 
  y|Y ) echo "rm -rf Plots"
		# rm -rf *
		;;
  n|N ) echo "no"
		;;
  * ) echo "invalid"
		;;
esac
