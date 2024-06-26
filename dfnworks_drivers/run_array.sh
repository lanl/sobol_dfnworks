#!/bin/bash

num_of_total_tries=16
N=8
i=1
str1="/lclscratch/murph/outfiles/dfn_outfile_"
str2=.txt

str3="/lclscratch/murph/compiled_data_files/output"
str4="/lclscratch/murph/output"
str5="/num_of_tries.txt"
filestr1="/SA_params.txt"
filestr2="/dfnTrans_output"
filestr3="/parsed_vtk"
totalexps=8
twenty=20

(
for exp in `seq $num_of_total_tries`;
do

	exp2=$((exp-1))
	exp2=$((exp2%totalexps))
	((i=i%N)); ((i++==0)) && wait
	if [ $i == 1 ]; then
                echo "Moving files and removing extra logs"
                for dec in `seq $N`;
                do
			exp1=$((exp2-dec+1))
				
			if [ -d "$str4$exp1" ]; then
                        	fullstr="$str3$exp1"
				if [ ! -d $fullstr ]; then
                        		mkdir $fullstr
				fi
                        	filetxt="$str4$exp1$filestr1"
                       		mv $filetxt $fullstr
				echo "moved $filetxt to $fullstr"
                        	filetxt2="$str4$exp1$filestr2"
                        	mv $filetxt2 $fullstr
				echo "moved $filetxt2 to $fullstr"
                        	filetxt3="$str4$exp1$filestr3"
                        	mv $filetxt3 $fullstr
				echo "moved $filetxt3 to $fullstr"
                        	filetxt4="$str4$exp1"
                        	rm -r $filetxt4
			else
				echo "output directory $str4$exp1 does not exist, so it was not moved."
			fi

                done
        fi

	if  [ ! -f "$str3$exp2$filestr1" ] || [ ! -d "$str3$exp2$filestr2" ]; then	
		echo "Beginning job $exp2"
		echo "This is attempt number"
		expr $exp / $twenty
		fullstr="$str3$exp2"
		if [ ! -d "$fullstr" ]; then
			mkdir $fullstr
		fi

		attemptstr=$str3$exp2$str5
		expr $exp / $twenty > $attemptstr
		fullstr="$str1$exp2$str2"
		python3.8 TPL_sim_experiment.py $exp2 > $fullstr &
	else
		echo "Job $exp2 has finished, so we will skip it."
	fi
done
)
