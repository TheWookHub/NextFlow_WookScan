#!/bin/bash

set -euo pipefail

# if [[ -d $PWD/TrimmedReads ]]
# then
#     echo "$PWD/TrimmedReads directory exists!"
#     echo ""
# else
#     mkdir $PWD/TrimmedReads
#     echo "Trimmed reads directory created at $PWD/TrimmedReads"
#     echo ""
# fi

preTrimFile=!{file_path}
myPrefix=`echo $preTrimFile | cut -d "." -f 1`
newFile=${myPrefix}_trimmed.fastq.gz
report_json=${myPrefix}_trimmed.json
report_html=${myPrefix}_trimmed.html
# fastp -t 100 -i "$FILE" -z 9 -o "$newFile" -R "$myPrefix" -j "$report_json" -h "$report_html"
fastp -t 100 -i "$preTrimFile"  -z 9 -o "$newFile" -R "$myPrefix" -j "$report_json" -h "$report_html"



# for FILE in *;
#     do
#         if [[ "$FILE" =~ .fastq.gz ]]
#         then
#             myPrefix=`echo $FILE | cut -d "." -f 1`
#             newFile=TrimmedReads/${myPrefix}_trimmed.fastq.gz
#             report_json=TrimmedReads/${myPrefix}_trimmed.json
#             report_html=TrimmedReads/${myPrefix}_trimmed.html
#             # echo $newFile
#             # echo $report_html
#             # echo $report_json
#             # echo $myPrefix
#             # fastp -t 100 -i "$FILE" -z 9 -o "$newFile" -R "$myPrefix" -j "$report_json" -h "$report_html"
# 	      fastp -t "$1" -i "$FILE" -z 9 -o "$newFile" -R "$myPrefix" -j "$report_json" -h "$report_html"
#         fi
# done;
