#!/bin/bash
 
 main=memFA20mem
 link=foo
 subsh=runqsub
 init=1019

# ifcst=("24" "48" "72" "96" "120" "144" "168" "192"\
#           "216" "240")
  ifcst=("00")
  echo "Now the matlab is:      ${main}.m"
  echo "Now the foo is:      ${link}.txt"
    cp ${main}.m copy_${main}.m
    cp ${link}.txt copy_${link}.txt
    cp ${subsh}.sh copy_${subsh}.sh

  for ((xx=1;xx<20;xx=xx+1)) ; do
    cp copy_${main}.m ${main}_${xx}.m
    cp copy_${link}.txt ${link}_${xx}.txt
    cp copy_${subsh}.sh ${subsh}_${xx}.sh
    echo "Now the xx is: ${xx}"
    sed -i "s#memFA20mem#${main}_${xx}#g" ${link}_${xx}.txt
#set initial date
    sed -i "s#date=1001#date=${init}#g" ${main}_${xx}.m
#set forecast lead time
#    sed -i "s#fcst=120#fcst=${ifcst[$xx]}#g" ${main}_${xx}.m
#    sed -i "s#fcst=120#fcst=00#g" ${main}_${xx}.m
#read in file
    sed -i "s#foo#${link}_${xx}#g" ${subsh}_${xx}.sh
#define output log file
    sed -i "s#out#out_${xx}#g" ${subsh}_${xx}.sh
#define task name
    sed -i "s#memFA#memFA${xx}#g" ${subsh}_${xx}.sh
#cp sub.sh sub$xx.sh
#sed -i "s#a.out#${xx}.out#g" sub$xx.sh
    qsub ${subsh}_${xx}.sh
    echo "qsub ok"
  done
