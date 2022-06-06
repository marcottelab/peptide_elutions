i=1;
while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    header="$line"
    echo "$header" >> seq"${i}".fasta
  else
    seq="$line"
    echo "$seq" >> seq"${i}".fasta
    ((i++))
  fi
done < $1


