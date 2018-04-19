for dir in $*
do
while read p
do
canceljob $p
done <$dir/jobs.log
done
