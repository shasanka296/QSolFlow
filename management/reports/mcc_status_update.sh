echo
echo "*****************************************************"
echo "**                                                 **"
echo "**               STATUS UPDATE                     **"
echo "**          $(date)           **"
echo "**                                                 **"
echo "*****************************************************"
echo
echo "******************* JOBS RUNNING ******************"
cd /home/rdu230/.status
squeue -u rdu230 | sed 1d | sort -k1 > new_list
squeue -u rdu230

awk '{print $1}' new_list > comp_new
echo
echo "************* JOBID OF COMPLETED JOBS ************"
echo $(date) >> completed_jobs
comm --nocheck-order -23 comp_old comp_new > job_id
grep -f job_id old_list >> completed_jobs
tail -5 completed_jobs
mv comp_new comp_old
mv new_list old_list
echo
echo "******************* NODES IDLE ******************************"
sinfo | grep idle
echo
echo
#echo "******************* NODES MIX ******************************"
#sinfo | grep mix
#echo
squeue -A coa_cmri235_uksr | grep ' R ' > running
squeue -A coa_cmri235_uksr | grep ' PD ' > pending
usr=$(squeue -A coa_cmri235_uksr | awk '{print$4}'| awk '!seen[$0]++' | sed '1d')
for user in $usr
do
        echo $user
        grep $user running > usr
        rome=$(grep 'rome' usr | awk '{sum += $7} END {print sum}')
        jumbo=$(grep 'jumbo' usr | awk '{sum += $7} END {print sum}')
        total=$((128*rome+128*jumbo))
        echo -e "$(grep $user names | awk '{print $2}')'s usage is ${RED}$total${NC} cores"
done
echo
rome=$(grep 'rome' running | awk '{sum += $7} END {print sum}')
jumbo=$(grep 'jumbo' running | awk '{sum += $7} END {print sum}')
total=$((128*rome+128*jumbo))
echo "Group usage is $total cores"
echo
echo "Jobs pending $(cat pending | wc -l)"
for user in $usr
do
        echo "$(grep $user names | awk '{print $2}') has $(grep $user pending| wc -l) pending jobs"
done
echo "******************** END OF UPDATE ****************"