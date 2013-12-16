for i in `seq -w 1 384`
do
  cp juqueen_template.job juqueen.${i}.job
  sed -i "s/BLA/${i}/g" juqueen.${i}.job
  llsubmit juqueen.${i}.job
done
