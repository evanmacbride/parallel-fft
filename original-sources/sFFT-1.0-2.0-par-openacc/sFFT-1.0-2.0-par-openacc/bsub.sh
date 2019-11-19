bsub -P ctp -U ctptest -R "select[ostypemajor==RHEL6]" -We 00:20 -x -n 1 -R "select[priolow]" -R "select[hardwarenickname==cobalt]" -eo err -oo out ./my_run.sh
