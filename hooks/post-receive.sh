#!/bin/sh
dest=/home/cayek/Projects/Thesis/lfmm

echo "post-receive hook"

while read oldrev newrev ref
do
    if [[ $ref =~ .*/master$ ]];
    then
        echo "Master ref received.  Deploying master branch to production..."
        git --work-tree=$dest --git-dir=/home/cayek/Gits/2017/lfmm.git checkout -f
	      cd $dest
	      make lfmm_install
        source activate MaThese
	      make lfmm_install
    else
        echo "Ref $ref successfully received.  Doing nothing: only the master branch may be deployed on this server."
    fi
done
