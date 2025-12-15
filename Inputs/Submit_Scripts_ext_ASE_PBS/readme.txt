PBS submit scripts for an infinit run with the external ASE engine

main_submit.sh is submitted by the user, and then automatically submits as many worker jobs as specified in the script
these worker jobs are automatically killed once the infinit task finishes

it is possible that worker jobs survive and idle if the main job dies before reaching the kill command (walltime, HPC crash, ...)

