

Require to add external fork to git
to run infretisrun from branch infpp from fork https://github.com/dz24/infretis/

maybe this works
git pull

and checkout branch

git checkout infpp

To see if it works you can run the example

bash runner.sh

and from inftools, run

inft plot_ens -pp

# To analyze results using pyretisanalyse -i retis.rst

create fakeretis.rst file.

# convert infretis_data.txt to 0**/pathensemble.txt files
python3 conv_inf_py.py

then 

pyretisanalyse -i fakeretis.rst


# For creating initial paths for any other system

use init_paths_pp.py
which require a initial trajectory and its assoicated order parameter values in a text file.
