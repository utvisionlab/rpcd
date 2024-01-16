## ABOUT ##

The RPCD toolbox contains functions for low-multilinear rank approximation of high-order high-dimensional tensors.

Please see the end of this README for a list of the algorithms available.

This toolbox is copyright (C) 2021 by Mohammad Hamed & Reshad Hosseini and is distributed under the terms of the GNU General Public License (GPL) version 3 (or later).

Contact: [Reshad Hosseini](mailto:reshad.hosseini@ut.ac.ir) or  [Mohammad Hamed](mailto:mohammad.hamed@uantwerp.be)


## Quick installation guide ##

* Unzip and copy the whole rpcd directory you just downloaded in a location of your choice on disk, say, in /my/directory/.

* Go to /my/directory/rpcd/ at the Matlab command prompt and execute 'addPaths'.

## Directory structure ##

<pre>

./ The top directory with README.md
 |data                      - data directory
   |--- airquality/       
   |--- brainq/         
   |--- coil-100/
   |--- hsi/
   |--- yale/
   load_coil.m              - modified version of "load_boats" provided by DTucker authors to load coil-100 dataset
 |src/                      - source codes
   |--- rpcd.m              - pure RPCD method
   |--- rpcd_plus.m         - RPCD method with precision update
 |test/                     - test codes
   |--- test_real.m         - a script for test on real data
   |--- test_synthesis.m    - a script for test on synthesis data

</pre>

## Data files ##

To test with the real data, download each dataset from the following link and put it in the designated directory.

* Yale database:
	http://www.cad.zju.edu.cn/home/dengcai/Data/Yale/Yale_64x64.mat

* Brainq, Boats, Air Quality, HSI:
	https://datalab.snu.ac.kr/dtucker/

* Coil-100:
	https://www.kaggle.com/jessicali9530/coil100


## HOW TO CITE ##

If you find RPCD useful in your work, please cite the following two papers:

<pre>

@article{hamed2021riemannian,
  title={Riemannian preconditioned coordinate descent for low multi-linear rank approximation},
  author={Hamed, Mohammad and Hosseini, Reshad},
  journal={arXiv preprint arXiv:2109.01632},
  year={2021}
}

</pre>
