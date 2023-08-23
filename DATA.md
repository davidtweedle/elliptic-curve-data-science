### Obtaining the data

---------------------

The files `ec_curvedata_X.gp' were obtained from the LMFDB.
We processed the files using the script 'write_aps()' found in [generate-data.gp](gpscript/generate-data.gp) as follows

| input filename | output filename | conductor range | bound for p |
|----------------|-----------------|-----------------|-------------|
| ec_curvedata_7500_10000.gp | first_example.data | 7500-10000 | $p\leq 100000 $|
| ec_curvedata_50000.gp | cond_50k_i.data for i=1,2,...,14 | 1-50000 | $p\leq 50000$ |
| ec_curvedata_200000.gp (first 5672 curves)| cond_200k_1.data |  100002-100800 | $p\leq 200000$ |

I generated some data from scratch for the following curves.
This was done using the functions construct_box_curves(), construct_twists_sqrt_n1(), construct_twists_sqrt_n3(), construct_twists_sqrt_N() from the file generate-data.gp.
If you want to do the same, you can try to use gp2c-run -g generate-data.gp and run the scripts from gp.

| curve family | CM |filename | parameter range | bound for p |
|--------------|----|-----|-----------------|-------------|
| $y^2 = x^3 + Ax+b$ | |  second_example.data | $(A,B)\in [-64,64]^2$| $p\leq 100000$ |
| $y^2 = x^3 + Dx$ | $\mathbb{Z}[i]$ | sqrt_n1_D_10000_p_100000.data | $-10000\leq D\leq 10000$| $p\leq 100000$ |
| $y^2 = x^3 + D$ | $\mathbb{Z}[(1+\sqrt{-3})/2]$ | sqrt_n3_D_10000_p_100000.data | $-10000\leq D\leq 10000$| $p\leq 100000$ |
| $y^2 = x^3-270D^2x+1512D^3$  | $\mathbb{Z}[\sqrt{-2}]$ |sqrt_n2_D_10000-p_100000.data | $-10000\leq D\leq 10000$| $p\leq 100000$ |
| $y^2 = x^3 - 35D^2x -98 D^3$ |$\mathbb{Z}[(1+\sqrt{-7})/2]$ | sqrt_n7_D_10000_p_100000.data | $-10000\leq D\leq 10000$ | $p\leq 100000$|

