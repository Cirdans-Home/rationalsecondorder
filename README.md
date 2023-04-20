# Efficient computation of the sinc matrix function for the integration of second-order differential equations

Rational Krylov Methods for the Integration of Second Order Differential Equations

## External codes

The code contained in this repository makes use of some [Chebfun functions](https://www.chebfun.org/). To add to your MATLAB 
environment simply run the following command from the MATLAB command window:
```
unzip('https://github.com/chebfun/chebfun/archive/master.zip')
movefile('chebfun-master', 'chebfun'), addpath(fullfile(cd,'chebfun')), savepath
```
The other fundamental ingredient is the [Rational Krylov toolbox](http://guettel.com/rktoolbox/), if you do not have it already
you can automatically download and install it, by simply copying and pasting the following two lines to your MATLAB command window:
```
unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
cd('rktoolbox'); addpath(fullfile(cd)); savepath
```
Among the methodologies available for calculating exponential sums, the [expmv code](https://github.com/higham/expmv) is used. 
This is added as a Git submodule to the repository. By cloning the repository the code is not added automatically, this can be 
done by going to the appropriate directory and doing `git pull`, or from the root directory with
```
git submodule init
git submodule update
```
Information about this code can be found at:
- A. H. Al-Mohy and N. J. Higham, "[Computing the action of the matrix exponential, with an application to exponential integrators](https://doi.org/10.1137/100788860)" SIAM J. Sci. Comput., 33(2):488--511, 2011.

## Contributors
- Lidia Aceto
- Fabio Durastante

## Cite as

The associated paper is available on arXiv, can be cited as:
```
@misc{aceto2023efficient,
      title={Efficient computation of the sinc matrix function for the integration of 
        second-order differential equations}, 
      author={Lidia Aceto and Fabio Durastante},
      year={2023},
      eprint={2304.09676},
      archivePrefix={arXiv},
      primaryClass={math.NA}
}
```

If you decide to use the code in a scientific publication, **also cite the related works** for the codes described in the previous section.
