The QuAcK software can be downloaded on GitHub as a Git repository
```
git clone https://github.com/pfloos/QuAcK.git
```

Then, one must define the variable `QUACK_ROOT`. For example, 
```
export QUACK_ROOT=$HOME/Work/QuAcK
```
You must also install [PySCF](https://pyscf.org) (for example using `pip`)
```
pip install pyscf
```

PySCF is used for the computation of one- and two-electron integrals (mainly).
