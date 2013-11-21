PipelineMiRNA
=============

Installation
------------

You may either check out the code from version control, or use the archives provided.

### Installing Perl dependencies

If you are using the archives:
```
perl Build.PL
./Build installdeps
```

If you are using cpanminus:
```
cat requirements.txt | sudo cpanm
```

If you are using Aptitude:
```
apt-get -s install $(cat apt-requirements.txt)
```
(Note it does not install ODF::lpOD − this one must be installed by hand).


### Installing programs

- The Shell script `install-programs.sh` should install pretty much all dependencies.
```
sh install-programs.sh programs/
```
- `RNAstemloop` has to be installed manually


### Installing data

The Pipeline relies on a few data files, most of which are under version control.

The remaining data is to be deployed using The Shell script `install-data.sh`
```
sh install-data.sh data/
```

Installation troubleshooting
----------------------------

- Make sure you have a `results` directory created
- Make sure the owner is www-data


