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

### Installing programs

- The Shell script `install-deps.sh` should install pretty much all dependencies.
```
sh install-deps.sh programs/
```
- `RNAstemloop` has to be installed manually

Installation troubleshooting
----------------------------

- Make sure you have a `results` directory created
- Make sure the owner is www-data


