#RIVERSIM

Simulation of river growth using model based on Laplace equation[1].
Mathematicaly, in this program we solve PDE equation using Finite Element Method(FEM). And as result program produces VTK file, which contains solution and its details, which furthere can be viewed in [__ParaView__](https://www.paraview.org/) 

For build process look at riversim repo riversimpy branch.

How to update package on PyPi([source](https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56)):
```bash
python setup.py sdist
twine upload dist/*
pip install riversimpy --upgrade
```
