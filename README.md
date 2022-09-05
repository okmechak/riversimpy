# RIVERSIM

Simulation of river growth using model based on Laplace equation[1].
Mathematicaly, in this program we solve PDE equation using Finite Element Method(FEM). And as result program produces VTK file, which contains solution and its details, which furthere can be viewed in [__ParaView__](https://www.paraview.org/) 

For build process look at riversim repo riversimpy branch.

How to update package on PyPi([source](https://medium.com/@joel.barmettler/how-to-upload-your-python-package-to-pypi-65edc5fe9c56)):
```bash
python3 setup.py sdist
twine upload dist/*
pip3 install riversim --upgrade
```

## Master thesis

LaTeX generation setup:

https://nevalsar.hashnode.dev/writing-latex-documents-with-ubuntu-and-visual-studio-code

```bash
sudo apt install texlive-science texlive-extra-utils latexmk
```

and install LaTeX Workshop extension for VS Code.