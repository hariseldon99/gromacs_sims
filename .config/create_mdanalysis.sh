#!/bin/bash
conda create -n mdanalysis-develop numpy cython
conda activate mdanalysis-develop
git clone --depth=50 https://github.com/MDAnalysis/mdanalysis.git
cd mdanalysis
pip install package/
