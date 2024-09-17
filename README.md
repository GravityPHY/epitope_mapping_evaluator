# Epitope Mapping Evaluator

## Install  

```
pip install -e .
mamba create -n "mapping" python=3.10.0
mamba install pytorch torchvision torchaudio cudatoolkit=11.8 -c pytorch
mamba install -c dglteam/label/th23_cu118 dgl
```

### apbs
Download the source code from https://github.com/Electrostatics/apbs and unzip. 
```angular2html
export LD_LIBRARY_PATH=$HOME/apbs/lib:${LD_LIBRARY_PATH} 
export PATH=$HOME/apbs/bin:${PATH}
```

### pdb2pqr
