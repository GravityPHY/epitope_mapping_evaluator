# epitope_mapping_evaluator

## Install  

```
pip install -e .
mamba create -n "mapping" python=3.10.0
mamba install pytorch torchvision torchaudio cudatoolkit=11.8 -c pytorch
mamba install -c dglteam/label/th23_cu118 dgl
```