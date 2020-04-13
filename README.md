# PointSite: a point cloud segmentation tool for identification of protein ligand binding atoms
This repository contains the author's implementation for the paper:<br>
**PointSite: a point cloud segmentation tool for identification of protein ligand binding atoms** [[bioRxiv]](https://www.biorxiv.org/content/early/2019/11/05/831131.full.pdf)<br>
Created by [Zhen Li](https://github.com/icemansina),  [Xu Yan](https://github.com/yanx27) and [Sheng Wang](https://github.com/realbigws)
![](https://raw.githubusercontent.com/PointSite/PointSite_Inference/master/example/pipline.png)


## Setup

Tested with CUDA 9.0, Ubuntu 18.04, Python 3.6 with [Conda](https://www.anaconda.com/) and PyTorch 1.1.

```
git clone --recursive https://github.com/PointSite/PointSite_Inference.git
cd PointSite_Inference/
./install.sh
```

WARNING: To install the package successfully, users shall use the latest version Anaconda, such as [Anaconda2-2019.10](https://repo.anaconda.com/archive/Anaconda2-2019.10-Linux-x86_64.sh).


## Inference
 ```
python inference.py 
--gpu: GPU index, if you have not GPU, just ignore it
--output: output root (required)
--data: data root, only support .xyz file (required)
--select_list: TXT file for selected protein name, default None
--num_vote: voting number in inference (default 25, larger number can archieve more stable and high performance)
```

## Running Example
```
conda activate pointsite_inference
python inference.py --output blind_out --data example/blind --select_list example/blind_list
conda deactivate
```

Note that the above input data (in '.XYZ' format) contain the ground-truth label of binding atoms. Run below script for identifying binding atoms on unlabeled data in '.PDB' files.
```
chmod +x ./pointsite_run.sh
./pointsite_run.sh example/blind_list example/blind blind_out `pwd`
```


## Visualization
You will get .obj file in output folder, please use [MeshLab](http://www.meshlab.net/) to visualize.
![](https://raw.githubusercontent.com/PointSite/PointSite_Inference/master/example/result.png)


## Training and Testing Data
Users may find the training data [here](https://drive.google.com/drive/folders/1iqe1741ohxh4b3p6wk3toDgS-IVTPXtq?usp=sharing); <br>
Users may find the test data [here](https://drive.google.com/drive/folders/1-y5JNjPgZvxdiY-U0w1oTQBJn5cexC02?usp=sharing); <br>
Users may find the evaluation results [here](https://drive.google.com/drive/folders/1siev2t-DIAtGRk8u1K7TSnTUtsw1Btnd?usp=sharing). 


## Reference
3D Semantic Segmentation with Submanifold Sparse Convolutional Networks, CVPR 2018 
[facebookresearch/SparseConvNet](https://github.com/facebookresearch/SparseConvNet/tree/master/)
