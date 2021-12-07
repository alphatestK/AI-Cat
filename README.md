# AI-Cat
AI-Cat test version 0.1

This software package implements AI-Cat method that takes an arbitary reaction to predict possible mechanism.

The package provides two major functions:

- Train a AI-cat model with a customized dataset.
- Predict mechanism of possible reactions with a pre-trained AI-Cat model.

The code is written for Python 2.6 or 2.7. The current version is only a development version for functional test and the more clean and user firendlly version will be updated later.

The following paper describes the details of the AI-Cat framework:

Selectivity Origin of Catalytic Glycerol Hydrogenolysis on Cu Resolved from Artificial Intelligence based Monte Carlo Tree Search

## Table of Contents

- [Prerequisites](#prerequisites)
- [Usage](#usage)
  - [Gnerate a customized dataset](#generate-a-customized-dataset)
  - [Train a AI-Cat model](#train-a-ai-cat-model)
  - [Predict reactions with a pre-trained AI-Cat model](#predict-reactions-with-a-pre-trained-ai-cat-model)
- [Authors](#authors)

## Prerequisites

This package requires:
- [Openbabel](https://openbabel.org/wiki/Main_Page)
- [Tensorflow](https://www.tensorflow.org/install)
- [Numpy](https://numpy.org/)

## Usage

To use the program, you first need to add the ```$YOURWORKDIR/lib``` to your ```/PYTHONPATH```.

### Generate a customized dataset 

To input reaction dataset to AI-Cat, you will need to define a customized dataset. Note that this is required for both training and predicting. 


You can create a customized dataset with the following files:

1. `allgoodsect.arc`: a .arc file for all reactions collected. Each reaction is recorded by 3 structures(IS/TS/FS).

2. `allname.arc`: a .arc file for minimum structures of different intermdiates. The GM sturcture of intermediates can be provided by this file.

3. `sampled.arc`: a .arc file for labelling the high value/trustworthy structures (if needed). 

4. `Arc2RPdata.py`: the main program for processing data.

5. `reactiontype.py`: the program for generating additional data (e.g. estimated ZPE info).

The structure of the `workdir` should be:

```
workdir
├── Arc2RPdata.py
├── allgoodsect.arc
├── allname.arc
├── sampled.arc
├── reactiontype.py
├── ...
```
There is an example of raw data in the repository: `data/rawdata/`.

You can generate the dataset by running:
```bash
python Arc2RPdata.py
```

After the `testdata` and `allpair.arc` are generated, you can generate the additonal database by running:
```bash
python reactiontype.py
```

There is an example of customized data in the repository: `data/customized data/`. 


### Train a AI-Cat model

Before training a new AI-Cat model, you will need to:

- [Gnerate a customized dataset](#generate-a-customized-dataset) to store the possible reaction patterns and kinetic info.

Then, in directory `train`, you can train a AI-Cat model for your customized dataset by:

```bash
ln -s ../data/customized data/data* .
python train_net1.py
python train_net2.py
```

You can set the number of training, validation, and test data by easilly changing the code.(future version will change it more convenient way.)
After training, you will get two files in `train` directory.

- `policymodel.h5`: stores the Reaction Pattern model.
- `infonet.h5`: stores the Kinetic info model.

### Predict reactions with a pre-trained AI-Cat model

Before predicting the material properties, you will need to:

- [Gnerate a customized dataset](#generate-a-customized-dataset) to store the possible reaction patterns, and link the generated `*.pkl` file.
- Obtain two [pre-trained AI-Cat models](#pre-trained model) named `policymodel.h5` and `infonet.h5`.

Then, in directory `pre-trained model` with all these files, you can set the searching target in `searchjob` file and predict the mechanism of possible reaction by running:

```bash
ln -s ../data/customized data/*.pkl .
python MCTS_ZPE.py
```

For instace, you can predict the mechanism of glycerol hydrogenolysis to 1,2-propanediol in `pre-trained/` by setting the mode as **1** and inputfile as `glycerolto12.arc`, in which 2 structures (the reactants and possible product) are provided.

And if you have no clear idea of the target product, you can do not set the products by giving an empty structure **2** or set the mulit targets in mode **2** or set vauge target in mode **3**.

After predicting, you will get the 3 files in directory:

- `result`: stores the predicted possible reaction mechanism of reactant to target products.
- `goodresult`: stores the predicted possible products generated from reactants in serach.
- `allstep` : stores the searching log.

## Authors

This software was primarily written by Peilin Kang who was advised by [Prof. Zhipan Liu](https://zpliu.fudan.edu.cn/). 

