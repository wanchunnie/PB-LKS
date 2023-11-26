# PB-LKS

This is the official code for our paper "PB-LKS: a python package for predicting Phage-Bacteria interaction through Local K-mer Strategy"
```
 +--------------------------------------------+
 | |---\  |----\       |       |  /   /----\  |
 | |    | |     |      |       | /   |        |
 | |___/  |____/   --  |       |/\    \----\  |
 | |      |     \      |       |  \         | |
 | |      |_____/      |_____  |   \  \____/  |
 +--------------------------------------------+
```



## Local environment setup
1. install conda(to manage environment)
2.  Change directory to the path of this project
      ```bash
      cd {your_path_to_PBLKS}
      ```
3. Run the following codes in your terminal
   ```shell
   conda create -n PB-LKS python=3.10
   conda deactivate (if base environment is activated)
   conda activate PB-LKS
   pip install -r requirements.txt
   ``` 

## How to use
Both a command line tool and a python package is provided in this project. Developers and researchers may select according to needs after setting local environment.

### Command line interface
A command line interface is implemented and can be used by running following commands in your terminal. 

   ```bash
   cd {your_path_to_PBLKS}
   # display help message of PB-LKS CLI
   python PB-LKS.py -h
   # run and display result with example input
   python PB-LKS.py -e
   ```
other detailed arguments and their usages are listed as follows
```Shell
usage: python PB-LKS.py [-h] [-p PHAGE] [-b BAC] [-e] [-xgb] [-fea] [-ba] [-o OUTPUT] [-d]

options:
  -h, --help            show this help message and exit
  -p PHAGE, --phage PHAGE
                        path/folder to your phage sequence file(in fasta format)
  -b BAC, --bac BAC     path/folder to your bacteria sequence file(in fasta format)
  -e, --example         run model with example input and quit
  -xgb, --xgboost       run prediction with xgboost(default model is based on RandomForest)
  -fea, --feature       show 10 features of most importance and exit
  -ba, --batch          run prediction in batch, -p,-b should be folder if set True
  -o OUTPUT, --output OUTPUT
                        the folder you want to save batch prediction results (results will be printed into terminal by default)

  -d, --detail          if set true, print detailed prediction results(single prediction), including important features and decision path.
```

### python package
To make it more convenient for researchers and developers to use our model in more personalized tasks, a python package is also implemented. 


example code:
```python
import PBLKS
bacteria_path = 'example1.fasta'
phage_path = 'example2.fasta'
# directly gives the result of prediction
result, prob =  PBLKS.predict(bacteria_path, phage_path)
# gets descriptors from files
# can be used to reconstruct or train other type of models  
feas = PBLKS.get_descriptor(bacteria_path, phage_path)

# predicts the interaction of each bacteria-phage pair in both folder
bac_folder = "your_folder_path"
phage_folder  = "your_folder_path"
out_put_dir = "your_dest_folder"
# output'll be directed to dest folder, creating 2 result files
PBLKS.predict_in_batch(bac_folder, phage_folder, output_dir=out_put_dir)
# printing results into terminal
PBLKS.predict_in_batch(bac_folder, phage_folder, output_dir=None)

# get decision path
decision_path = PBLKS.show_decision_path(bacteria_path, phage_path)
# get top 10 important features
# result is orgnized in dict like: {kmer: importance}
features = PBLKS.show_feature_importance(top_cnt=10)
```

## decision path visualization
Several interfaces to visualize how our model make decisions are implemented. 

Except calling   `-d` or `PBLKS.show_decision_path`, we implemented a script to visualize decision trees in our model.

Users can run following command in termianl
```bash
python make_tree_plots.py
```
This'll create a detailed visualization of every decision tree into `/PBLKS/Trees`.

## other models
An xgboost based model is also included in `/models`.

can be used by adding `-xgb` argument in CLI or adding parameter `use_xgboost=True` when calling methods in our package 


## Source code

Visit this git repository to get source code of PB-LKS:

	https://github.com/wanchunnie/PB-LKS
