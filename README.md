# Nepre_coordinate
A program can extract Neighborhood relationship coordinate data from pdb file
## Summary
This project include the following things:
(1) coordinate _AA.py is a program can extract Neighborhood relationship between Amino Acid.
(2) coordinate_chain.py is a program can extract Neighborhood relationship between chain and chain.
(3) testdataAA is the test data for program coordinate_AA.py
(4) testdatachain is the test data for program coordinate_chain.py
(5) The result about Neighborhood relationship  saved as a dict in npy file type,you can use program Readnpy.py to view the result.

## Usage

### Requirement
The program is implemented with **Python 2.7** and **Python 3.6**.

The following package is required: **numpy**
```
pip install numpy
```

### Running the program
Download the code from GitHub, then you are ready to go.

```
python coordinate_AA.py filename the_path_of_ file output_path
```

* For example, to predict the 2CBA.pdb and if the SSBONDPredict is saved in '/Users/ssb/Desktop', the commandline will be:

```
python coordinate_AA.py "2CBA.pdb" "/Users/ssb/Desktop" "/Users/ssb/Desktop"
```

* You can find the result named 2CBA.npy saved in the diectory "/User/ssb/Desktop".

* And If you want to view the result, you can use the program Readnpy,for example:

```
python Readnpy.py "/Users/ssb/Desktop/2CAB.npy"
```

* The usage about coordinate_chain.py is the same with coordinate_AA.py

## Copyright
SSBONDPredict is created by liulab of Beijing Compulational Science Research Center.

