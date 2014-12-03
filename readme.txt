Source code accompanying the following paper:

@inproceedings{Bauer:2014:FMF:2645710.2645735,
 author = {Bauer, Josef and Nanopoulos, Alexandros},
 title = {A Framework for Matrix Factorization Based on General Distributions},
 booktitle = {Proceedings of the 8th ACM Conference on Recommender Systems},
 series = {RecSys '14},
 year = {2014},
 isbn = {978-1-4503-2668-1},
 location = {Foster City, Silicon Valley, California, USA},
 pages = {249--256},
 numpages = {8},
 url = {http://doi.acm.org/10.1145/2645710.2645735},
 doi = {10.1145/2645710.2645735},
 acmid = {2645735},
 publisher = {ACM},
 address = {New York, NY, USA},
 keywords = {automatic differentiation, machine learning, matrix factorization, maximum likelihood},
} 

Installation and usage instructions:

Running the algorithm requires R, Rtools, Python and SymPy. For a Windows platform, these can be downloaded here:

http://cran.r-project.org/bin/windows/base/
http://cran.r-project.org/bin/windows/Rtools/
https://www.python.org/downloads/
https://github.com/sympy/sympy/releases

In order to increase speed, just-in-time compilation as well as C++ is used for the computationally most expensive part of the code. 
For this purpose, the Rcpp and the inline package have to be installed. This can be done by the following command in RGui:

install.packages("Rcpp")
install.packages("inline")

Likewise, the following libraries have to be installed:

install.packages("rSymPy")
install.packages("VGAM")
install.packages("SuppDists")

After this has been done, the program can be started by opening and running dmfmain.r or pasting the code in R (RGui.exe). 
The distributions discussed in the paper are provided as templates which can be adjusted. 
The default directory where the files of this archive are stored is C:\DMF and this can be changed at the beginning of the dmfmain file.
By default, the data set is stored in the data subfolder, in the format (user-id, item-id, count/rating), whereby user- and item-ids are integer. 
Data sets can be added by placing the corresponding files in the data folder and adding the read-in details (like header and separation) in the createuseritemdata.r file.
In the generatetrainingvalidationdata.r file one can specify for each data set how many data points are used for training. 
These will be randomly selected and the remaining ones are used for validation. The number of batches to use can be specified in the same file.

In the dmfmain file, one can specify at the beginning the distribution, which data set should be used, the matrix of (hyper-)parameters, the number of epochs and basic filtering (minimum and maximum count values as well as the minimum number of observations for each user and item) as it is described in the paper.
The parameter matrix contains the values for the learning rate, the regularization parameter, the momentum parameter, the number of features, the initial value parameter and additional parameters depending on the distribution. 
In the template additional values can be added by a comma. In the program grid search will be performed for all possible combinations of these parameter values.

During training in each epoch the current mean absolute error on the validation set is printed. 
For the current parameters these values can be saved by the R command save(validationerror, file = path) for further analysis and plotting, e.g., by R's plot function.

No warranty is given. Please contact josef.bauer@ku.de in case of questions.
