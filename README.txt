
Author: Yuanyuan Sun, Ling Yang
Date: 2015-05-27
Version: 1.0

This program implements our minCounter algorithm . We reference the basic code (more information introduced below) and have some modifications based on the basic code. 

basic code: http://sourceforge.net/projects/cuckoo-cpp/files/?source=navbar
   We obtain the basic code at the above link and download cuckoo-1.0.zip.
   Now we introduce the details of cuckoo-1.0.zip below:
   cuckoo-1.0/cuckoo.hpp : this file implements the basic algorithm of cuckoo hashing. This algorithm leverages the suitable selection to only select the position of kicking out elements. We always kick out the elements from the first hash table.
   cuckoo-1.0/cuckoo.pdf : the introduction to cuckoo.hpp. It illustrates details of cuckoo.hpp.
   Programming Language :C++
   Registered date :2011-05-30
   Author :Mariya Fomkina
   
   More information you can get from http://sourceforge.net/projects/cuckoo-cpp/?source=navbar and cuckoo-1.0/cuckoo.pdf.
   

minCounter code :
   MinCounter code implements this algorithm  based on the above basic code. 
   We consider the cuckoo.hpp in the basic code. Based on cuckoo.hpp, we mainly add some key parameters needed to minCounter and significantly modify the functions of add_new() , insert_i() and insert(). Besides the above modifications, we also add some necessary modifications at other functions to implement minCounter.
   The introduction of each .cpp or .h within minCounter below:
   cuckoo.hpp : this mainly implements the minCounter algorithm. More details have been mentioned at the above introduction.
   cuckoo-main.cpp : this includes the main function to examine minCounter. Moreover, it also includes three hash functions and some other functions.
   
   More information about minCounter algorithm you can get from the paper "MinCounter: An Efficient Cuckoo Hashing Scheme for Cloud Storage Systems" on MSST 2015.
   
This section introduces the ways to run the minCounter and the input format of parameters :
    compile: g++ *.cpp -o cuckoo
	run: ./cuckoo max_loop init_rate file_name
    input parameters: 
    	A) max_loop: when inserting a new element, if the times of kicking-out is equal to max_loop, we think there is an infinite loop.
    	B) init_rate: when initiating the hash table, we set the size of hash table be equal to the size of dataset  multiplied by init_rate.
	C) file_name: dataset file.
   
