#include "cuckoo.hpp"
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
using namespace std;

typedef pair<uint64_t,int> Data;//define the style of data

/**
 *hash  function
 */
unsigned long hash0(string str){
	unsigned long hash = 5381;
	for(int i=0;i<str.length();i++){
		char key_i=str.at(i);
		hash = ((hash << 5) + hash) + atoi( &key_i);  
	}
	return hash; 
}

unsigned long hash1(string str){
	unsigned long BitsInUnsignedInt = (unsigned long)(4 * 8);  
	unsigned long ThreeQuarters     = (unsigned long)((BitsInUnsignedInt  * 3) / 4);  
	unsigned long OneEighth         = (unsigned long)(BitsInUnsignedInt / 8);  
	unsigned long HighBits          = (unsigned long)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);  
	unsigned long hash              = 0;  
	unsigned long test              = 0;  

	for(int i = 0; i < str.length(); i++)  
	{  
		char key_i=str.at(i);
		hash = (hash << OneEighth) + atoi( &key_i);  
		if((test = hash & HighBits)  != 0)  
		{  
			hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));  
		}  
	}  
	return hash; 
}

unsigned long hash2(string str){
	 
	unsigned long seed = 131; // 31 131 1313 13131 131313 etc..  
      unsigned long hash = 0;  
      for(int i = 0; i < str.length(); i++)  
      {  
		 char key_i=str.at(i);
         hash = (hash * seed) +  atoi( &key_i);  
      }  
      return hash; 
}

template<class Key>
class Hasher{
public:
	size_t operator()(Key key, size_t hash_num,int table_length)const{
		unsigned long hash=0;
		string skey;
		 stringstream ss;
		 ss<<key;
		 ss>>skey;


		switch(hash_num) {
		case 0:
			hash=hash0(skey);
			hash=hash % table_length;
			break;
		case 1:
			hash = 0;  
			hash = hash1(skey);
			hash=hash % table_length;
			break;
		case 2:
			hash=0;
			hash = hash2(skey);
			do{
			    hash=hash % table_length;
			}while(hash>=table_length);
			break;
		}
		return hash;
	}
}; 

/**
 *equal  function
 */
template<class Key>
class Equal{
	public:
		bool operator()(Key key ,Key k){
			if(key==k)
				uint64_t t=k;
			return key==k?true:false;
		}
};


 /**
  * Generate a random number.
  */
  size_t genRandom(size_t rangeStart, size_t rangeEnd) {
      size_t r;
      do{
            r = rangeStart + (size_t)((rangeEnd + 1.0) * rand()/(RAND_MAX + 1.0));
      } while (r < rangeStart || r > rangeEnd);
      return r;
  }
  
 /**
  * transform hex to 10 system.
  */
  uint64_t changeNUM(string str,int len){
	  uint64_t sum=0;int temp=0;
	  for(int i=0;i<len;i++){
		  if(str[i]>='0'&&str[i]<='9')
			  temp=str[i]-48;
		  else if(str[i]>='a'&&str[i]<='f')
			  temp=str[i]-'a'+10;
		  else if(str[i]>='A'&&str[i]<='F')
			  temp=str[i]-'A'+10;
		  else
			  temp=0;
		  sum=sum*16+temp;
	  }
	  return sum;
  }

  
  /**
  * main function.
  */
int main(int argc, char** argv)
{
    if (argc != 4) {
        printf("\n[Usage]: %s max_loop init_rate dataset\n\n", argv[0]);
        printf("\tmax_loop\tthe iteration times bound\n"
                "\tinit_rate\tset the multiple based on dataset size\n"
                "\tdataset\t\tdataset filename\n\n");
        exit(0);
    }
	printf("argv[0]: %d \n",argv[0]);
	printf("argv[1]: %d \n",argv[1]);
	printf("argv[2]: %d \n",argv[2]);
	
	stringstream sl,sr,sf;
	size_t max_loop;double init_rate;
	sl<<argv[1];
	sl>>max_loop;// max_loop = argv[1];
	sr<<argv[2];
	sr>>init_rate;// init_rate = argv[2];
	printf("max_loop: %d \n",max_loop);
	printf("init_rate: %f \n",init_rate);
	
	
	
	size_t num=0;
	string name;
	sf<<argv[3];
	sf>>name;
	cout<<"name:"<<name<<endl;
	string frname=name+".data.txt";//get read file name


	FILE *fread;
	if(!(fread=fopen(frname.c_str(),"r+"))){
		printf("Failure to open test_dataset.txt！\n");
		return 0;
	}
	char s[32]; 
	if(fgets(s,32,fread)==NULL)//read the size of dataset 
		printf("reading the size of dataset is worry!!");
	stringstream ss;
	ss<<s;
	ss>>num;
	if(fgets(s,14,fread)==NULL)//read black 
		printf("reading the second row is worry!!");
	printf("reading the second row !!\n");


	
	string ir,ml;
	stringstream ssl,ssr;
	ssl<<init_rate;
	ssl>>ir;
	ssr<<max_loop;
	ssr>>ml;
	string dirname="minCounter//";
	dirname+="randnum_rate_";
	dirname+=ir;
	dirname+="//max_loop_";
	dirname+=ml;
	dirname+="//";
	dirname+=name;
	string fname=dirname+".minCounter_result.txt";// get write file name 
	FILE *fp;
	if(!(fp=fopen(fname.c_str(),"w+"))){
		printf("Failure to open test_result.text ！\n");
	}


	fprintf(fp,"minCounter : \n\n");
	fprintf(fp,"the initial rate of table : %f\n", INIT_TABLE);
	fprintf(fp,"the max loop of table : %d\n\n", MAX_LOOP);

	clock_t start,finish;//start time
    start=clock();

	cuckoo<uint64_t, int, Hasher<uint64_t>, Equal<uint64_t> > Cuckoo(num,max_loop,init_rate);
	//printf("constructe the cuckoo!!\n");

	fprintf(fp,"The total number of dataset :%d \n",num);
	fprintf(fp,"The total size of hash tables: %d \n",Cuckoo.length());
	fprintf(fp,"The size of each hash table :%d \n\n",Cuckoo.len_part());
	
	/************set test data************/
	Data k;
	int insert_res=-1;
	int insert_num=0;
	

	int i=0;
	while(fgets(s,32,fread)!=NULL){
		stringstream st;uint64_t key;
		st<<s;
		st>>key;
		
		k=make_pair(key,0);
		insert_res=Cuckoo.insert(k);//insert k to cuckoo hash
		if(insert_res!=-1)
			insert_num++;
		if(insert_res==-2){
			fprintf(fp,"when first insertion operation failed:\n");
	        fprintf(fp,"    The number of items within hash tables:%d \n",Cuckoo.first_kick_size());
	        double fkick_rate=(double)Cuckoo.first_kick_size()/(double)Cuckoo.length();
	        fprintf(fp,"    The utilization rate of hash tables::%lf \n",fkick_rate);
			Cuckoo.printTotal_counter(fp);
	        Cuckoo.getStat_counter(fp);
			cout<<"the first failure to insert key 's number of kick-out："<<Cuckoo.getTotal_counter()<<endl;
		}
		i++;
	}

	fclose(fread);

	finish=clock();

	fprintf(fp,"the excution time of randwalk : %f second \n",difftime(finish,start)/1000.00);  
	
	fprintf(fp,"When insertion operation completed:\n");
	fprintf(fp,"The number of items within hash tables:%d \n",Cuckoo.size());
	fprintf(fp,"The utilization rate of hash tables:%lf \n \n",Cuckoo.load_factor());
	fprintf(fp,"The number of inserting failed items:%d \n",Cuckoo.fail_insert_times());
	double fail_rate=(double)Cuckoo.fail_insert_times()/(double)num;
	fprintf(fp,"the success rate of insertion:%lf \n \n",1.0-fail_rate);

	Cuckoo.printTotal_counter(fp);
	Cuckoo.getStat_counter(fp);

	fclose(fp);

	/******TEST!!!!*********/
	cout<<"MAX_LOOP="<<MAX_LOOP<<endl;
	cout<<"the number of Key："<<num<<endl;
	cout<<"the number of inserted Key："<<insert_num<<endl;
	cout<<"When insertion operation completed,the utilization rate of hash tables:"<<Cuckoo.load_factor()<<endl;
	cout<<"when first insertion operation failed:"<<endl;
	cout<<"The number of items within hash tables:"<<Cuckoo.first_kick_size()<<endl;
	double fkick_rate=(double)Cuckoo.first_kick_size()/(double)Cuckoo.length();
	cout<<"The utilization rate of hash tables:"<<fkick_rate<<endl;
	cout<<"the total number of kick-out："<<Cuckoo.getTotal_counter()<<endl;
	
	return 0;
}

