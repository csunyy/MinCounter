#ifndef _CUCKOO_HPP_
#define _CUCKOO_HPP_

#include <stdio.h>
#include <cstring>
#include <assert.h>
#include <utility>
#include <stdlib.h>
#include <stack>
#include <iostream>
#include <fstream>
using namespace std;


const static size_t D = 3;
const static size_t MAX_LOOP = 50;
const static double STEP = 1.1;
const static double INIT_TABLE = 1.1;



/**
 * @param Key key type in hash.
 * @param Value value type in hash.
 * @param Hash function that gets data of type Key and some number
 * (which is the number of appropriate hash function) and returns size_t
 * (e.g. Hash(int, size_t)).
 * @param Equal predicator that compares two Key values.
 * @example cuckoo<int, int, Hasher, std::key_equal_to<int> > Cuckoo;
 */
template <class Key, class Value, class Hash, class Equal>
class cuckoo {
private:
  typedef pair<Key, Value> Data;
  /**
   * @variable The number of hash functions (thus arrays also)
   * that will be used in the program (can be >= 2).
   */
  size_t d_;
  /**
   * @variable The initial length of the whole structure.
   * When you know the approximate number of records to be used,
   * it is a good idea to take this value in 1.05-1.1 times more and
   * small value of step.
   */
  size_t init_length_;
  /**
   * @variable The maximum number of kick cycles during
   * insertion before rehash.
   */
  size_t max_loop_;
  /**
   * @variable The ratio of increasing the size of hash during rehash.
   * The less it is the less memory will be used but the more time is needed.
   */
  double step_;
  
  double init_rate_;
  /**
   * @variable The hash function object (template parameter by default).
   */
  Hash hasher_;
  /**
   * @variable The equal predicator object (template parameter by default).
   */
  Equal key_equal_;
  /**
   * @variable The array of vectors, each of which is hash array.
   */
  Data** data_;
  /**
   * @variable The array of vectors, each of which is the number of kick-out of each bucket .
   */
  int** kick_counter;
  /**
   * @variable The array of flags indicating existence of the element in hash.
   */
  char* exists_;
  /**
   * @variable The total length of all the hash arrays.
   */
  size_t len_;
  /**
   * @variable The length of every hash array.
   */
  size_t len_part_;
  /**
   * @variable The actual number of elements in cuckoo hash.
   */
  size_t size_;
  /**
   * @variable The flag that anounces that rehash was made recently.
   */
  bool is_rehashed_;

  /**
  *@variable the stack that store the elements which fail to insert in hash table 
  */
  stack<Data> fail_stack_;

  /**
  *@variable whether it is the first time occuring kick out
  */
  bool is_first_kick_;    

  /**
  *@variable the number of elements inserted in hash table when it occurs kick out at first time
  */
  size_t fkick_size_;  
  

public:
  friend class c_iterator;
  friend class const_iterator;

  /**
   * Used with const cuckoo objects.
   */
  class const_iterator {
  private:
    size_t pos;
    const cuckoo* hash;
    const_iterator(const size_t p, const cuckoo* h) : pos(p), hash(h) {}

  public:
    friend class cuckoo;

    /**
     * Default constructor.
     */
    const_iterator() : pos(0), hash(NULL) {}

    void operator=(const const_iterator& it) {
      pos = it.pos;
      hash = it.hash;
    }

    const_iterator(const const_iterator& it) {
      *this = it;
    }

    const_iterator& operator++() {
      assert(hash != NULL);
      if (pos != hash->len_) ++pos;
      while ((pos < hash->len_) && (!(hash->get_exists(pos)))) {
        ++pos;
      }
      return *this;
    }

    const_iterator operator++(int) {
      const_iterator res = *this;
      this->operator++();
      return res;
    }

    const Data& operator*() const {
      return (*hash).data_from(pos);
    }

    const Data* operator->() const {
      return &((*hash).data_from(pos));
    }

    bool operator==(const const_iterator &it) {
      return (pos == it.pos) && (hash == it.hash);
    }

    bool operator!=(const const_iterator &it) {
      return !(*this == it);
    }
  };

  /**
   * Used with non-const cuckoo objects.
   */
  class c_iterator {
  private:
    size_t pos;
    cuckoo* hash;
    c_iterator(size_t p, cuckoo* h) : pos(p), hash(h) {}

  public:
    friend class cuckoo;

    /**
     * Default constructor.
     */
    c_iterator() : pos(0), hash(NULL) {}

    /**
     * Convertion to const_iterator.
     */
    operator const_iterator() {
      return const_iterator(pos, hash);
    }

    void operator=(const c_iterator &it) {
      pos = it.pos;
      hash = it.hash;
    }

    c_iterator(const c_iterator &it) {
      *this = it;
    }

    c_iterator& operator++() {
      assert(hash != NULL);
      if (pos != hash->len_) ++pos;
      while ((pos < hash->len_) && !(hash->get_exists(pos))) {
        ++pos;
      }
      return *this;
    }

    c_iterator operator++(int) {
      c_iterator res = *this;
      this->operator++();
      return res;
    }

    Data& operator*() {
      return (*hash).data_from(pos);
    }

    Data* operator->() {
      return &((*hash).data_from(pos));
    }

    bool operator==(const c_iterator &it) {
      return (pos == it.pos) && (hash == it.hash);
    }

    bool operator!=(const c_iterator &it) {
      return !(*this == it);
    }
  };

private:

  /**
   * Check whether there is a Data element at pos.
   *
   * @param pos Position in hash arrays.
   * @return true if some element presents on this position, false otherwise.
   */
  bool get_exists(size_t pos) const {
    return (bool)(exists_[pos >> 3] & (1 << (pos - ((pos >> 3) << 3))));
  }

  /**
   * Set flag of existence of Data element at pos.
   *
   * @param pos Position in hash arrays.
   */
  void set_exists(size_t pos) {
    exists_[pos >> 3] = exists_[pos >> 3] | (1 << (pos - ((pos >> 3) << 3)));
  }

  /**
   * Unset flag of existence Data element at pos.
   *
   * @param pos Position in hash arrays.
   */
  void unset_exists(size_t pos) {
    exists_[pos >> 3] = exists_[pos >> 3] & (~(1 << (pos - ((pos >> 3) << 3))));
  }

  /**
   * Initialize all the variables of cuckoo.
   */
  void init() {
    len_part_ = init_length_ / d_ + 1;
    len_part_ = ((len_part_ + 7) >> 3) << 3;
    len_ = len_part_ * d_;
    data_ = (Data**)(malloc(d_ * sizeof(Data*)));
	kick_counter=(int**)(malloc(d_ * sizeof(int*)));
    for (size_t i = 0; i < d_; ++i) {
      data_[i] = (Data*)(malloc(len_part_ * sizeof(Data)));
	  memset(data_[i],0,len_part_ * sizeof(Data));
	  kick_counter[i] = (int*)(malloc(len_part_ * sizeof(int)));
	  memset(kick_counter[i],0,len_part_ * sizeof(int));
    }
    exists_ = (char*)(malloc((len_ >> 3) * sizeof(char)));
    for (size_t i = 0; i < len_ >> 3; ++i) exists_[i] = 0;
    size_ = 0;
    is_rehashed_ = false;
	is_first_kick_=false;
	fkick_size_=0;
  }

  /**
   * Copy information from Cuckoo to clean object.
   *
   * @param Object for information to be copied from.
   */
  void copy(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    d_ = Cuckoo.d_;
    init_length_ = Cuckoo.init_length_;
    max_loop_ = Cuckoo.max_loop_;
    step_ = Cuckoo.step_;
    hasher_ = Cuckoo.hasher_;
    key_equal_ = Cuckoo.key_equal_;
    init();
    const_iterator it = Cuckoo.begin();
    size_t final_len = Cuckoo.length();
    while (it.pos != final_len) {
      insert(*it);
      ++it;
    }
  }

  /**
   * Clear all the data from cuckoo.
   */
  void clear_all() {
    for (size_t i = 0; i < d_; ++i) {
      free(data_[i]);
	  free(kick_counter[i]);
    }
    free(data_);
	free(kick_counter);
    free(exists_);
  }

  /**
   * Get Data element from pos.
   *
   * @param pos Position in hash arrays.
   * @return Reference to Data element.
   */
  inline Data& data_from(size_t pos) const {
    return data_[pos / len_part_][pos % len_part_];
  }

  /**
   * Check whether element at pos has key k.
   *
   * @param k Key value.
   * @param pos Position in hash arrays.
   * @return true if element at pos is equal to k.
   */
  inline bool is_here(const Key& k, size_t pos) const {
    return get_exists(pos) && key_equal_(data_from(pos).first, k);
  }

  /**
   * Return hash function result for key k.
   *
   * @param k Key value.
   * @param hash_num The number of hash function from Hash family.
   * @return Hash value.
   */
  inline size_t hash(const Key& k, size_t hash_num) const {
    return hasher_(k, hash_num,len_part_) /*% len_part_*/;
  }

  /**
   * Increase size of exists_ up to len_temp_.
   *
   * @param len_temp_ New size of exists_.
   */
  void update_exists(size_t len_temp_) {
    char* t = (char*)(malloc((len_ >> 3) * sizeof(char)));
    char* s = t;
    char* f = s + (len_ >> 3);
    for (; s < f; ++s) {
      *s = 0;
    }
    for (size_t i = 0; i < d_; ++i) {
      memcpy(t + (i * (len_part_ >> 3)), exists_ + (i * (len_temp_ >> 3)), (len_temp_ >> 3));
    }
    std::swap(t, exists_);
    free(t);
  }

  /**
   * Increase size of data_ up to len_temp_.
   *
   * @param len_temp_ New size of data_.
   */
  void update_data(size_t len_temp_) {
    for (size_t i = 0; i < d_; ++i) {
      Data* t = (Data*)(malloc(len_part_ * sizeof(Data)));
      memcpy(t, data_[i], len_temp_*sizeof(Data));
      std::swap(t, data_[i]);
      free(t);
    }
  }

  /**
   * Increase size of kick_counter up to len_temp_.
   *
   * @param len_temp_ New size of kick_counter.
   */
  void update_kick_counter(size_t len_temp_){
	 for (size_t i = 0; i < d_; ++i) {
      int* t = (int*)(malloc(len_part_ * sizeof(int)));
	  memset(t,0,len_part_ * sizeof(int));
      memcpy(t, kick_counter[i], len_temp_*sizeof(int));
      std::swap(t, kick_counter[i]);
      free(t);
	 }
  }

  /**
   * Rehash all the cuckoo (i.e. change size and replace Data elements).
   */
  void rehash() {
    size_t len_temp_ = len_part_;

    len_part_ = (size_t)(len_part_ * step_);
    len_part_ = ((len_part_ + 7) >> 3) << 3;
    len_ = len_part_ * d_;

    update_exists(len_temp_);
    update_data(len_temp_);
	update_kick_counter(len_temp_);

    size_t n = 0;
    c_iterator it = begin();
    while (it != end()) {
      size_t i = it.pos / len_part_;
      size_t j = it.pos % len_part_;
      ++n;
      if (j != hash((*it).first, i)) {
        Data t = *it;
        remove(it);
        add_new(t);
        if (is_rehashed_) {
          it = begin();
          is_rehashed_ = false;
        }
      } else {
        ++it;
      }
    }
    is_rehashed_ = true;
  }

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
   * Add new Data element.  
   * Used the minCounter algorithm to add_new new element into minCounter cuckoo hash
   *
   * @param p New element.
   * @return 0 to finish a recursion if it was needed.
   */
  size_t add_new(Data p) {

	size_t temp=-2;
    for (size_t i = 0; i< max_loop_; ++i)
    {
		if(i==0){
		    temp = d_ +1; //record the random selected position at the last step 
		}
		for(size_t j = 0;j < d_; j++)
		{
			if(temp == j ) 
				j +=1;  //if j is the random position at the last step, skipping this position
			if(j>=d_) 
				break;
			size_t pos = hash(p.first,j);
			bool exists = get_exists(j * len_part_ +pos);
			if(!exists)
			{
				data_[j][pos] = p;
				set_exists(j * len_part_ +pos);
				++ size_;
				return 0;
			}
		}


		int lmax_counter_pos=0;
		for(size_t j = 1;j < d_; j++){
			if(temp == j )  //if j is the random position at the last step, skipping this position
				j +=1; 
			if(j>=d_) 
				break;
			if(lmax_counter_pos==temp&&lmax_counter_pos==0){//if temp=0£¬only compare the counter between j=1 and j=2
				lmax_counter_pos++;
				j++;
			}
			size_t pos = hash(p.first,j);
			size_t pos_m = hash(p.first,lmax_counter_pos);
			if(kick_counter[lmax_counter_pos][pos_m] > kick_counter[j][pos])
				lmax_counter_pos=j;

		}
		size_t r=lmax_counter_pos;
		temp = r;
        size_t pos = hash(p.first, r);
        std::swap(p,data_[r][pos]);

 		kick_counter[r][pos]++;

    }
	if(is_first_kick_==false){
		is_first_kick_=true;
		fkick_size_=size_;
		fail_stack_.push(p);
		return -2;
	}

	fail_stack_.push(p);
	return -1;
  }

  /**
   * Remove element from cuckoo.
   *
   * @param it c_iterator to element to be removed.
   * @return c_iterator to next element after removed.
   */
  c_iterator remove(c_iterator& it) {
    unset_exists(it.pos);
    --size_;
	kick_counter[it.pos / len_part_][it.pos % len_part_]=0;//clean the times of kick_out
    return ++it;
  }

public:
  /**
   * Default constructor.
   *
   * @param d The number of hash functions (thus arrays also)
   * that will be used in the program (can be >= 2).
   * @param init_length The initial length of the whole structure.
   * When you know the approximate number of records to be used,
   * it is a good idea to take this value in 1.05-1.1 times more and
   * small value of step.
   * @param max_loop The maximum number of kick cycles during
   * insertion before rehash.
   * @param step The ratio of increasing the size of hash during rehash.
   * The less it is the less memory will be used but the more time is needed.
   * @param hasher The hash function object (template parameter by default).
   * @param equal The equal predicator object (template parameter by default).
   */
  /*explicit cuckoo(size_t d = D, size_t init_length = INIT_LENGTH,
                  size_t max_loop = MAX_LOOP, double step = STEP,
                  const Hash& hasher = Hash(),
                  const Equal& equal = Equal())
    : d_(d), init_length_(init_length),
      max_loop_(max_loop), step_(step),
      hasher_(hasher), key_equal_(equal) {
    init();
  }*/

  explicit cuckoo(size_t init_length, size_t max_loop ,double init_rate,
                  size_t d = D,double step = STEP,
                  const Hash& hasher = Hash(),
                  const Equal& equal = Equal())
    :d_(d),max_loop_(max_loop), step_(step),
	  hasher_(hasher), key_equal_(equal){

		  init_length_=init_length*init_rate;
		  init();
  }
  
  explicit cuckoo(size_t init_length,size_t d = D,
                  size_t max_loop = MAX_LOOP, double step = STEP,
                  const Hash& hasher = Hash(),
                  const Equal& equal = Equal())
    :d_(d),max_loop_(max_loop), step_(step),
	  hasher_(hasher), key_equal_(equal){

		  init_length_=init_length*INIT_TABLE;
		  init();
  }

  /**
   * Destructor.
   */
  ~cuckoo() {
    clear_all();
  }

  /**
   * Copy information from cuckoo of the same type.
   *
   * @param Cuckoo The source of data
   */
  cuckoo<Key, Value, Hash, Equal>& operator=
  (const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    clear_all();
    copy(Cuckoo);
    return *this;
  }

  /**
   * Copy constructor.
   *
   * @param Cuckoo The source of information
   */
  cuckoo(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    copy(Cuckoo);
  }

  /**
   * Constructor from range [fisrt, last).
   *
   * @param first The begin of range c_iterator.
   * @param last The end of range c_iterator.
   */
  /*template <class InputIterator>
  cuckoo(InputIterator first, InputIterator last,
         const Hash& hasher = Hash(),
         const Equal& equal = Equal()) {
    d_ = D;
    init_length_ = INIT_LENGTH;
    max_loop_ = MAX_LOOP;
    step_ = STEP;
    hasher_ = hasher;
    key_equal_ = equal;
    init();
    size_t last_pos = last.pos;
    for (c_iterator it = first; it.pos != last_pos; ++it) {
      add_new(*it);
    }
  }*/

  /**
   * Update parameter of cuckoo, deleting all the data from it.
   *
   * @warning For test only!!!
   */
  /*void set_up(size_t d = D, size_t init_length = INIT_LENGTH,
              size_t max_loop = MAX_LOOP, double step = STEP) {
    clear_all();
    d_ = d;
    init_length_ = init_length;
    max_loop_ = max_loop;
    step_ = step;
    init();
  }*/

  /**
   * Operator==
   *
   * @param Cuckoo Object to be compared with.
   * @return true if objects are equal (they are references to the same object).
   */
  bool operator==(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    return (*this).data_ == Cuckoo.data_;
  }

  /**
   * @see operator==
   */
  bool operator!=(const cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    return !(*this == Cuckoo);
  }

  /**
   * Swap data with Cuckoo.
   *
   * @param Cuckoo Cuckoo of the same type
   */
  void swap(cuckoo<Key, Value, Hash, Equal>& Cuckoo) {
    std::swap(d_, Cuckoo.d_);
    std::swap(init_length_, Cuckoo.init_length_);
    std::swap(max_loop_, Cuckoo.max_loop_);
    std::swap(step_, Cuckoo.step_);
    std::swap(hasher_, Cuckoo.hasher_);
    std::swap(key_equal_, Cuckoo.key_equal_);
    std::swap(data_, Cuckoo.data_);
	std::swap(kick_counter, Cuckoo.kick_counter);//swap kick_counter
    std::swap(exists_, Cuckoo.exists_);
    std::swap(len_, Cuckoo.len_);
    std::swap(len_part_, Cuckoo.len_part_);
    std::swap(size_, Cuckoo.size_);
  }

  /**
   * Get c_iterator to begin of cuckoo.
   *
   * @return c_iterator to the first element of cuckoo.
   */
  inline c_iterator begin() {
    c_iterator it = c_iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }

  /**
   * Get const_iterator to begin of cuckoo (for const objects).
   *
   * @see begin()
   */
  inline const_iterator begin() const {
    const_iterator it = const_iterator(0, this);
    if (!get_exists(it.pos)) ++it;
    return it;
  }

  /**
   * Get c_iterator to end of cuckoo.
   *
   * @return c_iterator to the position AFTER last element of cuckoo.
   */
  inline c_iterator end() {
    return c_iterator(len_, this);
  }

  /**
   * Get const_iterator to begin of cuckoo (for const objects).
   *
   * @see end()
   */
  inline const_iterator end() const {
    return const_iterator(len_, this);
  }

  /**
   * Get value by key or created pair key-value.
   *
   * @param k Key value.
   * @return Reference to Value, associated with k.
   */
  Value& operator[](const Key& k) {
    return (*((this->insert(make_pair(k, Value()))).first)).second;
  }

  /**
   * Erase data at c_iterator.
   *
   * @param it c_iterator to Data element to be removed.
   */
  void erase(c_iterator it) {
    remove(it);
  }

  /**
   * Erase range of Data elements.
   *
   * @param first The begin of range c_iterator.
   * @param last The end of range c_iterator.
   */
  void erase(c_iterator first, c_iterator last) {
    while (first.pos != last.pos) {
      first = remove(first);
    }
  }

  /**
   * Erase element by key.
   *
   * @param k Key value.
   * @return 1 if element was erased and 0 if it didn't exist.
   */
  size_t erase(const Key& k) {
    c_iterator it = find(k);
    if (it.pos != len_) {
      remove(it);
      return 1;
    }
    return 0;
  }

  /**
   * Find element by key.
   *
   * @param k Key value
   * @return c_iterator to element or to end of cuckoo, if element
   * doesn't exist.
   */
  c_iterator find(const Key& k) {
    size_t dist = 0;
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      size_t position = pos + dist;
      if ((get_exists(position)) && key_equal_(data_[i][pos].first, k)) {
        return c_iterator(position, this);
      }
      dist += len_part_;
    }
    return end();
  }

  /**
   * Find element by key (for const objects).
   *
   * @param k Key value
   * @return Const_iterator to element or to end of cuckoo, if element
   * doesn't exist.
   */
  const_iterator find(const Key& k) const {
    size_t dist = 0;
    for (size_t i = 0; i < d_; ++i) {
      size_t pos = hash(k, i);
      size_t position = pos + dist;
      if (key_equal_(data_[i][pos].first, k) && (get_exists(position))) {
        return const_iterator(position, this);
      }
      dist += len_part_;
    }
    return end();
  }

  /**
   * Count number of elements with this key.
   *
   * @param k Key value.
   * @return 1 if element exists and 0 otherwise.
   */
  size_t count(const Key& k) const {
    return (find(k)).pos != len_;
  }

  /**
   * Find range of elements with key.
   *
   * @param k Key value.
   * @return Pair that determines the range [fisrt, last) or
   * pair with both iterators pointing to the end of cuckoo.
   */
  pair<c_iterator, c_iterator> equal_range(const Key& k) {
    c_iterator l = find(k);
    c_iterator r = l;
    return std::make_pair<c_iterator, c_iterator>(l, ++r);
  }

  /**
   * Find range of elements with key (for const objects).
   *
   * @param k Key value
   * @return Pair that determines the range [fisrt, last) or
   * pair with both c_iterator pointing to the end of cuckoo.
   */
  pair<const_iterator, const_iterator> equal_range(const Key& k) const {
    const_iterator l = find(k);
    const_iterator r = l;
    return std::make_pair<const_iterator, const_iterator>(l, ++r);
  }

  /**
   * Insert Data element to cuckoo.
   *
   * @param k The new Data element.
   * @return Pair with c_iterator to existing element and bool value,
   * which is true if element was inserted or false if it existed before.
   */
  
   pair<c_iterator, int> insert_i(const Data& k) {
    c_iterator res = find(k.first);
    if (res.pos != len_) {
      return make_pair(res, -1);
    }
    int is_success=add_new(k);
	if(is_success==-1)
		return make_pair(c_iterator(hash(k.first, 0), this), 0);
    return make_pair(c_iterator(hash(k.first, 0), this), 1);
  }

   int insert(const Data& k) {
        c_iterator res = find(k.first);
        if (res.pos != len_) {
            return  -1;
        }
        int is_success=add_new(k);
		if(is_success==-2)
			return -2;
	    if(is_success==-1)
		   return  0;
        return  1;
  }

  /**
   * Insert range of Data elements.
   *
   * @param first The begin of range c_iterator.
   * @param last The end of range c_iterator.
   */
  template <class InputIterator>
  void insert(InputIterator first, InputIterator last) {
    size_t last_pos = last.pos;
    for (c_iterator it = first; it.pos != last_pos; ++it) {
      insert(*it);
    }
  }

  /**
   * Clear all data from cuckoo.
   */
  void clear() {
    clear_all();
    init();
  }

  /**
   * Check whether cuckoo is empty.
   *
   * @return true of cuckoo is empty and false otherwise.
   */
  inline bool empty() const {
    return (size_ == 0);
  }

  /**
   * Show size of cuckoo.
   *
   * @return Size of cuckoo.
   */
  inline size_t size() const {
    return size_;
  }

  /**
   * Show length of cuckoo (actual number of elements).
   *
   * @return Length of cuckoo.
   */
  inline size_t length() const {
    return len_;
  }

  inline size_t len_part() const{
	  return len_part_;
  }
  /**
  *Show the load of cuckoo.
  *
  *@ return load factor of cuckoo.
  */
  inline double load_factor() const{
	  double load_factor_= (double)size_ / (double)len_;
	  return load_factor_;
  }

  
  /**
  *Show the number of elements failing to insert into cuckoo hash.
  *
  *@ return the size of fail_stack_ .
  */
  inline size_t fail_insert_times() const{
	  return this->fail_stack_.size();
  }

  
  /**
  *Show the number of elements within hash tables when first insertion operation failed.
  *
  *@ return the number of elements within hash tables .
  */
  inline size_t first_kick_size() const{
	  return fkick_size_;
  }
  

/********only for test***********/
  public:
  
	 void printData(FILE *fp){
		  int format_=0;
		  for(int i=0;i<this->d_;i++){
			  fprintf(fp,"the Data of table %d: \n",i);
			  format_=0;
			  for(int j=0;j<this->len_part_;j++){
				  format_++;
				  fprintf(fp,"%ld",data_[i][j].first);
				  if(format_==20){
					  fprintf(fp,"\n");
					  format_=0;
				  }
			  }
			  fprintf(fp,"\n");
		  }
	  }
	  void printKick_counter(FILE *fp){
		  int format_=0;
		  for(int i=0;i<this->d_;i++){
			  fprintf(fp,"The values of counters of buckets in table %d: \n",i);
			  format_=0;
			  for(int j=0;j<this->len_part_;j++){
				  format_++;
				  fprintf(fp,"%d ",kick_counter[i][j]);
				  if(format_==50){
					  fprintf(fp,"\n");
					  format_=0;
				  }
			  }
			  fprintf(fp,"\n");
		  }
	  }

	  void printTotal_counter(FILE *fp){
		  int total_counter_=0;
		  int counter[3]={0,0,0};
		  for(int i=0;i<this->d_;i++){
			  for(int j=0;j<this->len_part_;j++){
				  counter[i]=counter[i]+kick_counter[i][j];
			  }
			  fprintf(fp,"The kick-out times of table %d : %d \n",i,counter[i]);
		  }
		  total_counter_=counter[0]+counter[1]+counter[2];
		  fprintf(fp,"The total kick-out times of hash tables: %d \n\n",total_counter_);
	  }
	  
	  size_t getTotal_counter(){
		  int total_counter_=0;
		  int counter[3]={0,0,0};
		  for(int i=0;i<this->d_;i++){
			  for(int j=0;j<this->len_part_;j++){
				  counter[i]=counter[i]+kick_counter[i][j];
			  }
		  }
		  total_counter_=counter[0]+counter[1]+counter[2];
		  return total_counter_;
	  }

	  /******count the range of counter*******/
	  void getStat_counter(FILE *fp){
		  int count[8];
		  memset(count,0,8*sizeof(int));
		  for(int i=0;i<this->d_;i++){
			  for(int j=0;j<this->len_part_;j++){
				  if(kick_counter[i][j]>=0&&kick_counter[i][j]<2)
					  count[0]++;
				  else if(kick_counter[i][j]>=2&&kick_counter[i][j]<4)
					  count[1]++;
				  else if(kick_counter[i][j]>=4&&kick_counter[i][j]<8)
					  count[2]++;
				  else if(kick_counter[i][j]>=8&&kick_counter[i][j]<16)
					  count[3]++;
				  else if(kick_counter[i][j]>=16&&kick_counter[i][j]<32)
					  count[4]++;
				  else if(kick_counter[i][j]>=32&&kick_counter[i][j]<64)
					  count[5]++;
				  else if(kick_counter[i][j]>=64&&kick_counter[i][j]<128)
					  count[6]++;
				  else if(kick_counter[i][j]>=128)
					  count[7]++;
			  }
		  }
		  fprintf(fp,"The number of values between 0~2: %d \n",count[0]);
		  fprintf(fp,"The number of values between 2~4: %d \n",count[1]);
		  fprintf(fp,"The number of values between 4~8: %d \n",count[2]);
		  fprintf(fp,"The number of values between 8~16: %d \n",count[3]);
		  fprintf(fp,"The number of values between 16~32: %d \n",count[4]);
		  fprintf(fp,"The number of values between 32~64: %d \n",count[5]);
		  fprintf(fp,"The number of values between 64~128: %d \n",count[6]);
		  fprintf(fp,"The number of values larger than 128: %d \n\n",count[7]);
	  }
};

#endif /* _CUCKOO_HPP_ */
