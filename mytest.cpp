/*
  File: mytest.cpp
  Author: Andrew Greenberg
  Date: 4/18/25
  Desc: file containing and which runs all tests for the DnaDb class
 */

#include "dnadb.h"

// copy Random class (by Prof. Donyaee) to be able to generate random numbers
// also include necessary libraries for it
#include <math.h>
#include <algorithm>
#include <random>
#include <vector>
enum RANDOM {UNIFORMINT, UNIFORMREAL, NORMAL, SHUFFLE};

class Random {
public:
  Random(){}
  Random(int min, int max, RANDOM type=UNIFORMINT, int mean=50, int stdev=20) : m_min(min), m_max(max), m_type(type)
  {
    if (type == NORMAL){
      //the case of NORMAL to generate integer numbers with normal distribution
      m_generator = std::mt19937(m_device());
      //the data set will have the mean of 50 (default) and standard deviation of 20 (default)
      //the mean and standard deviation can change by passing new values to constructor
      m_normdist = std::normal_distribution<>(mean,stdev);
    }
    else if (type == UNIFORMINT) {
      //the case of UNIFORMINT to generate integer numbers
      // Using a fixed seed value generates always the same sequence
      // of pseudorandom numbers, e.g. reproducing scientific experiments
      // here it helps us with testing since the same sequence repeats
      m_generator = std::mt19937(10);// 10 is the fixed seed value
      m_unidist = std::uniform_int_distribution<>(min,max);
    }
    else if (type == UNIFORMREAL) { //the case of UNIFORMREAL to generate real numbers
      m_generator = std::mt19937(10);// 10 is the fixed seed value
      m_uniReal = std::uniform_real_distribution<double>((double)min,(double)max);
    }
    else { //the case of SHUFFLE to generate every number only once
      m_generator = std::mt19937(m_device());
    }
  }
  void setSeed(int seedNum){
    // we have set a default value for seed in constructor
    // we can change the seed by calling this function after constructor call
    // this gives us more randomness
    m_generator = std::mt19937(seedNum);
  }
  void init(int min, int max){
    m_min = min;
    m_max = max;
    m_type = UNIFORMINT;
    m_generator = std::mt19937(10);// 10 is the fixed seed value
    m_unidist = std::uniform_int_distribution<>(min,max);
  }
  void getShuffle(vector<int> & array){
    // this function provides a list of all values between min and max
    // in a random order, this function guarantees the uniqueness
    // of every value in the list
    // the user program creates the vector param and passes here
    // here we populate the vector using m_min and m_max
    for (int i = m_min; i<=m_max; i++){
      array.push_back(i);
    }
    shuffle(array.begin(),array.end(),m_generator);
  }

  void getShuffle(int array[]){
    // this function provides a list of all values between min and max
    // in a random order, this function guarantees the uniqueness
    // of every value in the list
    // the param array must be of the size (m_max-m_min+1)
    // the user program creates the array and pass it here
    vector<int> temp;
    for (int i = m_min; i<=m_max; i++){
      temp.push_back(i);
    }
    std::shuffle(temp.begin(), temp.end(), m_generator);
    vector<int>::iterator it;
    int i = 0;
    for (it=temp.begin(); it != temp.end(); it++){
      array[i] = *it;
      i++;
    }
  }

  int getRandNum(){
    // this function returns integer numbers
    // the object must have been initialized to generate integers
    int result = 0;
    if(m_type == NORMAL){
      //returns a random number in a set with normal distribution
      //we limit random numbers by the min and max values
      result = m_min - 1;
      while(result < m_min || result > m_max)
	result = m_normdist(m_generator);
    }
    else if (m_type == UNIFORMINT){
      //this will generate a random number between min and max values
      result = m_unidist(m_generator);
    }
    return result;
  }

  double getRealRandNum(){
    // this function returns real numbers
    // the object must have been initialized to generate real numbers
    double result = m_uniReal(m_generator);
    // a trick to return numbers only with two deciaml points
    // for example if result is 15.0378, function returns 15.03
    // to round up we can use ceil function instead of floor
    result = std::floor(result*100.0)/100.0;
    return result;
  }

  string getRandString(int size){
    // the parameter size specifies the length of string we ask for
    // to use ASCII char the number range in constructor must be set to 97 - 122
    // and the Random type must be UNIFORMINT (it is default in constructor)
    string output = "";
    for (int i=0;i<size;i++){
      output = output + (char)getRandNum();
    }
    return output;
  }

  int getMin(){return m_min;}
  int getMax(){return m_max;}
private:
  int m_min;
  int m_max;
  RANDOM m_type;
  std::random_device m_device;
  std::mt19937 m_generator;
  std::normal_distribution<> m_normdist;//normal distribution
  std::uniform_int_distribution<> m_unidist;//integer uniform distribution
  std::uniform_real_distribution<double> m_uniReal;//real uniform distribution
};

// hashCode is used to initialize databases, but own hash function is used in place of m_hash's
unsigned int hashCode(const string str) {
  unsigned int val = 0 ;
  const unsigned int thirtyThree = 33 ;  // magic number from textbook
  for (long unsigned int i = 0 ; i < str.length(); i++)
    val = val * thirtyThree + str[i] ;
  return val ;
}

// constant strings for the pass/fail statements in couts in tests
const string FAIL_STATEMENT = "*****TEST FAILED: ";
const string PASS_STATEMENT = "     TEST PASSED: ";

class Tester {
public:
  //Name: testDnaDbConstructorNormal
  //Case: normal
  //Expected result: after creating a DnaDb object, it will have the correct member variable
  // values (appropriate initializations and those of the passed values)
  bool testConstructorNormal() {
    // define constants to create the DnaDb object
    const int SIZE = 103; // prime number greater than MINPRIME
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;

    // define constants for expected member variable values of the DnaDb object
    const int DEFAULT_SIZE = 0;
    const int DEFAULT_NUM_DELETED = 0;          
    const int DEFAULT_TRANSFER_INDEX = 0;

    // create the DnaDb object
    DnaDb database(SIZE, HASH, PROBING);
      
    // check all member variables of the DnaDB object
    // fail the test if any variables have unexpected values
    if (database.m_hash != HASH or database.m_newPolicy != PROBING
	or database.m_currentTable == nullptr or database.m_currentCap != SIZE
	or database.m_currentSize != DEFAULT_SIZE
	or database.m_currNumDeleted != DEFAULT_NUM_DELETED
	or database.m_currProbing != PROBING or database.m_oldTable != nullptr
	or database.m_oldCap != DEFAULT_SIZE or database.m_oldSize != DEFAULT_SIZE
	or database.m_oldNumDeleted != DEFAULT_NUM_DELETED or database.m_oldProbing != DEFPOLCY
	or database.m_transferIndex != DEFAULT_TRANSFER_INDEX) {
      cout << FAIL_STATEMENT << "constructor failed the normal case\n";
      return false;
    }

    // check that the currentTable was initialized to have nullptrs in it
    for (int i = 0; i < DEFAULT_SIZE; i++) {
      if (database.m_currentTable[i] != nullptr) {
	cout << FAIL_STATEMENT << "constructor failed the normal case\n";
	return false;
      }
    }
    
    // pass if no failures with any tests
    cout << PASS_STATEMENT << "constructor passed the normal case\n";
    return true;
  }

  //Name: testHashFunctionNormal
  //Case: normal
  //Expected result: hashing a random sequence will result in a unique integer, which will,
  // when % cap, result in a fairly unique index. There should be minimal collisions
  bool testHashFunctionNormal() {
    const int CAPACITY = 103; // size of the array
    const int NUM_SEQUENCES = 206; // number of sequences to insert into the array
    int numSeq[CAPACITY] = {0}; // array to count how many sequences are in each slot
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    int randNum = 0; // used to generate random DNA sequences
    int value = 0; // used to store the hash function's return value, so later lines are cleaner
    int LEEWAY = 3; // ideally every bucket has size/capacity items in it. LEEWAY allows for
    // how much more items can be in a bucket and still pass (ie: size 10 and capacity 5
    // expects 2 things per bucket, but LEEWAY of 1 would also allow 3 things in buckets)
    
    // define constants to create the DnaDb object
    const int SIZE = 103; // prime number greater than MINPRIME
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;
    
    // dnadb object to be able to use hash()
    DnaDb dnadb(SIZE, HASH, PROBING);
    
    
    // random object to create random DNA sequences
    Random random(1, 4);

    // simulate adding all sequences to the array
    for (int j = 0; j < NUM_SEQUENCES; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	randNum = random.getRandNum();
	switch (randNum) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }

      // increment the array's slot for whichever index the current sequence would be mapped to
      value = dnadb.hash(currentSequence);
      numSeq[(value % (unsigned int)CAPACITY)] += 1;

      //reset current sequence
      currentSequence = "";
    }
    
    // ideally each slot in the array would have NUM_SEQUENCES / CAPACITY sequences in it.
    // if it is more than a couple above that, then fail. else, pass
    for (int i = 0; i < CAPACITY; i++) {
      if (numSeq[i] > (NUM_SEQUENCES / CAPACITY + LEEWAY)) {
	cout << FAIL_STATEMENT << "hashFunction failed the normal case\n";
	return false;
      }
    }
    
    cout << PASS_STATEMENT << "hashFunction passed the normal case\n";
    return true;
  }

  //Name: testInsertNormal
  //Case: normal
  //Expected result: after a decent number of insertions to a database (less than is required
  // to trigger a rehash), the database will contain those DNA and in an expected position
  bool testInsertNormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 5; // less than a fourth of the capacity, avoid rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;
        
    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations
    
    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // fail if the insert fails (all dna should be unique and with valid members)
      if (not dnadb.insert(tempDna)) {
	cout << FAIL_STATEMENT << "insert failed the normal case\n";
	return false;
      }
    }

    // check that all expected DNA objects are present in the database
    for (int i = 0; i < SIZE; i++) {
      // if any inserted DNA aren't in the database, fail
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "insert failed the normal case\n";
	return false;
      }

      // if any inserted DNA has an incorrect location, fail
      if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	cout << FAIL_STATEMENT << "insert failed the normal case\n";
	return false;
      }
    }

    // check that the size member variable updated correctly after the insertions
    if (dnadb.m_currentSize != SIZE) {
      cout << FAIL_STATEMENT << "insert failed the normal case\n";
      return false;
    }
    
    // failed no tests
    cout << PASS_STATEMENT << "insert passed the normal case\n";
    return true;
  }

  //Name: testGetDNANormalCollide
  //Case: normal
  //Expected result: after inserting a decent number of DNA objects into a database, without
  // triggering a rehash and with some colliding and some non-colliding keys/sequences, the
  // getDNA function will return the correct items
  bool testGetDNANormalCollide() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 103; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    DNA tempDna; // used to add dna to the database
    const int SIZE = 11; // less than a fourth of the capacity, avoid rehashing. Number of
                         // keys with only 3 unique indexes, will be probed
    bool unique = true; // represents if currend data in DNA is unique
    
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;
        
    // arrays to store DNA objects' data
    string dnaSeq[] = {"CTTCTATCAT", "TTTCCCTTGC", "CCGACGGAAT", "GCCAGACCTG", "CCTCCCAGTC",
		     "TAGGTCCGAC", "CTGGTTCCAA", "ACTCCTGATG", "CTTAAAGACT", "TCCTCGGTGC",
		     "CCTACAGGAT"}; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random object for generating DNA object
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new location is unique, else generate again (decrement j)
      unique = true;
      for (int h = 0; h < j; h++) {
	if (dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // search for all the added DNA
    // fail the test if any aren't found
    for (int i = 0; i < SIZE; i++) {    
      if (dnadb.getDNA(dnaSeq[i], dnaLoc[i]).m_sequence != dnaSeq[i]
	  or dnadb.getDNA(dnaSeq[i], dnaLoc[i]).m_location != dnaLoc[i]) {
	cout << FAIL_STATEMENT << "getDNA failed the normal (collide) case\n";
	return false;
      }
    }
      
    // all DNA were found, passed
    cout << PASS_STATEMENT << "getDNA passed the normal (collide) case\n";
    return true;
  }

  //Name: testGetDNANormal
  //Case: normal
  //Expected result: after a decent number of insertions (without triggering a rehash), getDNA
  // will be able to locate all inserted DNA
  bool testGetDNANormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    const int SIZE = CAPACITY / 5; // less than a fourth of the capacity, avoid rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = DOUBLEHASH;
    bool unique = true; // represents if currend data in DNA is unique
    
    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true;
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }
    

    // search for all the added DNA
    // fail the test if any aren't found
    for (int i = 0; i < SIZE; i++) {
      if (dnadb.getDNA(dnaSeq[i], dnaLoc[i]).m_sequence != dnaSeq[i]
	  or dnadb.getDNA(dnaSeq[i], dnaLoc[i]).m_location != dnaLoc[i]) {
	cout << FAIL_STATEMENT << "getDNA failed the normal case\n";
	return false;
      }
    }

    // all DNA were found, passed
    cout << PASS_STATEMENT << "getDNA passed the normal case\n";
    return true;
  }
  
  //Name: testGetDNAError
  //Case: error
  //Expected result: getDna will return an empty DNA object / not find the target DNA when
  // searching for a non-existant target DNA
  bool testGetDNAError() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    string fakeSequence = "AACTGACTATGCGAT"; // is greater than seqLength, so won't be in dnadb
    int fakeLoc = MINLOCID - 1; // invalid location, so won't be in dnadb
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    const int SIZE = CAPACITY / 5; // less than a fourth of the capacity, avoid rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;
    bool unique = true; // represents if currend data in DNA is unique
    
    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true;
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }


    // search for all the added DNA
    // pass the test if the fake/non-existant DNA isn't found
      if (dnadb.getDNA(fakeSequence, fakeLoc).m_sequence == ""
	  and dnadb.getDNA(fakeSequence, fakeLoc).m_location == 0) {
	cout << PASS_STATEMENT << "getDNA passed the error case\n";
	return true;
      }
	
    // if the fake DNA was found, fail
    cout << FAIL_STATEMENT << "getDNA failed the error case\n";
    return false;
  }

  //Name: testInsertNormalRehashPartial
  //Case: normal
  //Expected result: after inserting a decent number of DNA to a database, enough to trigger
  // a rehashe, the current and old table will contain all the data that was inserted
  bool testInsertNormalRehashPartial() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 101; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    const int SIZE = 52; // more than half of the capacity to trigger two rehashes
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;
    bool unique = true; // represents if currend data in DNA is unique

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true;
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // fail if the insert fails (all dna should be unique and with valid members)
      if (not dnadb.insert(tempDna)) {
	cout << FAIL_STATEMENT << "insert failed the normal (rehash, partial) case\n";
	return false;
      }
    }

    // check that all expected DNA objects are present in the database
    for (int i = 0; i < SIZE; i++) {
      // if any inserted DNA aren't in the database, fail
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "insert failed the normal (rehash, partial) case\n";
	return false;
      }

      // if any inserted DNA has an incorrect location, fail
      // item could be in the current or old table, but if the above test passed, then it
      // must exist in one of them - check which by checking its sequence
      if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])] != nullptr
	  and dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_sequence == dnaSeq[i]) {
	if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "insert failed the normal (rehash, partial) case\n";
	  cout << endl;
	  return false;
	}
      }
      else { // located in old table
	if (dnadb.m_oldTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "insert failed the normal (rehash, partial) case\n";
	  return false;
	}
      }
    }

    // check that the size member variable updated correctly after all the insertions
    // data is random and hash function isn't perfect, so can't guarantee the value of any
    // of the sizes of the tables, but the old table should (likely) still contain dna,
    // and the current table should have between 0 and the size
    if (dnadb.m_currentSize >= SIZE or dnadb.m_currentSize == 0 or dnadb.m_oldSize == 0) {
      cout << FAIL_STATEMENT << "insert failed the normal (rehash, partial) case\n";
      return false;
    }

    // failed no tests
    cout << PASS_STATEMENT << "insert passed the normal (rehash, partial) case\n";
    return true;
  }

  //Name: testInsertError
  //Case: error
  //Expected result: insert will fail after being given an invalid DNA
  bool testInsertError() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    DNA fakeDna("AAA", (MINLOCID - 1), false); // have invalid location
    const int SIZE = 0; // nothing will be added to the database
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // insert the fake DNA into the database
    // modify the temp DNA object, then add it to the database
    // it should fail

    // fail if the insert fails (all dna should be unique and with valid members)
    if (dnadb.insert(fakeDna)) {
      cout << FAIL_STATEMENT << "insert failed the error case\n";
      return false;
    }
    
    // check that nothing has changed within the database
    // check that the size member variable updated correctly after the insertions
    if (dnadb.m_currentSize != SIZE) {
      cout << FAIL_STATEMENT << "insert failed the error case\n";
      return false;
    }

    // failed no tests
    cout << PASS_STATEMENT << "insert passed the error case\n";
    return true;
  }
  
  //Name: testRemoveNormalUnique
  //Case: normal
  //Expected result: after inserting some DNA with non-colliding keys, removing some of those
  // keys will "lazy delete" them (set their m_used as false)
  bool testRemoveNormalUnique() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 101; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = 10;
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE] = {"GAACCGGAAC", "ATTTAACAGA", "ACGCATATGA", "CGCAACTTCT",
      "TTTGCCTCCA", "CCGAGGAAAA", "CTAAACTACC", "GCATCGGACT",
      "CCTCCCAGTC", "ACTATACCTT"}; // sequences, manually set to avoid
    // collisions
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new location is unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete half the data, which isn't enough to trigger a rehash
    for (int i = 0; i < SIZE; i += 2) {
      // set up tempDna for current target of deletion
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // remove the target, fail if anything goes wrong
      if (dnadb.remove(tempDna) == false) {
	cout << FAIL_STATEMENT << "remove failed the normal (unique) case\n";
	return false;
      }
    }

    // see if size is correct
    if (dnadb.m_currentSize != (SIZE / 2)) {
      cout << FAIL_STATEMENT << "remove failed the normal (unique) case\n";
      return false;
    }

    // check if num deleted is accurate
    if (dnadb.m_currNumDeleted != (SIZE / 2)) {
      cout << FAIL_STATEMENT << "remove failed the normal (unique) case\n";
      return false;
    }

    // confirm that deleted data can't be found
    for (int i = 0; i < SIZE; i += 2) {
      if (getDNAIndex(dnadb, dnaSeq[i]) != -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (unique) case\n";
	return false;
      }
    }

    // confirm that non-deleted data can be found
    for (int i = 1; i < SIZE; i += 2) {
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (unique) case\n";
	return false;
      }
    }

    // passed all the tests
    cout << PASS_STATEMENT << "remove passed the normal (unique) case\n";
    return true;
  }
                         
  //Name: testRemoveNormalCollide
  //Case: normal
  //Expected result: after inserting some DNA with colliding keys, removing some of those
  // keys will "lazy delete" them (set their m_used as false)
  bool testRemoveNormalCollide() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 101; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = 10;
    const hash_fn HASH = hashCode;
    const prob_t PROBING = DOUBLEHASH;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE] = {"CTTCTATCAT", "TTTCCCTTGC", "CCGACGGAAT", "GCCAGACCTG",
			   "CCTCCCAGTC", "TAGGTCCGAC", "CTGGTTCCAA", "ACTCCTGATG",
			   "CTTAAAGACT", "TCCTCGGTGC"}; // sequences, manually set to have
                                                        // mulitple collisions
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new location is unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete half the data, which isn't enough to trigger a rehash
    for (int i = 0; i < SIZE; i += 2) {
      // set up tempDna for current target of deletion
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // remove the target, fail if anything goes wrong
      if (dnadb.remove(tempDna) == false) {
	cout << FAIL_STATEMENT << "remove failed the normal (collide) case\n";
	return false;
      }
    }

    // see if size is correct
    if (dnadb.m_currentSize != (SIZE / 2)) {
      cout << FAIL_STATEMENT << "remove failed the normal (collide) case\n";
      return false;
    }

    // check if num deleted is accurate
    if (dnadb.m_currNumDeleted != (SIZE / 2)) {
      cout << FAIL_STATEMENT << "remove failed the normal (collide) case\n";
      return false;
    }

    // confirm that deleted data can't be found
    for (int i = 0; i < SIZE; i += 2) {
      if (getDNAIndex(dnadb, dnaSeq[i]) != -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (collide) case\n";
	return false;
      }
    }

    // confirm that non-deleted data can be found
    for (int i = 1; i < SIZE; i += 2) {
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (collide) case\n";
	return false;
      }
    }

    // passed all the tests
    cout << PASS_STATEMENT << "remove passed the normal (collide) case\n";
    return true;
  }

  //Name: testRemoveNormalRehash
  //Case: normal
  //Expected result: after inserting and then removing DNA from the database, enough to trigger
  // a rehash
  bool testRemoveNormalRehash() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 307; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    DNA tempDna; // used to add dna to the database
    int seqLength = 10; // number of chars in a sequence
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = 150;
    const int NUM_TO_DELETE = 124; // number of DNA to remove, is more than 80% so a 
                                   // rehash will be triggered
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object data
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 as there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";
      
      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new location is unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete a little more than 80% of the data, to trigger a rehash
    for (int i = 0; i < NUM_TO_DELETE; i++) {
      // set up tempDna for current target of deletion
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // remove the target, fail if anything goes wrong
      if (dnadb.remove(tempDna) == false) {
	cout << FAIL_STATEMENT << "remove failed the normal (rehash) case\n";
	return false;
      }
    }
    
    // see if size is correct
    if (dnadb.m_currentSize > (SIZE - NUM_TO_DELETE)) {
      cout << FAIL_STATEMENT << "remove failed the normal (rehash) case\n";
      return false;
    }

    // confirm that deleted data can't be found
    for (int i = 0; i < NUM_TO_DELETE; i++) {
      if (getDNAIndex(dnadb, dnaSeq[i]) != -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (rehash) case\n";
	return false;
      }
    }

    // confirm that non-deleted data can be found
    for (int i = NUM_TO_DELETE; i < SIZE; i++) {
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "remove failed the normal (rehash) case\n";
	return false;
      }
    }

    // passed all the tests
    cout << PASS_STATEMENT << "remove passed the normal (rehash) case\n";
    return true;
  }

  //Name: testRemoveError
  //Case: error
  //Expected result: attempting to remove a non-existant DNA will fail
  bool testRemoveError() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 307; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    DNA tempDna; // used to add dna to the database
    int seqLength = 10; // number of chars in a sequence
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = 150;
    const int NUM_TO_DELETE = 124; // number of DNA to remove, is more than 80% so a
                                   // rehash will be triggered
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object data
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 as there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations
    
    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new location is unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all but one of the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < (SIZE - 1); i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete a little more than 80% of the data, to trigger a rehash
    for (int i = 0; i < NUM_TO_DELETE; i++) {
      // set up tempDna for current target of deletion
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      // remove the target, fail if anything goes wrong
      if (dnadb.remove(tempDna) == false) {
	cout << FAIL_STATEMENT << "remove failed the error case\n";
	return false;
      }
    }

    // attempt to remove the not-added DNA, should fail
    tempDna.m_sequence = dnaSeq[SIZE];
    tempDna.m_location = dnaLoc[SIZE];
    
    if (dnadb.remove(tempDna)) {
      cout << FAIL_STATEMENT << "remove failed the error case\n";
      return false;
    }
    
    // see if size is correct
    if (dnadb.m_currentSize > (SIZE - NUM_TO_DELETE)) {
      cout << FAIL_STATEMENT << "remove failed the error case\n";
      return false;
    }

    // confirm that deleted data can't be found
    for (int i = 0; i < NUM_TO_DELETE; i++) {
      if (getDNAIndex(dnadb, dnaSeq[i]) != -1) {
	cout << FAIL_STATEMENT << "remove failed the error case\n";
	return false;
      }
    }

    // confirm that non-deleted data can be found
    for (int i = NUM_TO_DELETE; i < (SIZE - 1); i++) {
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "remove failed the error case\n";
	return false;
      }
    }

    // passed all the tests
    cout << PASS_STATEMENT << "remove passed the error case\n";
    return true;
  }
  
  //Name: testLambdaNormal
  //Case: normal
  //Expected result: after a decent number of insertions, and then some deletions, but neither
  // being enough to trigger a rehash, lambda will return the load factor
  bool testLambdaNormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 5; // less than a fourth of the capacity, avoid rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];
      
      dnadb.insert(tempDna);
    }

    // confirm that lambda will return the correct float after populating the hash table
    if (dnadb.lambda() != (SIZE / (float)CAPACITY)) {
      cout << FAIL_STATEMENT << "lambda failed the normal case\n";
      return false;
    }

    // delete half of the data using tempDna
    for (int i = 0; i < SIZE; i += 2) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];
      
      dnadb.remove(tempDna);
    }

    // lambda's value should be the same as before, as deleted buckets count as occupied
    if (dnadb.lambda() != (SIZE / (float)CAPACITY)) {
      cout << FAIL_STATEMENT << "lambda failed the normal case\n";
      return false;
    }

    // passed all the tests
    cout << PASS_STATEMENT << "lambda passed the normal case\n";
    return true;
  }

  //Name: testLambdaNormalRehash
  //Case: normal
  //Expected result: after a decent number of insertions and deletions, with a rehash being
  // triggered, lambda will still be accurate with the current table
  bool testLambdaNormalRehash() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 2 + 1; // more than a half the capacity, trigger rehashing
    int numOccupied = 0; // used to count used buckets in the rehashed current table
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert half the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < (SIZE / 2); i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete half of the data, and confirm that lambda isn't incorrect
    for (int i = 0; i < (SIZE / 2); i += 2) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.remove(tempDna);
    }

    // insert the rest of the DNA into the database, should trigger a rehash
    for (int i = (SIZE / 2); i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];
      
      dnadb.insert(tempDna);
    }
    
    // count all occupied buckets in the new current table
    for (int i = 0; i < dnadb.m_currentCap; i++) {
      if (dnadb.m_currentTable[i] != nullptr) {
      numOccupied++;
      }
    }
    

    // see if lambda returns the correct value
    if (dnadb.lambda() == (numOccupied / (float)dnadb.m_currentCap)) {
      cout << PASS_STATEMENT << "lambda passed the normal (rehash) case\n";
      return true;
    }

    //else, failed the test
    cout << FAIL_STATEMENT << "lambda failed the normal (rehash) case\n";
    return false;
  }

  //Name: testDeletedRatioNormal
  //Case: normal
  //Expected result: after a decent number of insertions and deletions, the deleted ratio
  // will be accurate
  bool testDeletedRatioNormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 5; // less than a fourth of the capacity, avoid rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete half of the data using tempDna
    for (int i = 0; i < SIZE; i += 2) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.remove(tempDna);
    }

    // test that deletedRatio isn't innacurate
    // expecting (.5 * Size) / (size), so just .5
    if (dnadb.deletedRatio() != 0.5) {
      cout << FAIL_STATEMENT << "deletedRatio failed the normal case\n";
      return false;
    }

    // passed all the tests
    cout << PASS_STATEMENT << "deletedRatio passed the normal case\n";
    return true;
  }

  //Name: testDeletedRatioNormalRehash
  //Case: normal
  //Expected result: after a decent number of insertions and deletions, enough to trigger a
  // rehash, the deleted ratio will remain accurate and only consider the current table
  bool testDeletedRatioNormalRehash() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 751; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 2 + 1; // more than a half the capacity, trigger rehashing
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert half the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < (SIZE / 2); i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // delete half of the data, and confirm that lambda isn't incorrect
    for (int i = 0; i < (SIZE / 2); i += 2) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.remove(tempDna);
    }

    // insert the rest of the DNA into the database, should trigger a rehash
    for (int i = (SIZE / 2); i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // see if deleted ratio is correct
    // should be 0 - as no deleted dna would be transfered to the current table
    if (dnadb.deletedRatio() == 0.0) {
      cout << PASS_STATEMENT << "deletedRatio passed the normal (rehash) case\n";
      return true;
    }

    //else, failed the test
    cout << FAIL_STATEMENT << "deletedRatio failed the normal (rehash) case\n";
    return false;
  }

  //Name: testUpdateLocIdNormal
  //Case: normal
  //Expected result: a decently-populated database, after a rehash and a call to updateLocId,
  // the DNA will have an updated location
  bool testUpdateLocIdNormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 307; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 2 + 1; // trigger a rehash
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database, expect the last one (will be used as a new
    // location for one of the DNA)
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < (SIZE - 1); i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // manually trigger a rehash
    dnadb.rehashData();
    
    // update the location of the first dna inserted with the last generated location
    // using tempDna as the target
    tempDna.m_sequence = dnaSeq[0];
    tempDna.m_location = dnaLoc[0];
    
    dnadb.updateLocId(tempDna, dnaLoc[SIZE - 1]);

    // update the first dna location to avoid having to edit the following tests
    dnaLoc[0] = dnaLoc[SIZE - 1];
    
    // check that all expected DNA objects are present in the database
    for (int i = 0; i < (SIZE - 1); i++) {
      // if any inserted DNA aren't in the database, fail
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "updateLocId failed the normal case\n";
	return false;
      }

      // if any inserted DNA has an incorrect location, fail
      // item could be in the current or old table, but it must exist in one of them
      // check which by checking its sequence
      if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])] != nullptr
	  and dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_sequence == dnaSeq[i]) {
	if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "updateLocId failed the normal case\n";
	  cout << endl;
	  return false;
	}
      }
      else { // located in old table
	if (dnadb.m_oldTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "updateLocId failed the normal case\n";
	  return false;
	}
      }
    }

    // failed no tests
    cout << PASS_STATEMENT << "updateLocId passed the normal case\n";
    return true;
  }

  //Name: testUpdateLocIdError
  //Case: error
  //Expected result: attempting to change the location of an id to an invalid or duplicate
  // location will fail, and attempting to change a non-existant dna will as well
  bool testUpdateLocIdError() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 307; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    string fakeSeq = "ACTG"; // sequence not of seqLength, so will not already exist
    int fakeLoc = MINLOCID + 1; // fake location, which will be assured to be unique
    bool unique = true; // stores if a DNA sequence or location is unique or not
    const int SIZE = CAPACITY / 2 + 1; // trigger a rehash
    const hash_fn HASH = hashCode;
    const prob_t PROBING = LINEAR;

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true; // assume the new DNA is unique
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j] or dnaLoc[h] == fakeLoc) {
	  unique = false;
	}
      }
      // redo this generation if its a repeat value
      if (not unique) {
	j--;
      }
    }

    // insert all the DNA into the database, expect the last one (will be used as a new
    // location for one of the DNA)
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // fail if trying to change a non-existant DNA doesn't cause updateLocId to fail
    // update tempDna to a fake sequence/DNA first
    tempDna.m_sequence = fakeSeq;
    if (dnadb.updateLocId(tempDna, fakeLoc)) {
      cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
      return false;
    }

    // fail if it's possible to update a location with a non-unique location
    tempDna.m_sequence = dnaSeq[0];
    tempDna.m_location = dnaLoc[0];
    if (dnadb.updateLocId(tempDna, dnaLoc[1])) {
      cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
      return false;
    }

    // fail if it's possible to update a location with an invalid location
    if (dnadb.updateLocId(tempDna, (MINLOCID - 1))) {
      cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
      return false;
    }

    // check that all expected DNA objects are present in the database
    for (int i = 0; i < SIZE; i++) {
      // if any inserted DNA aren't in the database, fail
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
	return false;
      }

      // if any inserted DNA has an incorrect location, fail
      // item could be in the current or old table, but it must exist in one of them
      // check which by checking its sequence
      if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])] != nullptr
	  and dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_sequence == dnaSeq[i]) {
	if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
	  cout << endl;
	  return false;
	}
      }
      else { // located in old table
	if (dnadb.m_oldTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "updateLocId failed the error case\n";
	  return false;
	}
      }
    }
 
    // failed no tests
    cout << PASS_STATEMENT << "updateLocId passed the error case\n";
    return true;
  }

  //Name: testChangeProbPolicyNormal
  //Case: normal
  //Expected result: after calling changing the probing policy, the current table will use that
  // policy after a rehash has been triggered
  bool testChangeProbPolicyNormal() {
    // define constants to create a database and generate sample DNA for it
    const int CAPACITY = 101; // size of the array
    string currentSequence = ""; // current sequence being added to the array
    int seqLength = 10; // the number of chars in a sequence
    DNA tempDna; // used to add dna to the database
    const int SIZE = 52; // more than half of the capacity to trigger a rehash
    const hash_fn HASH = hashCode;
    const prob_t PROBING = QUADRATIC;
    const prob_t NEW_POLICY = DOUBLEHASH; // probing policy to be changed into
    bool unique = true; // represents if currend data in DNA is unique

    // arrays to store DNA objects' data
    string dnaSeq[SIZE]; // sequences
    int dnaLoc[SIZE] = {0}; // locations

    // create the database
    DnaDb dnadb(CAPACITY, HASH, PROBING);

    // create the Random objects for generating DNA object
    Random randomSeq(1, 4); // for sequences; 1 thru. 4 because there's A, C, T and G
    Random randomLoc(MINLOCID, MAXLOCID); // for locations

    // generate data for the DNAs
    for (int j = 0; j < SIZE; j++) {
      // for the length of a sequence, add an A, C, T or G to the current sequence
      for (int i = 0; i < seqLength; i++) {
	switch (randomSeq.getRandNum()) {
	case 1:
	  currentSequence += "A";
	  break;
	case 2:
	  currentSequence += "C";
	  break;
	case 3:
	  currentSequence += "G";
	  break;
	case 4:
	  currentSequence += "T";
	}
      }
      // add the sequence to the array
      dnaSeq[j] = currentSequence;
      currentSequence = "";

      // add a random location to its array
      dnaLoc[j] = randomLoc.getRandNum();

      // check that the new sequence and location are unique, else generate again (decrement j)
      unique = true;
      for (int h = 0; h < j; h++) {
	if (dnaSeq[h] == dnaSeq[j] or dnaLoc[h] == dnaLoc[j]) {
	  unique = false;
	}
      }
      if (not unique) {
	j--;
      }
    }

    // attempt to change the probing policy
    dnadb.changeProbPolicy(NEW_POLICY);
    
    // insert all the DNA into the database
    // modify the temp DNA object, then add it to the database
    for (int i = 0; i < SIZE; i++) {
      tempDna.m_sequence = dnaSeq[i];
      tempDna.m_location = dnaLoc[i];

      dnadb.insert(tempDna);
    }

    // check that all expected DNA objects are present in the database
    for (int i = 0; i < SIZE; i++) {
      // if any inserted DNA aren't in the database, fail
      if (getDNAIndex(dnadb, dnaSeq[i]) == -1) {
	cout << FAIL_STATEMENT << "changeProbPolicy failed the normal case\n";
	return false;
      }

      // if any inserted DNA has an incorrect location, fail
      // item could be in the current or old table, but if the above test passed, then it
      // must exist in one of them - check which by checking its sequence
      if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])] != nullptr
	  and dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_sequence == dnaSeq[i]) {
	if (dnadb.m_currentTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "changeProbPolicy failed the normal case\n";
	  cout << endl;
	  return false;
	}
      }
      else { // located in old table
	if (dnadb.m_oldTable[getDNAIndex(dnadb, dnaSeq[i])]->m_location != dnaLoc[i]) {
	  cout << FAIL_STATEMENT << "changeProbPolicy failed the normal case\n";
	  return false;
	}
      }
    }

    // check that the size member variable updated correctly after all the insertions
    // data is random and hash function isn't perfect, so can't guarantee the value of any
    // of the sizes of the tables, but the old table should (likely) still contain dna,
    // and the current table should have between 0 and the size
    if (dnadb.m_currentSize >= SIZE or dnadb.m_currentSize == 0 or dnadb.m_oldSize == 0) {
      cout << FAIL_STATEMENT << "changeProbPolicy failed the normal case\n";
      return false;
    }
    
    // confirm that the current table has the new policy
    if (dnadb.m_currProbing != NEW_POLICY) {
      cout << FAIL_STATEMENT << "changeProbPolicy failed the normal case\n";
    }      

    // failed no tests
    cout << PASS_STATEMENT << "changePropPolicy passed the normal case\n";
    return true;  
  }
  
  //********HELPER TESTER FUNCTIONS********//

  //Name:getDNAIndex
  //Preconditions: is given a valid DnaDb object and sequence
  //Postconditions: returns the index of the DNA object in the hash table. returns -1 if the
  // DNA cannot be found
  int getDNAIndex(DnaDb &dnadb, string sequence) {
    // start searching indices
    int index = -1; // index of the hash table to insert the dna into
    int i = 0; // used by the probing functions to adjust OG hash index if already occupied
    
    // determine the index of the target dna
    while (index == -1) {
      // calculate the index of the target dna
      // equation is based on the probing policy
      switch (dnadb.m_currProbing) {
      case QUADRATIC:
	index = ((dnadb.hash(sequence) + i * i) % (unsigned int)dnadb.m_currentCap);
	break;
	
      case DOUBLEHASH:
	index = ((dnadb.hash(sequence) % (unsigned int)dnadb.m_currentCap)
		 + i * ((unsigned int)11 - (dnadb.hash(sequence) % (unsigned int)11)))
	  % (unsigned int)dnadb.m_currentCap;
	break;
	
      case LINEAR:
	index = ((dnadb.hash(sequence) + i) % (unsigned int)dnadb.m_currentCap);
      }
            
      // if the index holds DNA, and it doesn't contain the target DNA, continue
      if (dnadb.m_currentTable[index] != nullptr and dnadb.m_currentTable[index]->m_used
	  and dnadb.m_currentTable[index]->m_sequence != sequence) {
	index = -1;
      }
      
      // if the index holds DNA, and it is the target DNA, return that index
      else if (dnadb.m_currentTable[index] != nullptr and dnadb.m_currentTable[index]->m_used
	       and dnadb.m_currentTable[index]->m_sequence == sequence) {
	return index;
      }
     
      // if the index used to have DNA, continue searching
      else if (dnadb.m_currentTable[index] != nullptr and
	       dnadb.m_currentTable[index]->m_used == false) {
	index = -1;
      }

      // if the index doesn't have DNA, the target is not in the current table
      // check the old table
      else {
	// only check the old table if it exists!
	if (dnadb.m_oldTable != nullptr) {
	  index = -1; // reset index, go thru while loop
	  i = 0; // reset i for the old table search
	  while (index == -1) {
	    // calculate the index of the target dna
	    // equation is based on the probing policy
	    switch (dnadb.m_oldProbing) {
	    case QUADRATIC:
	      index = ((dnadb.hash(sequence) + i * i) % (unsigned int)dnadb.m_oldCap);
	      break;
	      
	    case DOUBLEHASH:
	      index = ((dnadb.hash(sequence) % (unsigned int)dnadb.m_oldCap)
		       + i * ((unsigned int)11 - (dnadb.hash(sequence) % (unsigned int)11)))
		% (unsigned int)dnadb.m_oldCap;
	      break;
	      
	    case LINEAR:
	      index = ((dnadb.hash(sequence) + i) % (unsigned int)dnadb.m_oldCap);
	    }
	    
	    // if the index holds DNA, and it doesn't contain the target DNA, continue
	    if (dnadb.m_oldTable[index] != nullptr and dnadb.m_oldTable[index]->m_used
		and dnadb.m_oldTable[index]->m_sequence != sequence) {
	      index = -1;
	    }
	    
	    // if the index holds DNA, and it is the target DNA, return that index
	    else if (dnadb.m_oldTable[index] != nullptr and dnadb.m_oldTable[index]->m_used
		     and dnadb.m_oldTable[index]->m_sequence == sequence) {
	      // only count it as found if it wasn't already transfered, otherwise it could
	      // be that it was deleted, but only after a transfer
	      if (dnadb.m_transferIndex <= index) {
		return index;
	      }
	      else {
		index = -1;
	      }
	    }
	    
	    // if the index used to have DNA, continue searching
	    else if (dnadb.m_oldTable[index] != nullptr and
		     dnadb.m_oldTable[index]->m_used == false) {
	      index = -1;
	    }
	    
	    // if the index doesn't have DNA, the target is not in the database
	    else {
	      return -1;
	    }
	    i++;
	  }
	}
      }
      i++;
    }
    
    // command should not reach here. Return -1 if so, however
    return -1;
  }
};

int main() {
  Tester tester;

  
  tester.testConstructorNormal();

  tester.testHashFunctionNormal();

  tester.testInsertNormal();  
  tester.testInsertNormalRehashPartial();
  tester.testInsertError();
  
  tester.testGetDNANormalCollide();
  tester.testGetDNANormal();
  tester.testGetDNAError();
  
  tester.testRemoveNormalUnique();
  tester.testRemoveNormalCollide();
  tester.testRemoveNormalRehash();
  tester.testRemoveError();
  
  tester.testLambdaNormal();
  tester.testLambdaNormalRehash();

  tester.testDeletedRatioNormal();
  tester.testDeletedRatioNormalRehash();

  tester.testUpdateLocIdNormal();
  tester.testUpdateLocIdError();

  tester.testChangeProbPolicyNormal();
    
  return 0;
}
