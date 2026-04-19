// CMSC 341 - Spring 25 - Project 4

/*
  File: dnadb.cpp
  Author: Andrew Greenberg
  Date: 4/18/25
  Desc: file containing member function definitions for the DnaDb class
 */

#include "dnadb.h"

//Name: DnaDb (constructor)
//Preconditions: passed size should be a prime number between MINPRIME and MAXPRIME
//Postconditions: initializes all member variables and creates memory for the current table
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
  // initialize the size of the array
  // check that the user passed a prime number, if they didn't, then make size the next prime
  if (isPrime(size) and size >= MINPRIME and size <= MAXPRIME) {
    m_currentCap = size;
  }
  else {
    m_currentCap = findNextPrime(size);
  }

  // initialize the remaining member variables
  // current hash table
  m_hash = hash;
  m_newPolicy = probing;
  
  m_currProbing = probing;
  m_currentTable = new DNA*[m_currentCap];

  m_currentSize = 0;
  m_currNumDeleted = 0;

  // set the DNA pointers in the current table to nullptr
  for (int i = 0; i < m_currentCap; i++) {
    m_currentTable[i] = nullptr;
  }

  m_oldTable = nullptr;
  m_oldCap = 0;
  m_oldSize = 0;
  m_oldNumDeleted = 0;
  m_oldProbing = QUADRATIC;

  m_transferIndex = 0;
}

DnaDb::~DnaDb(){
  // if applicable, go through the currentTable and delete all DNA after resetting their data
  if (m_currentTable != nullptr) {
    // for the number of possible entries
    // check if there is DNA in the slot, delete it if so
    for (int i = 0; i < m_currentCap; i++) {
      if (m_currentTable[i] != nullptr) {
	m_currentTable[i]->m_sequence = "";
	m_currentTable[i]->m_location = 0;
	m_currentTable[i]->m_used = false;

	delete m_currentTable[i];
	m_currentTable[i] = nullptr;
      }
    }
  }
  // if applicable, go through the oldTable and delete all DNA after resetting their data
  if (m_oldTable != nullptr) {
    // for the number of possible entries
    // check if there is DNA in the slot, delete it if so
    for (int i = 0; i < m_oldCap; i++) {
      if (m_oldTable[i] != nullptr) {
	m_oldTable[i]->m_sequence = "";
	m_oldTable[i]->m_location = 0;
	m_oldTable[i]->m_used = false;
	
	delete m_oldTable[i];
	m_oldTable[i] = nullptr;
      }
    }
  }
  
  // deallocate memory for the tables/arrays, if they hold data
  if (m_currentTable != nullptr) {
    delete[] m_currentTable;
    m_currentTable = nullptr;
  }
  if (m_oldTable != nullptr) {
    delete[] m_oldTable;
    m_oldTable = nullptr;
  }
  
  // reset all remaining member variables
  m_hash = nullptr;
  m_newPolicy = DEFPOLCY;
  m_currProbing = DEFPOLCY;
  m_currentCap = 0;
  m_currentSize = 0;
  m_currNumDeleted = 0;
  m_oldCap = 0;
  m_oldSize = 0;
  m_oldNumDeleted = 0;
  m_oldProbing = DEFPOLCY;
  m_transferIndex = 0;
}

//Name: changeProbPolicy
//Preconditions: N/A
//Postconditions: if the new policy is different from the current policy, then m_newPolicy
// will be updated and the next beginning of a rehash will have the new probing policy
void DnaDb::changeProbPolicy(prob_t policy){
  if (policy != m_currProbing) {
    m_newPolicy = policy;
  }
}

//Name: insert
//Preconditions: the DNA object must not be a duplicate and its ID value must be valid
//Postconditions: the DNA object is inserted into the current table, transfers 25% of the
// current table to the new table, and returns true. If insertion fails, returns false
bool DnaDb::insert(DNA dna){
  // check that the DNA is unique and contains a valid ID, else cancel and return false
  if ((getDNA(dna.m_sequence, dna.m_location).m_sequence != "")
      or (dna.m_location < MINLOCID or dna.m_location > MAXLOCID)) {
    return false;
  }
  
  int index = -1; // index of the hash table to insert the dna into
  int i = 0; // used by the probing functions to adjust OG hash index if already occupied
  
  // determine the index of the new dna
  while (index == -1) {
    // calculate the index of the new dna
    // equation is based on the probing policy
    switch (m_currProbing) {
    case QUADRATIC:
      index = ((hash(dna.m_sequence) + i * i) % (unsigned int)m_currentCap);
      break;
      
    case DOUBLEHASH:
      index = ((hash(dna.m_sequence) % (unsigned int)m_currentCap)
	       + i * ((unsigned int)11 - (hash(dna.m_sequence) % (unsigned int)11)))
	% (unsigned int)m_currentCap;
      break;
      
    case LINEAR:
      index = ((hash(dna.m_sequence) + i) % (unsigned int)m_currentCap);
    }
    
    // if the index is open, exit the loop, otherwise try again and increment i
    if (m_currentTable[index] != nullptr) {
      index = -1; // stay in the while loop
    }
    i++;
  }
  
  // add the dna object to the hash table and correct the size
  m_currentTable[index] = new DNA(dna.m_sequence, dna.m_location, true);
  m_currentSize++;
  
  // rehash if the load factor is greater than 0.5
  if (lambda() > 0.5) {
    rehashData();
  }
  return true;
}


//Name: remove
//Preconditions: the target DNA exists in the database
//Postconditions: the target DNA will be tagged as deleted (m_used set to false) and will be
// officially deleted in a rehash or with the destruction of the database. A rehash will be
// triggered if the number of deleted buckets is >=80% of the number of occupied buckets.
// returns true if successful, else false
bool DnaDb::remove(DNA dna){
  bool result = false; // reflects if the removal was successful (true) or not (false)
  
  // locate the item, return false if not found, else true after changing its m_used to false

  // note the indices of where the dna is in the current and old table
  int index_current = locateDNA(m_currentTable, dna.m_sequence, m_currProbing, m_currentCap);
  int index_old = locateDNA(m_oldTable, dna.m_sequence, m_oldProbing, m_oldCap);

  // if the dna doesn't exist, can't remove it
  if (index_current == -1 and index_old == -1) {
    result = false;
  }
  // if the dna is in the current table, remove it
  if (index_current != -1) {
    m_currentTable[index_current]->m_used = false;
    m_currentSize--;
    m_currNumDeleted++;
    result = true;
  }
  // if it was in the old table, remove it
  if (index_old != -1) {
    m_oldTable[index_old]->m_used = false;
    m_oldSize--;
    m_oldNumDeleted++;
    result = true;
  }
  
  // check if a rehash is necessary (deleted buckets are >80% of occupied buckets
  // rehash if it is
  if (deletedRatio() > 0.8) {
    rehashData();
  }
  
  // conclude by returning the result of the removal
  return result;
}

//Name: getDNA
//Preconditions: the passed sequence and location belong to the same DNA
//Postconditions: if found, the DNA object is returned, otherwise an empty object is
const DNA DnaDb::getDNA(string sequence, int location) const{
  // define an empty DNA object to return if the target isn't found
  DNA emptyDNA;

  // note the indexes of where the dna is in the current and old table
  int index_current = locateDNA(m_currentTable, sequence, m_currProbing, m_currentCap);
  int index_old = locateDNA(m_oldTable, sequence, m_oldProbing, m_oldCap);

  // if it was in the current table, return the dna object
  if (index_current != -1) {
    return *(m_currentTable[index_current]);
  }
  // if it was in the old table, return it
  else if (index_old != -1) {
    return *(m_oldTable[index_old]);
  }
  
  // if it wasn't in either, return the empty object
  return emptyDNA;
}

//Name: updateLocId
//Preconditions: given DNA exists in the current or old table and the new location is unique
//Postconditions: if found, the DNA object's location will be changed and the function will
// return true, else it fails and returns false. Fails if the target has been soft/easy deleted
bool DnaDb::updateLocId(DNA dna, int location){
  bool uniqueLoc = true; // states if the new location is unique. must be true to update

  // check the current and old tables to see if any DNA exist with the location to update with
  // currentTable must exist, if it somehow doesn't then no DNA exists to update
  if (m_currentTable == nullptr) {
    return false;
  }

  // if the new location is invalid, fail
  if (location < MINLOCID or location > MAXLOCID) {
    return false;
  }
  
  // for each bucket in the current table, check if the new location is unique
  for (int i = 0; i < m_currentCap; i++) {
    if (m_currentTable[i] != nullptr and m_currentTable[i]->m_location == location) {
      uniqueLoc = false;
    }
  }

  // confirm that an old table exists
  if (m_oldTable != nullptr) {
    // check the old table's DNA to confirm that the new location is unique
    for (int i = 0; i < m_oldCap; i++) {
      if (m_oldTable[i] != nullptr and m_oldTable[i]->m_location == location) {
	uniqueLoc = false;
      }
    }
  }

  // if the new location is a repeat, then the update must fail
  if (uniqueLoc == false) {
    return false;
  }

  // note the indexes of where the dna is in the current and old table
  int index_current = locateDNA(m_currentTable, dna.m_sequence, m_currProbing, m_currentCap);
  int index_old = locateDNA(m_oldTable, dna.m_sequence, m_oldProbing, m_oldCap);

  // if it was in the current table, update its location
  if (index_current != -1) {
    m_currentTable[index_current]->m_location = location;
    return true;
  }
  // if it was in the old table, update it
  else if (index_old != -1) {
    m_oldTable[index_old]->m_location = location;
    return true;
  }

  //update failed
  return false;
}


//Name: lambda
//Preconditions: m_currentSize, m_currNumDeleted and m_currentCap must be correct
//Postconditions: returns the load factor (almost size / capacity) of the current table. Note
// that a deleted node does count as an occupied bucket
float DnaDb::lambda() const {
  float numOccupied = m_currentSize + m_currNumDeleted; // number of live and deleted nodes in
                                                        // the current table
  return (numOccupied / m_currentCap);
}

//Name: deletedRatio
//Preconditions: m_currentSize and m_currNumDeleted must be correct
//Postconditions: returns the ratio of deleted DNA to occupied buckets
float DnaDb::deletedRatio() const {
  float numOccupied = m_currentSize + m_currNumDeleted; // number of live and deleted nodes in
                                                        // the current table
  return (m_currNumDeleted / numOccupied);
}


/*********HELPER FUNCTIONS*************/

//Name: hash
//Preconditions: is given the unique sequence of DNA
//Postconditions: will return an int to be % to map to a specific index of an array. Given the
// same string input, it will give the same int output
int DnaDb::hash(const string sequence) const {
  unsigned int value = 0;
  const int MAGIC_PRIME = 33; // magic number from the textbook

  // for each character in the sequence, use its ASCII value to add to value
  // the order of and which character it is does impact the final value
  for (long unsigned int i = 0; i < sequence.length(); i++) {
    value = (value * MAGIC_PRIME) + sequence[i];
  }
  return value;
}

//Name: rehashData
//Preconditions: a rehash should be in order (either the load factor is greater than 0.5,
// the number of deleted buckets is more than 80% of the total buckets, or a rehash was already
// triggered
//Postconditions: 25% at a time, the DNA in the old table will be transfered to the current
// table or, if it is the first rehash, it will switch current and old table, then transfer 25%
void DnaDb::rehashData() {
  int tempTransferIndex = m_transferIndex; // to store the starting value of the transfer index
  
  // if this is the start of a new rehash, accept the new probing policy, and swap the
  // old and current table and respective variables 
  if (m_transferIndex == 0) {
    m_oldProbing = m_currProbing;
    m_currProbing = m_newPolicy;
    
    m_oldCap = m_currentCap;
    m_currentCap = findNextPrime(4 * (m_currentSize - m_currNumDeleted));
  
    m_oldSize = m_currentSize;
    m_currentSize = 0;
    
    m_oldNumDeleted = m_currNumDeleted;
    m_currNumDeleted = 0;

    m_oldTable = m_currentTable;
    m_currentTable = nullptr;
    m_currentTable = new DNA*[m_currentCap];
        
    // initialize DNA pointers to nullptr
    for (int j = 0; j < m_currentCap; j++) {
      m_currentTable[j] = nullptr;
    }
  }
  
  // transfer the next 25% of the old table over to the current table
  // capacity of the table is prime, so 4 won't divide evenly. As such, add one to
  // effectively ceil every divide
  for (; m_transferIndex <= (tempTransferIndex + (m_oldCap / 4) + 1);
       m_transferIndex++) {
  
    // rounded when calculating the upper bound of the index, so proceed iff index is valid
    if (m_transferIndex < m_oldCap) {
  
      // if the current index contains live data, transfer it to the currentTable
      if (m_oldTable != nullptr and m_oldTable[m_transferIndex] != nullptr
	  and m_oldTable[m_transferIndex]->m_used) {

	// insert to curr table - find index and then just point to this, then yeah
	int index = -1; // index of the hash table to insert the dna into
	int i = 0; // used by the probing functions to adjust OG hash index if already used
	
	// determine the index of the new dna
	while (index == -1) {

	  // calculate the index of the new dna
	  // equation is based on the probing policy
	  switch (m_currProbing) {
	  case QUADRATIC:
	    index = ((hash(m_oldTable[m_transferIndex]->m_sequence) + i * i)
		      % (unsigned int)m_currentCap);
	    break;
	    
	  case DOUBLEHASH:
	    index = ((hash(m_oldTable[m_transferIndex]->m_sequence)
		       % (unsigned int)m_currentCap)
		      + i * ((unsigned int)11 - (hash(m_oldTable[m_transferIndex]->m_sequence)
						 % (unsigned int)11)))
	      % (unsigned int)m_currentCap;
	    break;
	    
	  case LINEAR:
	    index = ((hash(m_oldTable[m_transferIndex]->m_sequence) + i)
		      % (unsigned int)m_currentCap);
	  }

	  // if the index is open, exit the loop, otherwise try again and increment i
	  if (m_currentTable[index] != nullptr and m_currentTable[index]->m_used) {
	    index = -1; // stay in the while loop
	  }
	  i++;
	}

	// add the dna object to the hash table and correct the size
	m_currentTable[index] = new DNA(m_oldTable[m_transferIndex]->m_sequence,
					 m_oldTable[m_transferIndex]->m_location, true);
	m_currentSize++;
      }      
    }
  }
  // if done transfering, then reinitialize oldTable and reset m_transferIndex
  if (m_transferIndex >= m_oldCap) {
    //copied code from destructor - if it ain't broke don't fix it
    if (m_oldTable != nullptr) {
      // for the number of possible entries
      // check if there is DNA in the slot, delete it if so
      for (int i = 0; i < m_oldCap; i++) {
	if (m_oldTable[i] != nullptr) {
	  m_oldTable[i]->m_sequence = "";
	  m_oldTable[i]->m_location = 0;
	  m_oldTable[i]->m_used = false;

	  delete m_oldTable[i];
	  m_oldTable[i] = nullptr;
	}
      }
    }
    // deallocate memory for the tables/arrays, if they hold data
    if (m_oldTable != nullptr) {
      delete[] m_oldTable;
    }
    //end of destructor copied code
    m_oldTable = nullptr;
    m_oldCap = 0;
    m_oldSize = 0;
    m_oldNumDeleted = 0;
    m_oldProbing = QUADRATIC;
    
    m_transferIndex = 0;
  }
}

//Name: locateDNA
//Preconditions: N/A
//Postconditions: will search the array for the DNA given its sequence. Returns its index if
// found, -1 otherwise. If the DNA is found, but was deleted (m_used == false), counts as not
// being found.
int DnaDb::locateDNA(DNA** table, string seq, prob_t probing, int capacity) const {
  // exit if the table is nullptr, no DNA found
  if (table == nullptr) {
    return -1;
  }
  
  // start searching indices
  int index = -1; // index of the hash table to where the dna would belong
  int i = 0; // used by the probing functions to adjust OG hash index if already occupied
  
  // determine the index of where the target dna should be
  do {
    // calculate the index of the target dna
    // equation is based on the probing policy
    switch (probing) {
    case QUADRATIC:
      index = ((hash(seq) + i * i) % (unsigned int)capacity);
      break;

    case DOUBLEHASH:
      index = ((hash(seq) % (unsigned int)capacity)
	       + i * ((unsigned int)11 - (hash(seq) % (unsigned int)11)))
	% (unsigned int)capacity;
      break;

    case LINEAR:
      index = ((hash(seq) + i) % (unsigned int)capacity);
    }

    // if the index holds DNA, and it doesn't contain the target DNA, continue
    if (table[index] != nullptr and table[index]->m_used
	and table[index]->m_sequence != seq) {
      index = -1;
    }

    // if the index holds DNA, and it is the target DNA, update its location
    else if (table[index] != nullptr and table[index]->m_used
	               and table[index]->m_sequence == seq) {
      return index;
    }

    // if the index used to have DNA, continue searching
    else if (table[index] != nullptr and
	     table[index]->m_used == false) {
      index = -1;
    }

    // DNA wasn't found. return -1
    else {
      return -1;
    }
    
    i++;
  } while (index == -1);

  // control shouldn't reach here. failed to locate the dna if it did, however
  return -1;
}



void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
	  cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}
