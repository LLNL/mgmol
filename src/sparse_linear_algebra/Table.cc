// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*--------------------------------------------------------------------
  Table   : create a hash table.
  ~Table     : free the memory for the table.
  get_value      : get the corresponding value for a given key.
  get_size  : get the total number of entries in the table.
  insert      : put the pair (key, value) into the table.

  ------------------------------------------------------------------*/
#include "Table.h"
#include "MPIdata.h"

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>

// Timer   Table::destructor_tm_("Table::destructor");
// Timer   Table::reset_tm_("Table::reset");
// Timer   Table::get_value_tm_("Table::get_value");
// Timer   Table::insert_tm_("Table::insert");

#ifdef COUNT_LINKS
#include "MGmol_MPI.h"
#endif

/**
 * hash table constructor
 *
 * @param tsize  The number of entries stored in the table. NOTE: tsize = 0 is
 * OK.
 *
 */
Table::Table(const int tsize)
{
#ifdef COUNT_LINKS
    maxlink_   = 0;
    init_size_ = tsize;
#endif

    int i, lwr;
#ifdef USE_POWERS2
    static int powers[] = { 9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
    static int powersof2[]
        = { 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536, 131072, 262144 };
    for (i = 1; lwr = (powersof2[i] * 2) / 3, lwr < tsize; i++)
        ;
    const int space     = powersof2[i - 1];
    space_power_minus1_ = pow(2, powers[i - 1]) - 1;
#else
    static int primes[] = { 389, 389, 769, 1543, 3079, 6151, 12289, 24593,
        49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469,
        12582917, 25165843, 50331653, 100663319, 201326611, 402653189,
        805306457, 1610612741 };

    /* find the maximum prime less than tsize */
    for (i = 1; lwr = (primes[i] * 2) / 3, lwr < tsize; i++)
        ;
    const int space = primes[i - 1];
#endif

    /* allocate the space for private member variable */
    size_     = 0;
    space_    = space;
    capacity_ = space;
    // allocate array of slot pointers and initialize slot pointers to NULL
    Slots_ = new Slot*[space_];
    for (i = 0; i < space; i++)
    {
        Slots_[i] = nullptr;
    }

    // Allocate memory for storing data
    struct Slot* newstorage = new Slot[capacity_];
    slot_storage_.push_back(newstorage);
    // initialize slot_ptr_
    slot_ptr_ = newstorage;
}

/**
 * Get the corresponding value for a given key.
 *
 * If return NULL, then the entry with key is not in the table,
 * otherwise return a pointer to the value.
 *
 * @param key   The key value.
 *
 * @return A pointer to the value.
 */
void* Table::get_value(int key)
{
    //  get_value_tm_.start();
    const int index = computeIndex(key);

#ifdef COUNT_LINKS
    short counter = 0;
#endif

    struct Slot* const slot = Slots_[index];
    struct Slot* p;
    for (p = slot; p; p = p->link)
    {
        if (p->key == key)
        {
#ifdef COUNT_LINKS
            maxlink_ = (maxlink_ > counter) ? maxlink_ : counter;
#endif
            return &p->value;
        }
#ifdef COUNT_LINKS
        counter++;
#endif
    }
    return nullptr;
    //  get_value_tm_.stop();
}

/**
 * Insert key into the table.
 * Here, the associated value is assumed to be the current size
 * of the table.
 * NOTE: Does NOT check if the key already exists. It will reinsert key if it
 * already exists.
 *
 * @param key 	The key of the pair.
 *
 * @return 0 on success.
 */

int Table::insert(int key)
{
//  insert_tm_.start();
// reallocate storage if needed
#ifdef USE_POWERS2
    if ((size_ & space_power_minus1_) == 0)
#else
    if (size_ % capacity_ == 0)
#endif
    {
        struct Slot* newstorage = new Slot[capacity_];
        slot_storage_.push_back(newstorage);
        slot_ptr_ = newstorage;
    }

    slot_ptr_->key   = key;
    slot_ptr_->value = size_;

    const int index = computeIndex(key);
    slot_ptr_->link = Slots_[index];
    Slots_[index]   = slot_ptr_;
    size_++;
    slot_ptr_++;

    //  insert_tm_.stop();
    return 0;
}

/**
 * Insert multiple keys into the table.
 * Here, the associated value is assumed to be the current size
 * of the table. If the key already exists, then nothing is done.
 *
 * @param key 	The key of the pair.
 *
 * @return 0 on success.
 */
int Table::insert(const std::vector<int>& keys)
{
    //  insert_tm_.start();

    // reallocate storage if needed
    const int diff = static_cast<int>(keys.size()) > (capacity_ - size_)
                         ? (capacity_ - size_)
                         : 0;
    std::vector<int>::const_iterator key;
    if (diff)
    {
        // first fill up current storage
        key                                     = keys.begin();
        std::vector<int>::const_iterator endkey = keys.begin() + diff;
        while (key != endkey)
        {
            struct Slot* p  = slot_ptr_;
            const int index = computeIndex(*key);
            p->key          = *key;
            p->value        = size_;
            p->link         = Slots_[index];
            Slots_[index]   = p;
            size_++;
            slot_ptr_++;
            key++;
        }

        struct Slot* newstorage = new Slot[capacity_];
        slot_storage_.push_back(newstorage);
        slot_ptr_ = newstorage;
    }

    // insert remaining entries
    key = keys.begin() + diff;
    while (key != keys.end())
    {
        struct Slot* p  = slot_ptr_;
        const int index = computeIndex(*key);
        p->key          = *key;
        p->value        = size_;
        p->link         = Slots_[index];
        Slots_[index]   = p;
        size_++;
        slot_ptr_++;
        key++;
    }

    //  insert_tm_.stop();
    return 0;
}

/**
 * Put the pair (key, value) into the table.
 *
 * @param key 	The key of the pair.
 * @param val 	The value of the pair.
 *
 * @return 0 on success.
 */
int Table::insert(int key, int value)
{
    struct Slot* p;
    //  insert_tm_.stop();
    const int index         = computeIndex(key);
    struct Slot* const slot = Slots_[index];
    for (p = slot; p; p = p->link)
    {
        if (p->key == key)
        {
            break;
        }
    }
    if (p == nullptr)
    { /* new entry */
        // reallocate storage if needed
#ifdef USE_POWERS2
        if ((size_ & space_power_minus1_) == 0)
#else
        if (size_ % capacity_ == 0)
#endif
        {
            struct Slot* newstorage = new Slot[capacity_];
            slot_storage_.push_back(newstorage);
            slot_ptr_ = newstorage;
        }

        p             = slot_ptr_;
        p->key        = key;
        p->value      = value;
        p->link       = Slots_[index];
        Slots_[index] = p;
        size_++;
        slot_ptr_++;
    }
    else
    {
        p->key   = key;
        p->value = value;
    }
    //  insert_tm_.stop();
    return 0;
}

/**
 * Destructor: Free the memory for the table object.
 *
 *
 * @return 0 on success.
 */
Table::~Table()
{
#ifdef COUNT_LINKS
    MGmol_MPI& mmpi(*(MGmol_MPI::instance()));
    //  mmpi.allreduce(&maxlink_, 1, MPI_MAX);
    int maxsize = size_;
    //  mmpi.allreduce(&maxsize, 1, MPI_MAX);
    if (onpe0 && (maxsize > space_ || maxlink_ > 9))
        std::cout << "maxlink=" << maxlink_ << ", maxsize=" << maxsize
                  << ", space=" << space_ << ", init_size_=" << init_size_
                  << std::endl;
#endif

    // destructor_tm_.start();
    freeStorage();
    delete[] Slots_;
    // destructor_tm_.stop();
}

void Table::reset()
{
    //  reset_tm_.start();
    for (int i = 0; i < space_; i++)
    {
        Slots_[i] = nullptr;
    }
    size_ = 0;
    freeStorage();
    struct Slot* newstorage = new Slot[capacity_];
    slot_storage_.push_back(newstorage);
    // initialize slot_ptr_
    slot_ptr_ = newstorage;
    //  reset_tm_.stop();
}

// free allocated storage
void Table::freeStorage()
{
    std::vector<Slot*>::iterator it = slot_storage_.begin();
    while (it != slot_storage_.end())
    {
        delete[] * it;
        it++;
    }
    slot_storage_.clear();
}
