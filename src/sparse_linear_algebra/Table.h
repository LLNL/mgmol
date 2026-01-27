// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

/*!
C++ header file for Hash table class.
Adapted from pARMS code. (original version by Z. Li)
-DOK
*/

#ifndef _TABLE_H_
#define _TABLE_H_

#include "Timer.h"

#include <vector>

#define USE_POWERS2 1
// #define COUNT_LINKS 1

typedef struct Slot
{
    struct Slot* link; //!< pointer to a next entry in this slot.
    int key; //!< key in a pair (key,value)
    int value; //!< value in a pair (key, value)
} Slot;

/*! \struct parms_Table_
  \brief parms_Table_ structure.
 */
class Table
{
    // static  Timer  destructor_tm_;
    // static  Timer  reset_tm_;
    // static  Timer  get_value_tm_;
    // static  Timer  insert_tm_;

    Slot** Slots_; //!< array of a pointer array to struct Slot
    int space_; //!< size of array Slots.
    int space_power_minus1_;
    int size_; //!< number of pairs in the table.

    Slot* slot_ptr_; //!< pointer to slot in memory
    std::vector<Slot*> slot_storage_; //!< storage for slots
    int capacity_; //!< size of slot_storage

    void freeStorage();

    int computeIndex(const int key)
    {
#ifdef USE_POWERS2
        return (int)(key & space_power_minus1_);
#else
        return (int)(key % space_);
#endif
    }

#ifdef COUNT_LINKS
    short maxlink_;
    short init_size_;
#endif

public:
    Table(const int); // constructor
    ~Table(); // Destructor

    void* get_value(int); /* get value corresponding to key */
    int insert(int, int); /* insert (key, value) pair */
    int insert(const std::vector<int>& keys); /* insert a set of keys */
    int insert(int); /* insert (key, current_size) pair */
    void reset(); /* reset size of table -- table will be overwritten by
                     subsequent inserts */

    /* return current size of Table */
    int get_size() const { return size_; }
    static void printTimers(std::ostream& /*os*/)
    {
        //         insert_tm_.print(os);
        //         get_value_tm_.print(os);
        //         reset_tm_.print(os);
        //         destructor_tm_.print(os);
    }
};

#endif
