// Copyright (c) 2017, Lawrence Livermore National Security, LLC and
// UT-Battelle, LLC.
// Produced at the Lawrence Livermore National Laboratory and the Oak Ridge
// National Laboratory.
// LLNL-CODE-743438
// All rights reserved.
// This file is part of MGmol. For details, see https://github.com/llnl/mgmol.
// Please also read this link https://github.com/llnl/mgmol/LICENSE

#ifndef MGMOL_DotProductManagerFACTORY_H
#define MGMOL_DotProductManagerFACTORY_H

#include "DotProductDiagonal.h"
#include "DotProductSimple.h"
#include "DotProductWithDM.h"
#include "DotProductWithInvS.h"

template <class T>
class DotProductManagerFactory
{
public:
    static DotProductManager<T>* create(const short type)
    {
        DotProductManager<T>* dot_product_manager = nullptr;
        switch (type)
        {
            case 0:
                dot_product_manager = new DotProductDiagonal<T>();
                break;
            case 1:
                dot_product_manager = new DotProductWithInvS<T>();
                break;
            case 2:
                dot_product_manager = new DotProductWithDM<T>();
                break;
            case 3:
                dot_product_manager = new DotProductSimple<T>();
                break;
            default:
                std::cerr << "DotProductManager* create() --- option invalid\n";
        }
        return dot_product_manager;
    }
};

#endif
