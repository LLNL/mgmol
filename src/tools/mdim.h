/*=============================================================================

                            M D I M

 =============================================================================

           ALLOCATION MEMOIRE DE STRUCTURES INDICEES

    DIM1   : allocation d'un vecteur de cellules de type donné
    DIM2   : allocation d'un tableau de cellules de type donné

    La mise à zéro est automatiquement faite lors de l'allocation
    mémoire. S'il faut imposer une valeur constante, utiliser les
    procédures  *init  de BLAS_1 qui imposent une valeur. (Pour les
    matrices, se rappeler que les éléments sont aussi stockés en
    ligne, même si l'accès se fait par deux indices.)


    Auteurs :    P. Anderson, G. Anderson
    -------      Advanced C Tips and Techniques
                 Hayden Books
                 Indianapolis, 1989

                 Rachid Touzani
                 Département de mathématique
                 Ecole Polytechnique Fédérale de Lausanne
                 1015 Lausanne (Suisse)

  ===================================(Mars 1992)==============================*/

#ifndef __MDIM_H
#define __MDIM_H

#include <cassert>
#include <cstdlib>

/*-----------------------------------------------------------------------------

               ALLOCATION MEMOIRE POUR UN TABLEAU A UN INDICE

   DIM1 (p, longueur, type)

   Entrée
   ------
   p        : Nom du tableau à allouer
   longueur : Longueur de ce tableau
   type     : Type des composantes de ce tableau

   Note
   ----
   La macro alloue l'espace mémoire demande si celui-ci est disponible.
   Sinon, un message d'erreur est affiché. L'espace mémoire est rempli
   de zéros.

  -----------------------------------------------------------------------------*/

#define DIM1(p, longueur, type)                                                \
    {                                                                          \
        p = new type[longueur];                                                \
        if (p == (type*)NULL)                                                  \
        {                                                                      \
            exit(EXIT_FAILURE);                                                \
        }                                                                      \
        assert(longueur > 0);                                                  \
    }

/*-----------------------------------------------------------------------------

               ALLOCATION MEMOIRE POUR UN TABLEAU A DEUX INDICES

   DIM2 (p, lig, col, type)

   Entree
   ------
   p      : Nom du tableau à allouer
   lig    : Nombre de lignes de ce tableau
   col    : Nombre de colonnes de ce tableau
   type   : Type des composantes de ce tableau

   Note
   ----
   La macro alloue l'espace mémoire demandé si celui-ci est disponible.
   Sinon, un message d'erreur est affiché. L'espace mémoire obtenu est
   rempli de zeros.

  -----------------------------------------------------------------------------*/

#define DIM2(p, lig, col, type)                                                \
    {                                                                          \
        register type* pdata;                                                  \
        int iiiiiil;                                                           \
        DIM1(pdata, (lig) * (col), type);                                      \
        DIM1(p, lig, type*);                                                   \
        for (iiiiiil = 0; iiiiiil < lig; iiiiiil++)                            \
        {                                                                      \
            p[iiiiiil] = pdata;                                                \
            pdata += (col);                                                    \
        }                                                                      \
    }

#define D2FREE(p)                                                              \
    {                                                                          \
        delete[] * p;                                                          \
        delete[] p;                                                            \
    }

/*-----------------------------------------------------------------------------

               ALLOCATION MEMOIRE POUR UN TABLEAU A TROIS INDICES

   DIM3 (p, lig, col, couche, type)

   Entree
   ------
   p      : Nom du tableau à allouer
   lig    : Nombre de lignes de ce tableau
   col    : Nombre de colonnes de ce tableau
   couche : Nombre de couches de ce tableau
   type   : Type des composantes de ce tableau

   Note
   ----
   La macro alloue l'espace mémoire demandé si celui-ci est disponible.
   Sinon, un message d'erreur est affiché. L'espace mémoire obtenu est
   rempli de zeros.

   Ajoute par J.-L. Fattebert, juillet 94.

  -----------------------------------------------------------------------------*/

#define DIM3(p, lig, col, couche, type)                                        \
    {                                                                          \
        register type* pdata;                                                  \
        register type** ppoint;                                                \
        int iiiiiil;                                                           \
        int jjjjjjl;                                                           \
        DIM1(pdata, (lig) * (col) * (couche), type);                           \
        DIM1(ppoint, (lig) * (col), type*);                                    \
        DIM1(p, lig, type**);                                                  \
        for (iiiiiil = 0; iiiiiil < lig; iiiiiil++)                            \
        {                                                                      \
            for (jjjjjjl = 0; jjjjjjl < (col); jjjjjjl++)                      \
            {                                                                  \
                ppoint[jjjjjjl] = pdata;                                       \
                pdata += (couche);                                             \
            }                                                                  \
            p[iiiiiil] = ppoint;                                               \
            ppoint += (col);                                                   \
        }                                                                      \
    }

#define D3FREE(p)                                                              \
    {                                                                          \
        delete[](**p);                                                         \
        delete[](*p);                                                          \
        delete[](p);                                                           \
    }

#define DIM4(p, dim4, dim3, dim2, dim1, type)                                  \
    {                                                                          \
        register type* pdata;                                                  \
        register type** ppoint;                                                \
        register type*** pcouche;                                              \
        int iiiiiil;                                                           \
        int jjjjjjl;                                                           \
        int kkkkkkl;                                                           \
        DIM1(pdata, (dim4) * (dim3) * (dim2) * (dim1), type);                  \
        DIM1(ppoint, (dim4) * (dim3) * (dim2), type*);                         \
        DIM1(pcouche, (dim4) * (dim3), type**);                                \
        DIM1(p, (dim4), type***);                                              \
        for (iiiiiil = 0; iiiiiil < dim4; iiiiiil++)                           \
        {                                                                      \
            for (jjjjjjl = 0; jjjjjjl < dim3; jjjjjjl++)                       \
            {                                                                  \
                for (kkkkkkl = 0; kkkkkkl < dim2; kkkkkkl++)                   \
                {                                                              \
                    ppoint[kkkkkkl] = pdata;                                   \
                    pdata += dim1;                                             \
                }                                                              \
                pcouche[jjjjjjl] = ppoint;                                     \
                ppoint += dim2;                                                \
            }                                                                  \
            p[iiiiiil] = pcouche;                                              \
            pcouche += dim3;                                                   \
        }                                                                      \
    }

#define D4FREE(p)                                                              \
    {                                                                          \
        delete[](***p);                                                        \
        delete[](**p);                                                         \
        delete[](*p);                                                          \
        delete[](p);                                                           \
    }

#define DIM5(p, dim5, dim4, dim3, dim2, dim1, type)                            \
    {                                                                          \
        register type* pdata;                                                  \
        register type** ppoint;                                                \
        register type*** pcouche;                                              \
        register type**** p2;                                                  \
        int iiiiiil;                                                           \
        int jjjjjjl;                                                           \
        int kkkkkkl;                                                           \
        int lllllll;                                                           \
        DIM1(pdata, (dim5) * (dim4) * (dim3) * (dim2) * (dim1), type);         \
        DIM1(ppoint, (dim5) * (dim4) * (dim3) * (dim2), type*);                \
        DIM1(pcouche, (dim5) * (dim4) * (dim3), type**);                       \
        DIM1(p2, (dim5) * (dim4), type***);                                    \
        DIM1(p, (dim5), type****);                                             \
        for (iiiiiil = 0; iiiiiil < dim5; iiiiiil++)                           \
        {                                                                      \
            for (jjjjjjl = 0; jjjjjjl < dim4; jjjjjjl++)                       \
            {                                                                  \
                for (kkkkkkl = 0; kkkkkkl < dim3; kkkkkkl++)                   \
                {                                                              \
                    for (lllllll = 0; lllllll < dim2; lllllll++)               \
                    {                                                          \
                        ppoint[lllllll] = pdata;                               \
                        pdata += dim1;                                         \
                    }                                                          \
                    pcouche[kkkkkkl] = ppoint;                                 \
                    ppoint += dim2;                                            \
                }                                                              \
                p2[jjjjjjl] = pcouche;                                         \
                pcouche += dim3;                                               \
            }                                                                  \
            p[iiiiiil] = p2;                                                   \
            p2 += dim4;                                                        \
        }                                                                      \
    }

#define D5FREE(p)                                                              \
    {                                                                          \
        delete[](****p);                                                       \
        delete[](***p);                                                        \
        delete[](**p);                                                         \
        delete[](*p);                                                          \
        delete[](p);                                                           \
    }

#endif
